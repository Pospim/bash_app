#!/usr/bin/env python

"""
ELM Motif Search Script with Logging
====================================
This script searches for linear motifs in a protein sequence using ELM.

Pipeline:
1) Submits an amino acid sequence to ELM (via GET request).
2) Handles rate limits and retries with backoff.
3) Parses the returned TSV data for motif hits.
4) Filters results and saves to a user-defined output file.
"""
from pathlib import Path
import argparse
import json
import time
import requests
import logging
import os
import sys
import pandas as pd
from check_input import get_seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"

ELM_API_URL = "http://elm.eu.org/start_search/"
MAX_SEQUENCE_LENGTH = 2000     # ELM GET endpoint limit
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(BASE_DIR)
ELM_GO_FILE = f"{SCRIPT_DIR}/meta/elm_go_terms.tsv"

# high confidence classes
ALLOWED_CLASSES = {
    "LIG",  # Ligand binding motifs
    "DOC",  # Docking motifs
    "MOD",  # Post-translational modifications
    "TRG"   # Targeting signals
}

def load_go_mapping(go_terms_file=ELM_GO_FILE):
    """
    Loads ELM -> GO mappings from a TSV file.

    Returns:
        dict: { ELM_class: [GO:terms, evidence_score==1] } # ELM GO are not confident
    """
    #logging.info(f"Loading GO terms from: {go_terms_file}")
    go_path = Path(go_terms_file)

    if go_path.is_file():
        df = pd.read_csv(go_terms_file, sep='\t')

        go_mapping = df.groupby('ELM')['GOTerm'].apply(list).to_dict()
        return go_mapping

    else:
        logging.ERROR(f"Provided go annoation file is invalid. {go_terms_file}")
        return []


def submit_elm(sequence, max_attempts=5, delay=60):
    """
    Submits a protein sequence to ELM motif search.
    Retries on 429 errors or connection failures.
    """
    url = f"{ELM_API_URL}{sequence}"
    headers = {"Accept": "text/tab-separated-values"}

    for attempt in range(1, max_attempts+1):
        try:
            logging.info(f"Submitting ELM search (Attempt {attempt}/{max_attempts})...")
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                logging.info("ELM search completed successfully.")
                return response.text

            elif response.status_code == 429:
                logging.warning("Rate limit hit (429 Too Many Requests). Waiting 60 seconds before retrying...")
                time.sleep(delay)

            else:
                logging.error(f"Unexpected response: {response.status_code}")
                response.raise_for_status()

        except requests.RequestException as e:
            logging.error(f"Connection error on attempt {attempt}: {e}")

            if attempt < max_attempts:
                logging.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                raise RuntimeError("ELM search failed after maximum attempts.")

    raise RuntimeError("ELM search failed after all retry attempts.")

def parse_elm_tsv(tsv_text):
    """
    Parses the TSV output from ELM and filters motif hits.
    """
    lines = tsv_text.strip().split('\n')

    motifs = set()
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split('\t')

        motif = {
            "elm_identifier": cols[0],
            "start": int(cols[1]),
            "stop": int(cols[2]),
            "is_annotated": cols[3].lower() == "true",
            "is_phiblastmatch": cols[4].lower() == "true",
            "is_filtered": cols[5].lower() == "true",
            "structure": cols[6].lower() == "true",
            "topodomfilter": cols[7].lower() == "true",
            "taxonfilter": cols[8].lower() == "true"
        }

        if motif["is_filtered"]:
            continue  # Skip motifs filtered by ELM
        elm_id = motif["elm_identifier"]

        if elm_id.split('_')[0] in ALLOWED_CLASSES:
            motifs.add(elm_id)

    return motifs

def get_elm_to_go(motifs: set, go_mapping: dict):

    elm_to_go = {}

    for elm_id in motifs:
        go_terms = go_mapping.get(elm_id, [])
        if go_terms:
            terms_with_score = []

            for term in go_terms:
                term_with_score = (term, 1)
                terms_with_score.append(term_with_score)

            elm_to_go[elm_id] = terms_with_score

    logging.info(f"Mapped {len(elm_to_go)} ELM motifs to GO terms.")
    return  elm_to_go

def save_go_mapping(elm_to_go, output_json, output_csv):
    try:
        with open(output_json, "w") as f:
            json.dump(elm_to_go, f, indent=4)

    except Exception as E:
        logging.error(f"Failed to save GO terms to json: {e}")
        sys.exit(1)

    try:
        with open(output_csv, "w") as f:
            for _, value in elm_to_go.items():
                for val in value:
                    f.write(f"{val[0]}\t{val[1]}\n")
    except Exception as e:
        logging.error(f"Failed to save GO terms to csv: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run ELM motif search on a FASTA sequence.")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file.")
    parser.add_argument("--output_json", default="go_terms.json", help="Output file for ELM:GO terms dict.")
    parser.add_argument("--output_csv", default="go_terms.csv", help="Output file for GO terms - evidence score.")
    parser.add_argument("--go_map", default=ELM_GO_FILE, help="Path to elm to GO mapping")
    parser.add_argument("--max_attempts", type=int, default=5, help="Max attempts for ELM submission (default=5).")
    parser.add_argument("--delay", type=int, default=60, help="Delay in seconds between retries (default=60).")
    args = parser.parse_args()
    timestamp = int(time.time())

    try:
        sequence = get_seq(args.fasta, max_len=MAX_SEQUENCE_LENGTH)

    except Exception as e:
        logging.error(f"Error reading FASTA file: {e}")
        sys.exit(1)

    try:
        tsv_res = submit_elm(
            sequence=sequence,
            max_attempts=args.max_attempts,
            delay=args.delay
        )
        go_map = load_go_mapping(args.go_map)

        if not go_map:
            return []
    except Exception as e:
        logging.error(f"Failed to query ELM database: {e}")
        return []

    try:
        motifs = parse_elm_tsv(tsv_text=tsv_res)
        if not motifs:
            logging.info(f"No ELM results found")
            return []

        elm_to_go = get_elm_to_go(motifs, go_map)

        if not (elm_to_go and all(elm_to_go.values())):
            logging.info(f"No ELM results found")
            return []
        print(f"ELM_TO_GO\n{elm_to_go}")

    except Exception as e:
        logging.error(f"Failed to parse query results: {e}")
        return []

    try:
        save_go_mapping(elm_to_go=elm_to_go,
                        output_json=args.output_json,
                        output_csv=args.output_csv
        )
    except Exception as e:
        logging.error(f"Failed to save ELM results : {e}")
        sys.exit(1)

    elapsed = int(time.time() - timestamp)
    logging.info(f"Elapsed time: {format_time(elapsed)}")

if __name__ == "__main__":
    main()

# TOO RESTRICTIVE
# if not motif["is_annotated"] and not motif["is_phiblastmatch"]:
#    continue  # Keep only high-confidence motifs

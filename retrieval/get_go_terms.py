import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
import json
import argparse
import logging
import sys
import os
import pandas as pd
from map_to_uniprot import process_id_file

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(BASE_DIR)
BASE_URL = "http://rest.uniprot.org/uniprot/"
EVIDENCE_MAPPING = f"{SCRIPT_DIR}/meta/go_evidence_map.csv"

EVIDENCE_SCORES = {
    "IDA": 5,
    "EXP": 5,
    "IGI": 5,
    "IMP": 5,
    "IPI": 5,
    "IC":  4,
    "TAS": 4,
    "IBA": 4,
    "ISO": 3,
    "IEP": 3,
    "RCA": 3,
    "IGC": 3,
    "ISS": 2,
    "IEA": 2,
    "NAS": 2,
    "NR":  1,
    "ND":  0
}

def get_go_terms(uniprot_id: str, evidence_map: dict) -> list[tuple]:

    xml_url = f'{BASE_URL}{uniprot_id}.xml'
    try:
        with urllib.request.urlopen(xml_url) as response:
            xml = response.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        logging.error(f"Failed to retrieve GO data for {uniprot_id}: {e}")
        return []

    try:
        root = ET.fromstring(xml)
        namespace = "{http://uniprot.org/uniprot}"

        go_entries = []
        for db_ref in root.findall(f".//{namespace}dbReference[@type='GO']"):
            go_id = db_ref.get("id", "")
            evidence_code = None

            for prop in db_ref.findall(f"{namespace}property"):
                if prop.attrib["type"] == "evidence":
                    evidence_code = prop.attrib["value"]
                #elif prop.attrib["type"] == "term":
                #    term = prop.attrib["value"]
            evidence = evidence_map.get(evidence_code, "IEA")
            evidence_score = EVIDENCE_SCORES.get(evidence, 0)
            go_entries.append((
                go_id,
                #term,   # 'P' (Process), 'F' (Function), 'C' (Component)
                evidence_score
            ))
        return go_entries

    except Exception as e:
        logging.error(f"Failed to parse GO terms for {uniprot_id}: {e}")
        return []

def get_go_terms_batch(uniprot_ids: list, evidence_map: str = EVIDENCE_MAPPING) -> dict:

    map_df = pd.read_csv(evidence_map, sep="\t")
    map_dict = map_df.set_index("ECO_map")["evidence"].to_dict()
    uniprot_ids = list(set(uniprot_ids))
    go_terms = {}

    for uniprot_id in uniprot_ids:
        go_terms[uniprot_id] = get_go_terms(uniprot_id, map_dict)
    return go_terms

def validate_uniprot_ids(uniprot_ids: list) -> dict:
    """
    Validate UniProt IDs by checking if they return GO terms.
    :param uniprot_ids: List of UniProt IDs.
    :return: Dictionary of valid IDs and their GO terms.
    """
    valid_ids = {}
    for uniprot_id in uniprot_ids:
        try:
            go_terms = get_go_terms(uniprot_id)
            if go_terms:
                valid_ids[uniprot_id] = go_terms
            else:
                print(f"UniProt ID {uniprot_id} has no GO terms.")
        except RuntimeError as e:
            print(f"Error validating UniProt ID {uniprot_id}: {e}")
    return valid_ids

def load_uniprot_ids(file_path: str) -> list:
    """
    Load UniProt IDs from a file.

    :param file_path: Path to the file containing UniProt IDs.
    :return: List of UniProt IDs.
    """
    with open (file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

def save_go_terms(go_terms: dict, output_json: str, output_csv: str):
    """
    Save GO terms to a JSON file.

    :param go_terms: Dictionary of GO terms.
    :param output_file: Output file path.
    """
    try:
        with open(output_json, "w") as f:
            json.dump(go_terms, f, indent=4)

    except Exception as E:
        logging.error(f"Failed to save GO terms to json: {e}")
        sys.exit(1)

    try:
        with open(output_csv, "w") as f:
            for _, value in go_terms.items():
                for val in value:
                    f.write(f"{val[0]}\t{val[1]}\n")

    except Exception as e:
        logging.error(f"Failed to save GO terms to csv: {e}")
        sys.exit(1)

def load_go_terms_dict(file_path: str) -> dict:
    with open(file_path, "r") as f:
        return json.load(f)

def main(input_file: str, output_json: str, output_csv: str):

    uniprot_ids = process_id_file(input_file)
    if not uniprot_ids:
        sys.exit(1)

    go_terms = get_go_terms_batch(uniprot_ids)

    if go_terms:
        save_go_terms(go_terms, output_json, output_csv)
        #logging.info(f"GO terms saved to '{output_file}'.")
    else:
        logging.warning("No GO terms found for the provided UniProt IDs.")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve GO terms for UniProt IDs.")
    parser.add_argument("--file", required=True, help="Path to the file containing UniProt IDs.")
    parser.add_argument("--output_json", default="go_terms.json", help="Output file for Uniprot:GO terms dict.")
    parser.add_argument("--output_csv", default="go_terms.csv", help="Output file for GO terms - evidence score.")
    args = parser.parse_args()

    try:
        main(args.file, args.output_json, args.output_csv)
    except Exception as e:
        logging.error(f"Error during GO term retrieval: {e}")
        sys.exit(1)

#!/usr/bin/env python

"""
Foldseek Script with Logging
============================
This script aligns sequences or structures using FoldSeek.
It accepts a FASTA sequence string (already validated elsewhere) and:
1) Converts the sequence to PDB using ESMFold.
2) Submits the resulting PDB to FoldSeek for searching.
3) Polls for job completion.
4) Downloads and parses results.
5) Filters results by user-defined thresholds.
"""
import re
import argparse
import time
import tarfile
import requests
import pandas as pd
import logging
import sys
from io import BytesIO
from check_input import get_seq


"""
For pdb100:

Use SIFTS (Structure Integration with Function, Taxonomy, and Sequences) or the PDBe API to map PDB IDs to UniProt IDs.
For afdb50 and afdb-proteome:

Check if IDs are UniProt-like (AF-P12345-F1) or require mapping using the AlphaFold DB.
For afdb-swissprot:

Direct use of the UniProt IDs is safe.
"""

DBS = ['afdb50', 'afdb-swissprot', 'afdb-proteome']
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"

def get_esmfold_pdb(sequence_str: str) -> bytes:
    """
    Converts a sequence to PDB format using ESMFold API.

    Args:
        sequence_str (str): Amino acid sequence.
    Returns:
        bytes: PDB file content returned by ESMFold.
    """
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    logging.info("Converting sequence to PDB using ESMFold...")
    resp = requests.post(url, headers=headers, data=sequence_str, verify=True)
    resp.raise_for_status()
    return resp.content

def submit_ticket(file_data: bytes, file_name: str, db: str = "afdb50") -> dict:
    """
    Submit a job to FoldSeek.
    Args:
        file_data (bytes): File content.
        file_name (str): File name.
        db (str): Database to use.

    Returns:
        dict: FoldSeek job submission response.
    """
    payload = [("q", (file_name, file_data, "application/octet-stream")),
           ("mode", (None, "3diaa")),
           ("database[]", (None, db))]

    url = "https://search.foldseek.com/api/ticket"
    logging.info("Submitting job to FoldSeek...")
    resp = requests.post(url, files=payload)
    resp.raise_for_status()
    return resp.json()

def poll_job(job_id: str, interval: int = 5, max_attempts: int = 100) -> dict:
    """
    Poll FoldSeek job status until completion.
    Args:
        job_id (str): Job ID.
        interval (int): Polling interval in seconds.
        max_attempts (int): Maximum polling attempts.

    Returns:
        dict: Job status data.

    Raises:
        TimeoutError: If the job does not complete within the max attempts.
    """
    url = f"https://search.foldseek.com/api/ticket/{job_id}"
    for attempt in range(1, max_attempts+1):
        resp = requests.get(url)
        resp.raise_for_status()
        data = resp.json()
        status = data.get("status", "")

        if status == "COMPLETE":
            logging.info("FoldSeek job completed.")
            return data
        elif status == "ERROR":
            logging.error(f"FoldSeek job failed: {data}")
            raise RuntimeError(f"Job error: {data}")
        elif status == "RATELIMIT":
            logging.warning("Rate limit encountered. Waiting 10 seconds...")
            time.sleep(10)
        else:
            logging.info(f"Attempt {attempt}: status={status}. Retrying in {interval} seconds...")
            time.sleep(interval)
    raise TimeoutError("Job did not complete within the maximum attempts.")

def download_results(job_id: str) -> bytes:
    """
    Download results from a FoldSeek job.
    Args:
        job_id (str): Job ID.
    Returns:
        bytes: Compressed results file.
    """
    url = f"https://search.foldseek.com/api/result/download/{job_id}"
    logging.info("Downloading FoldSeek results...")
    resp = requests.get(url, stream=True)
    resp.raise_for_status()
    return resp.content

def parse_results(tar_bytes: bytes) -> pd.DataFrame:
    """
    Extract FoldSeek results into a DataFrame.
    Args:
        tar_bytes (bytes): Compressed results file.

    Returns:
        pd.DataFrame: Parsed results as a DataFrame.
    """
    frames = []
    with tarfile.open(fileobj=BytesIO(tar_bytes), mode='r:gz') as tar:
        for member in tar:
            if member.isreg():
                extracted = tar.extractfile(member)
                if extracted is not None:
                    frames.append(pd.read_csv(extracted, sep='\t', header=None))

    if not frames:
        logging.warning("No data found in FoldSeek results.")
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)

def get_fold(pdb_data: str, pdb_file: str="converted_by_esmfold.pdb", db: str = 'afdb50',
            max_eval: float = 10, min_ident: float = 30.0, k_max: int = 10) -> list:
    """
    Execute the full FoldSeek pipeline, starting from a FASTA string.

    Steps:
      1. Submit to FoldSeek.
      2. Poll until complete.
      3. Download + parse results.
      4. Filter hits by e-value and identity, then return top k_max IDs.
    Args:
        pdb_data (str): Fasta sequence converted to PDB
        pdb_file (str): Converted PDB file
        sequence_str (str): FASTA sequence (already checked).
        db (str, optional): Database to query (default: "afdb50").
        max_eval (float, optional): Maximum e-value threshold for filtering (default: 10).
        min_identity (float, optional): Minimum identity threshold (default: 30.0).
        k_max (int, optional): Maximum number of results to retain (default: 10).

    Returns:
        list: A list of UniProt identifiers (strings) extracted from top hits.
    """

    logging.info(f"Running FoldSeek on database: {db}...")
    job_resp = submit_ticket(pdb_data, pdb_file, db)
    job_id = job_resp.get("id")
    logging.info(f"FoldSeek job ID: {job_id}")

    poll_job(job_id)
    results = download_results(job_id)
    df = parse_results(results)

    filtered = df[(df[10] <= max_eval) & (df[2] >= min_ident)]
    filtered = filtered.sort_values(by=[10, 2], ascending=[True, False]).head(k_max)

    fs_results = filtered[1].tolist()

    uniprot_ids = [re.search(r"AF-([A-Z0-9]+)-", entry).group(1) for entry in fs_results if "AF-" in entry]
    logging.info(f"Filtered results: {len(uniprot_ids)} hits remain after filtering.")
    return uniprot_ids

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit a file to FoldSeek for alignment.")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file containing the query sequence.")
    parser.add_argument("--dbs",nargs='+',default=["afdb50"],help="One or more FoldSeek databases to query (default: afdb50).")
    parser.add_argument("--max_eval", type=float, default=10, help="Maximum e-value for filtering (default: 10).")
    parser.add_argument("--min_ident", type=float, default=60.0, help="Minimum identity percentage (default: 60).")
    parser.add_argument("--kMax", type=int, default=10, help="Maximum hits to retain (default: 10).")
    parser.add_argument("--max_len", type=int, default=400, help="Maximum sequence length (default: 400).")
    parser.add_argument("--output_file", required=True,  help="Path to file where results will be saved.")
    args = parser.parse_args()

    try:
        sequence = get_seq(args.fasta, max_len=args.max_len)

    except Exception as e:
        logging.error(f"Error reading FASTA file: {e}")
        sys.exit(1)

    valid_dbs = [db for db in args.dbs if db in DBS]

    if not valid_dbs:
        logging.warning(f"Invalid databases. Available options are: {', '.join(DBS)}.")
        logging.info("Using default database 'afdb50'.")
        valid_dbs = ["afdb50"]

    pdb_data = get_esmfold_pdb(sequence)
    pdb_file = "converted_by_esmfold.pdb"

    all_results = []
    timestamp = int(time.time())

    for db in valid_dbs:
        try:
            results = get_fold(
                pdb_data=pdb_data,
                pdb_file=pdb_file,
                db=db,
                max_eval=args.max_eval,
                min_ident=args.min_ident,
                k_max=args.kMax
            )
            #if results: print(results)
            all_results.extend(results)
        except Exception as e:
            logging.error(f"Error during FoldSeek pipeline on DB '{db}': {e}")

    elapsed = int(time.time()) - timestamp
    logging.info(f"Elapsed: {format_time(elapsed)}")

    all_results = list(set(all_results))

    try:
        with open(args.output_file, 'w') as outf:
            for result in all_results:
                outf.write(result.strip() + '\n')

    except IOError as e:
        logging.error(f"Error writing results to {args.output_file}: {e}")
        sys.exit(1)

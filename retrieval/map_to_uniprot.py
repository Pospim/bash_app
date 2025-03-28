#!/usr/bin/env python
import sys
import requests
import time
import re
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

UNIPROT_BASE_URL = "https://rest.uniprot.org/idmapping"
SUBMIT_URL = f"{UNIPROT_BASE_URL}/run"
STATUS_URL = f"{UNIPROT_BASE_URL}/status"
RESULTS_URL = f"{UNIPROT_BASE_URL}/results"

ID_TYPE_MAPPING = {
    "UniProtKB": r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^A0A[A-Z0-9]{7}$",
    "RefSeq_Protein": r"^(XP|WP|YP|NP|ZP)_[0-9]+(\.[0-9]+)?$",
    "PDB": r"^[0-9][A-Za-z0-9]{3}(_[A-Za-z0-9])?$",
    #"EMBL-GenBank-DDBJ": r"^[A-Z]{1,3}[0-9]{5,6}(\.[0-9]+)?$"
}

def clean_pdb(pdb_id):
    """Remove chain info and ensure PDB ID is uppercased."""
    return pdb_id.split("_")[0].upper()

def classify_id(id_str):
    for db_type, pattern in ID_TYPE_MAPPING.items():
        if re.match(pattern, id_str):
            return db_type
    return None

def extract_uniprot(status):
    uniprot_mapping = {}

    for result in status.get("results", []):
        original = result.get("from")
        mapped = result.get("to", {}).get("primaryAccession")
        uniprot_mapping[original] = mapped if mapped else None

    for failed in status.get("failedIds", []):
        uniprot_mapping[failed] = None

    return uniprot_mapping

def submit_mapping(ids, from_db, to_db="UniProtKB"):
    #logging.info(f"Submitting mapping request from {from_db} to {to_db} for {len(ids)} IDs.")
    response = requests.post(SUBMIT_URL, data={
        "from": from_db,
        "to": to_db,
        "ids": ",".join(ids)
    })
    if response.status_code != 200:
        logging.error(f"Failed to submit mapping job. Status code: {response.text}")
        return None
    return response.json()["jobId"]

def check_mapping_status(job_id, max_retries=3, poll_interval=5, max_wait_time=200):
    """
    Check the status of the mapping job with retry logic and exponential backoff.

    Args:
        job_id (str): The job ID returned by the mapping submission.
        max_retries (int): Maximum number of retries on failure.
        poll_interval (int): Base interval (in seconds) between status checks.
        max_wait_time (int): Maximum time (in seconds) to wait for job completion.

    Returns:
        dict or None: Mapping status if successful, None otherwise.
    """
    start_time = time.time()
    retry_cnt = 0

    while True:
        try:
            response = requests.get(f'{STATUS_URL}/{job_id}')
            if response.status_code == 200:

                status = response.json()

                if 'results' in status and status['results']:
                    #logging.info("Mapping job completed")
                    return status

                elapsed = time.time() - start_time

                if elapsed >= max_wait_time:
                    logging.error(f"Mapping job timed out after waiting {max_wait_time} seconds.")
                    return None

                time.sleep(poll_interval)
                continue

            else:
                logging.warning(f"Failed to check job status: {response.status_code}")
                retry_cnt += 1

                if retry_cnt > max_retries:
                    logging.error(f"Exceeded maximum retries ({max_retries}) for job status check.")
                    return None

                logging.info(f"Retrying in {poll_interval} (attempt {retry_cnt}/{max_retries})")
                time.sleep(poll_interval)

        except requests.exceptions.RequestException as e:
            retry_cnt += 1
            logging.error(f"Network error occurred: {str(e)}")

            if retry_cnt > max_retries:
                    logging.error(f"Exceeded maximum retries ({max_retries}) for job status check.")
                    return None
            logging.info(f"Retrying in {poll_interval} (attempt {retry_cnt}/{max_retries})")
            time.sleep(poll_interval)

def get_mapping_results(job_id):
    response = requests.get(f'{RESULTS_URL}/{job_id}')
    if response.status_code != 200:
        logging.error("Failed to get results: {response.text}")
        return None
    return response.json()

def map_ids_to_uniprot(mixed_ids):
    classified_ids =  {}
    final_mapping = {}

    # 1️) Classify IDs by their type (RefSeq, PDB, UniProt, etc.)
    for id in mixed_ids:
        id_clean = id.strip()
        id_type = classify_id(id_clean)
        if id_type == "PDB":
            id_clean = clean_pdb(id_clean)

        if not id_type:
            logging.warning(f"Could not classify ID: {id_clean}")
            final_mapping[id_clean] = None
        else:
            classified_ids.setdefault(id_type, []).append(id_clean)

    # 2) Process each ID type
    for source_db, ids in classified_ids.items():
        if source_db == "UniProtKB":
            # No mapping needed
            for uid in ids:
                final_mapping[uid] = uid
            continue

        logging.info(f"Mapping {len(ids)} IDs from {source_db} to UniProtKB...")
        job_id = submit_mapping(list(set(ids)), source_db, "UniProtKB")

        if job_id:
            status = check_mapping_status(job_id)
            if status:
                uniprot_mapping = extract_uniprot(status)
                final_mapping.update(uniprot_mapping)
        else:
            logging.error(f"Mapping failed for {source_db} ID")
            continue

    return final_mapping

def process_id_file(input_file: str) -> list:
    logging.info("Mapping protein IDs to UniProt IDs...")

    try:
        with open(input_file, 'r') as f:
            mixed_ids = list(set([line.strip() for line in f if line.strip()]))

        if not mixed_ids:
            logging.error("Protein IDs file is empty.")
            sys.exit(1)

        mapping = map_ids_to_uniprot(mixed_ids)
        uniprot_ids = sorted(set(uid for uid in mapping.values() if uid))

        with open(input_file, 'w') as f:
            for uid in uniprot_ids:
                f.write(f"{uid}\n")

        if uniprot_ids:
            logging.info(f"{len(uniprot_ids)} IDs successfully mapped to UniProt IDs")
            #print(uniprot_ids)
            return uniprot_ids

        else:
            logging.warning(f"0 IDs mapped to UniProt IDs, exiting...")
            sys.exit(1)

    except Exception as e:
        logging.error(f'Error mapping protein IDs to UniProt IDs: {e}')
        return []

if __name__ == "__main__":
    input_file = "protein_ids.txt"
    process_id_file(input_file)
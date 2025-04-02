import requests
import time
import argparse
import sys
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

NAMESPACE = {
    '1': 'Biological process',
    '2': 'Cellular component',
    '3': 'Molecular function'
}

def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"

BASE_URL = "http://revigo.irb.hr/"

def load_csv(file_path):
    path = Path(file_path)
    if not path.is_file():
        logging.error(f"File does not exist: {file_path}")
        sys.exit(1)

    try:
        terms_with_score = pd.read_csv(file_path, sep='\t')
        if terms_with_score.empty:
            return []

        if terms_with_score.shape[1] != 2:
            logging.error(f"GO - score file must have two columns: {file_path}\n{e}")
            sys.exit(1)

    except Exception as e:
        logging.error(f"Failed to load GO - score file: {file_path}\n{e}")
        sys.exit(1)

    return terms_with_score

def submit_revigo(terms_with_score: pd.DataFrame, cutoff:str ='0.7',
                  val_type:str ='Higher', species:str ='0',
                  measure:str ='SIMREL', max_attempts:int =5, delay:int =60):
    url = f"{BASE_URL}StartJob"
    if terms_with_score.empty:
        logging.error(f"GO - score file failed to load")
        sys.exit(1)

    payload = {'cutoff': cutoff,
               'valueType': val_type,
               'measure': measure,
               'speciesTaxon': species,
               'goList': terms_with_score
    }
    for attempt in range(1, max_attempts+1):
        try:
            logging.info(f"Submitting REVIGO job (Attempt {attempt}/{max_attempts})...")
            r = requests.post(url, data=payload)
            r.raise_for_status()

            job_id = r.json().get('jobid')

            if job_id:
                logging.info(f"REVIGO job submitted successfully. Job ID: {job_id}")
                return job_id

        except requests.RequestException as e:
            logging.error(f"POST request failed on attempt {attempt}: {e}")

            if attempt < max_attempts:
                logging.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                raise RuntimeError("ELM search failed after maximum attempts.")
    else:
        raise RuntimeError("ELM search failed after all retry attempts.")

def get_results(job_id: str):
    url = f"{BASE_URL}QueryJob"
    params = {'jobid': job_id, 'type': 'jstatus'}

    running = 1
    try:
        while running != 0:
            logging.info(f"Checking REVIGO job status...")

            r = requests.get(url=url, params=params)
            r.raise_for_status()
            response = r.json()
            running = response.get('running')

            if running != 0:
                logging.info(f"Progress: {response.get('message')}")
            else:
                logging.info(f"Job completed with status: {response}")

            time.sleep(1)

        return response

    except Exception as e:
        logging.error(f"Failed to get results from ReviGO: {e}")
        sys.exit(1)

def parse_results(job_id: str, namespaces: list[str], result_types: list[str]) -> list[dict]:
    url = f"{BASE_URL}QueryJob"
    results = []

    try:
        for ontology in namespaces:
            result = {}
            for type in result_types:
                logging.info(f"Fetching REVIGO results (ontology={NAMESPACE[ontology]}, type={type})...")
                params = {'jobid': job_id,
                          'namespace': ontology,
                          'type': type
                }
                r = requests.get(url=url, params=params)
                r.raise_for_status()
                result[type] = r.text

            results.append(result)
        return results

    except requests.RequestException as e:
        logging.error(f"Failed to fetch REVIGO results: {e}")
        sys.exit(1)

ontologies=['1', '2', '3']
result_types = ['Table', 'jScatterplot']
file = "/home/pospim/Desktop/work/GOLizard/bash_app/tmp/go_terms.csv"
terms = load_csv(file_path=file)

id = submit_revigo(terms_with_score=terms)
print(id)

ret_val = get_results(id)
print(ret_val)

results = parse_results(id, ontologies, result_types=result_types)
print(results)
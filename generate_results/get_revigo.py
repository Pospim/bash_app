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

BASE_DIR = Path(__file__).parent.parent
BASE_URL = "http://revigo.irb.hr/"
HEADERS = {"Content-Type": "application/x-www-form-urlencoded"}

# Maps numeric IDs to readable names:
NAMESPACE = {
    '1': 'biological_process',
    '2': 'cellular_component',
    '3': 'molecular_function'
}
RESULT_TYPES = ["jTable", "jScatterplot", "jCytoscape"]

def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"

def load_csv(file_path: str) -> str:
    """
    Reads a local TSV file with two columns: (GO_ID, Score).
    Returns a space-separated string, e.g.:
      'GO:0008150 0.0002\\nGO:0000003 0.0453\\n...'
    This string is suitable for the 'goList' parameter in REVIGO.
    """
    path = Path(file_path)
    if not path.is_file():
        logging.error(f"File does not exist: {file_path}")
        sys.exit(1)

    try:
        terms_with_score = pd.read_csv(file_path, sep='\t')
        if terms_with_score.empty:
            logging.error(f"GO list file is empty: {file_path}")
            return []

        if terms_with_score.shape[1] != 2:
            logging.error("GO - score file must have exactly two columns (GO Term, Score).")
            sys.exit(1)

    except Exception as e:
        logging.error(f"Failed to load GO - score file: {file_path}\n{e}")
        sys.exit(1)

    # Return a space-separated string for REVIGO
    return terms_with_score.to_csv(sep=" ", index=False, header=False)

def submit_revigo(terms_with_score: str,
                  cutoff:str = '0.7',
                  val_type:str = 'Higher',
                  species:str = '0',
                  measure:str = 'SIMREL',
                  max_attempts:int = 5,
                  delay:int = 60) -> str:
    """
    Submits a job to REVIGO with retry logic, returning the job ID upon success.
    """
    url = f"{BASE_URL}StartJob"
    if not terms_with_score.strip():
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
            r = requests.post(url, headers=HEADERS, data=payload, timeout=30)
            r.raise_for_status()

            job_data = r.json()
            job_id = job_data.get('jobid', -1)

            if job_id != -1:
                logging.info(f"REVIGO job submitted successfully. Job ID: {job_id}")
                return str(job_id)
            else:
                msg = job_data.get('message', "unknown error from REVIGO")
                raise RuntimeError(f"REVIGO error: {msg}")

        except requests.RequestException as e:
            logging.error(f"POST request failed on attempt {attempt}: {e}")

            if attempt < max_attempts:
                logging.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                raise RuntimeError("ELM search failed after maximum attempts.")

    raise RuntimeError("ELM search failed after all retry attempts.")

def wait_for_completion(job_id: str, output_dir:str, max_wait: int = 60):
    """
    Polls the jstatus endpoint until the job is complete or we hit max_wait seconds.
    """
    url = f"{BASE_URL}QueryJob"
    params = {'jobid': job_id, 'type': 'jstatus'}

    logging.info(f"Waiting for REVIGO job {job_id} to complete (timeout={max_wait}s)...")
    for _ in range(max_wait):
        try:
            r = requests.get(url, params=params, timeout=10)
            r.raise_for_status()
            status = r.json()

            if status.get('running') == 0:
                logging.info(f"Job completed with status: {status}")
                results_link = f"{BASE_URL}/Results?jobid={job_id}"
                out_path = Path(output_dir) / "revigo_link.txt"
                out_path.parent.mkdir(parents=True, exist_ok=True)
                with open(out_path, 'w') as f:
                    f.write(results_link)
                print("----------------------------------------------------------")
                print(f"RESULTS: {results_link}\n       !!!available for 15 minutes!!!")
                print("----------------------------------------------------------")

                return
            #else:
                #msg = status.get('message', "No progress message.")
                #logging.info(f"Progress: {msg}")
        except requests.RequestException as e:
            logging.error(f"Failed to query REVIGO job status: {e}")

        time.sleep(1)

    raise TimeoutError(f"REVIGO job {job_id} did not complete within {max_wait} seconds.")

def parse_results(job_id: str,
                  namespaces: list[str],
                  result_types: list[str],
                  output_dir: str):
    """
    Fetches the REVIGO results for each requested namespace & output type.
    Writes each result to {output_dir}/ontology/result_type.json
      {
        '1': {'jTable': "...", 'jScatterplot': "..."},
        '2': {...},
        ...
      }
    """
    url = f"{BASE_URL}QueryJob"

    time.sleep(5)

    for ns_id in namespaces:
        ontology = NAMESPACE.get(ns_id)
        for out_type in result_types:
            params = {'jobid': job_id,
                      'namespace': ns_id,
                      'type': out_type
            }
            try:
                r = requests.get(url=url, params=params, timeout=60)
                r.raise_for_status()

                if 'The Job has an errors, no data available' in r.text.lower():
                    logging.warning(f"REVIGO returned error for {ontology}/{out_type}: {r.text.strip()}")
                    continue
                else:
                    out_path = Path(output_dir) / ontology / f"{out_type}.json"
                    out_path.parent.mkdir(parents=True, exist_ok=True)

                    with open(out_path, 'w') as f:
                        f.write(r.text)
                    #logging.info(f"Saved {out_path}")
            except requests.RequestException as e:
                logging.error(f"Failed to fetch results for {ontology}/{out_type}: {e}")
                continue


def main():
    parser = argparse.ArgumentParser(description="Runs REVIGO for given GO : Score list.")
    parser.add_argument("--go_terms", required=True, help="Path to input GO_terms file.")
    parser.add_argument("--output_dir", default=f"{BASE_DIR}/results", help=f"Output dir for all Ontology sub-directories default={BASE_DIR}/results).")
    parser.add_argument("--result_type",nargs='+',default=RESULT_TYPES,help="One or more REVIGO result types.")
    parser.add_argument("--max_attempts", type=int, default=5, help="Max attempts for REVIGO submission (default=5).")
    parser.add_argument("--max_wait", type=int, default=120, help="Max wait time for REVIGO submission (default=120).")
    parser.add_argument("--cutoff", type=float, default=0.7, help="Cutoff for REVIGO algorithm (default=0.7).")
    #parser.add_argument("--delay", type=int, default=60, help="Delay in seconds between retries (default=60).")
    # TODO List of types + ontologies

    args = parser.parse_args()
    timestamp = int(time.time())

    ontologies = ontologies = ['1', '2', '3']
    result_types = [result for result in args.result_type if result in RESULT_TYPES]

    if not result_types:
        logging.warning(f"Invalid REVIGO result types. Available options are: {', '.join(RESULT_TYPES)}.")
        logging.info(f"Using default {', '.join(RESULT_TYPES)}.")
        result_types = RESULT_TYPES

    try:
        # 1) Load the local GO data
        terms_str = load_csv(file_path=args.go_terms)

        # 2) Submit job
        job_id = submit_revigo(terms_with_score=terms_str)

        # 3) Wait for completion
        wait_for_completion(job_id, output_dir=args.output_dir, max_wait=args.max_wait)

        # 4) Fetch results
        parse_results(job_id, ontologies, result_types, output_dir=args.output_dir)

        elapsed = int(time.time() - timestamp)
        logging.info(f"Elapsed time: {format_time(elapsed)}")
        print(f"[INFO] Saved REVIGO results => {args.output_dir}")

    except Exception as e:
        logging.error(f"Error getting REVIGO results")
        sys.exit(1)


if __name__ == "__main__":
    main()
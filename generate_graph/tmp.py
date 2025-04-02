import requests
import logging
import time
import sys

BASE_URL = "http://revigo.irb.hr/"

NAMESPACE_MAP = {
    "Biological process": "1",
    "Cellular component": "2",
    "Molecular function": "3"
}

OUTPUT_TYPES = ["jTable", "jScatterplot"]  # You can add more: jTreeMap, jClouds, etc.

HEADERS = {"Content-Type": "application/x-www-form-urlencoded"}


def submit_revigo_job(go_list: str,
                      cutoff: float = 0.7,
                      value_type: str = "PValue",
                      species: int = 0,
                      measure: str = "SIMREL",
                      remove_obsolete: bool = True) -> str:
    """Submits a job to REVIGO and returns the job ID."""
    url = f"{BASE_URL}StartJob"
    data = {
        "goList": go_list,
        "cutoff": cutoff,
        "valueType": value_type,
        "speciesTaxon": species,
        "measure": measure,
        "removeObsolete": str(remove_obsolete).lower()
    }

    logging.info("Submitting REVIGO job...")
    response = requests.post(url, headers=HEADERS, data=data)
    response.raise_for_status()

    job_data = response.json()
    if job_data["jobid"] == -1:
        raise RuntimeError(f"REVIGO error: {job_data['message']}")

    logging.info(f"REVIGO job submitted. Job ID: {job_data['jobid']}")
    return str(job_data["jobid"])


def wait_for_completion(job_id: str, max_wait: int = 60):
    """Polls until the job is completed or timeout occurs."""
    url = f"{BASE_URL}QueryJob"
    params = {"jobid": job_id, "type": "jstatus"}

    logging.info("Waiting for REVIGO job to complete...")
    for i in range(max_wait):
        resp = requests.get(url, params=params)
        resp.raise_for_status()
        status = resp.json()
        logging.info(f"Status: {status}")

        if status["running"] == 0:
            return
        time.sleep(1)
    raise TimeoutError(f"REVIGO job {job_id} timed out after {max_wait} seconds.")


def fetch_revigo_results(job_id: str,
                         namespaces: list[str],
                         result_types: list[str]) -> dict:
    """Fetches all requested result types for all GO namespaces."""
    results = {}
    for name in namespaces:
        ns_id = NAMESPACE_MAP[name]
        ns_results = {}
        for result_type in result_types:
            logging.info(f"Fetching {result_type} for {name}")
            url = f"{BASE_URL}QueryJob"
            params = {
                "jobid": job_id,
                "namespace": ns_id,
                "type": result_type
            }
            r = requests.get(url, params=params)
            r.raise_for_status()
            if 'error' in r.text.lower():
                logging.warning(f"Error for {name} / {result_type}: {r.text.strip()}")
                ns_results[result_type] = None
            else:
                ns_results[result_type] = r.text
        results[name] = ns_results
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    go_term_list = "GO:0008150 0.001\nGO:0009987 0.05\nGO:0050896 0.01"
    job_id = submit_revigo_job(go_term_list)
    wait_for_completion(job_id)

    results = fetch_revigo_results(
        job_id,
        namespaces=["Biological process", "Cellular component", "Molecular function"],
        result_types=["jTable", "jScatterplot"]
    )

    # Print or parse TSV / JSON outputs
    for ns, data in results.items():
        for t, content in data.items():
            print(f"\n===== {ns} / {t} =====\n")
            print(content[:1000])  # Print preview

import tempfile
import logging
import argparse
import time
import os
import subprocess
import shlex
import sys
from get_blast_api import format_time, validate_hit, extract_hits, cluster_sequences
# --------------------------------------------------
# Configuration
# --------------------------------------------------
PROGS = ['blastp', 'tblastn']

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def check_local_database(db_path: str) -> bool:
    """
    Checks if a local BLAST database is valid by running 'blastdbcmd -info'.

    Args:
        db_path (str): Path (excluding file extensions) to the local BLAST DB
                       (e.g., '/path/to/blastdb/swissprot').

    Returns:
        bool: True if the database is recognized and not corrupt, False otherwise.
    """
    db_dir = os.path.dirname(db_path)
    if db_dir and not os.path.isdir(db_dir):
        logging.error(f"Directory does not exist: {db_dir}")
        return False

    cmd =f'blastdbcmd -db {db_path} -info' # Build cmd command
    args = shlex.split(cmd)
    try:
        result = subprocess.run(args, check=True, capture_output=True, text=True)
        logging.debug(f'DB check output: {result.stdout.strip()}')
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to verify DB '{db_path}': {e.stderr.strip()}")
        return False

    except FileNotFoundError:
        logging.error("blastdbcmd is not installed or not in PATH.")
        return False
# --------------------------------------------------
# Local BLAST pipeline
# --------------------------------------------------
def run_blast_local(fasta_file: str,
                    db_path: str,
                    k_max: int = 10,
                    max_eval: float = 1e-5,
                    min_ident: float = 30.0,
                    cluster_thresh: float = 0.9,
                    prog: str = "blastp",
                    ) -> list:
    """
    Run a local BLAST search against a downloaded DB, parse hits, apply filters,
    cluster results, and return representative IDs.

    Args:
        fasta (str): Amino acid sequence to BLAST.
        db_path (str): Path to local BLAST database (excluding file extensions).
        k_max (int): Max hits to keep after filtering by rank.
        max_eval (float): Maximum allowed e-value.
        min_identity (float): Minimum percent identity.
        cluster_thresh (float): Clustering threshold (0-1).
        prog (str): BLAST program to run (e.g., 'blastp', 'blastn', etc.).
    Returns:
        list[str]: List of representative hit IDs after clustering.
    """
    out_xml = None
    if not check_local_database(db_path):
        return []
    try:
        # 1) Create a temporary file for BLAST output (XML)
        with tempfile.NamedTemporaryFile(suffix='.xml', delete=False) as out_fh:
            out_xml=out_fh.name

        # 2) Run local BLAST command
        logging.info(f"Running local BLAST on DB: {db_path} ...")
        cmd_args = [
            prog,
            '-query', fasta_file,
            '-db',db_path,
            '-evalue', str(max_eval*100), # Filtering later
            '-outfmt', '5', # XML
            '-out',out_xml
        ]
        result = subprocess.run(cmd_args, capture_output=True, text=True)
        if result.returncode != 0:
            logging.warning(f'Blast stderr: {result.stderr}')
            return []

        # 3) Parse the resulting XML
        with open(out_xml, 'r') as xml_handle:
            hits_df = extract_hits(xml_handle)
        if hits_df.empty:
            logging.warning("No hits found in local BLAST search.")
            logging.info("......................")
            return []

        # Running Blast locally might have different order
        hits_df = hits_df.sort_values("evalue", ascending=True).reset_index(drop=True)

        # 4) Filter hits using your validate_hit function
        hits_df = hits_df[hits_df.apply(
        lambda x: validate_hit(
            rank=x.name + 1,
            k_max=k_max,
            eval=x.evalue,
            max_eval=max_eval,
            identity=x.identity,
            min_ident=min_ident),
        axis=1
        )]
        if hits_df.empty:
            logging.warning("No hits passed the filtering criteria.")
            logging.info("......................")
            return []
        logging.info(f"{len(hits_df)} hits passed filtering. Proceeding to clustering.")

        # 5) Cluster hits
        clustered_df = cluster_sequences(hits_df, threshold=cluster_thresh)
        logging.info(f"Clustering reduced the hits to {len(clustered_df)} unique representatives.")

        return clustered_df['id'].tolist()

    except Exception as e:

        logging.error(f"Error in run_local_blast_tempfile: {e}")
        return []

    finally:
        if out_xml and os.path.exists(out_xml):
            os.remove(out_xml)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST localy to get a list of similar genes.")
    parser.add_argument("--fasta", required=True, help="Amino acid sequence for BLAST query.")
    parser.add_argument("--db_paths",nargs='+',required=True,help="One or more local BLAST databases (no file extension).")
    parser.add_argument("--prog", default="blastp", help="BLAST program (default: blastp).")
    parser.add_argument("--kMax", type=int, default=10, help="Maximum rank of hits to consider (default: 10).")
    parser.add_argument("--max_eval", type=float, default=1e-5, help="Maximum allowed e-value (default: 1e-5).")
    parser.add_argument("--min_ident", type=float, default=40.0, help="Minimum identity percentage (default: 40.0).")
    parser.add_argument("--cluster", type=float, default=0.9, help="Clustering threshold (default: 0.9).")
    parser.add_argument("--output_file", required=True, help="Path to file where results will be saved" )
    args = parser.parse_args()

    # Validate the selected program
    if args.prog not in PROGS:
        logging.error(f"Invalid program '{args.prog}'. Available options are: {', '.join(PROGS)}.")
        logging.info("Using default program 'blastp'.")
        args.prog = "blastp"

    timestamp = int(time.time())

    all_results = []
    for db in args.db_paths:
        try:
            results = run_blast_local(
                fasta_file=args.fasta,
                db_path=db,
                k_max=args.kMax,
                max_eval=args.max_eval,
                min_ident=args.min_ident,
                cluster_thresh=args.cluster,
                prog=args.prog,
            )
            #if results: print(results)
            all_results.extend(results)
        except Exception as e:
            logging.error(f"Error during BLAST pipeline on DB '{db}': {e}")

    elapsed = int(time.time()) - timestamp
    logging.info(f"Elapsed: {format_time(elapsed)}")

    all_results = list(set(all_results))
    try:
        with open(args.output_file, 'w') as outf:
            for result in all_results:
                outf.write(result.strip() + '\n')
        logging.info(f"Blast local results saved => {args.output_file}")
    except IOError as e:
        logging.error(f"Error writing results to {args.output_file}: {e}")
        sys.exit(1)

#!/usr/bin/env python

import numpy as np
import pandas as pd
import edlib
import logging
import argparse
import sys
import time
import re
import os
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from check_input import get_seq
from pathlib import Path

# --------------------------------------------------
# Configuration
# --------------------------------------------------
DBS = ['swissprot', 'refseq_protein', 'nr', 'genbank', 'gnomon', 'pdb']
PROGS = ['blastp', 'tblastn']

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
# --------------------------------------------------
# Utility functions
# --------------------------------------------------
def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"

# 1. Function to validate a hit
def validate_hit(rank=None, k_max=None, eval=None, max_eval=None, identity=None, min_ident=None):
    """
    Validate a BLAST hit based on rank, e-value, and identity thresholds.
    Args:
        rank (int): Rank of the hit.
        k_max (int): Maximum allowable rank.
        eval (float): E-value of the hit.
        max_eval (float): Maximum allowable e-value.
        identity (float): Sequence identity.
        min_identity (float): Minimum allowable identity.

    Returns:
        bool: True if hit is valid, False otherwise.
    """
    return not ((rank is not None and k_max is not None and rank >= k_max) or
                (eval is not None and max_eval is not None and eval > max_eval) or
                (identity is not None and min_ident is not None and identity < min_ident))

def seq_identity(seq1: str, seq2:str)->float:
    """
    Calculate sequence identity between two sequences.
    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.

    Returns:
        float: Fraction of identical characters.
    """
    if not seq1 or not seq2:
        logging.warning("Empty sequence detected in identity calculation.")
        return 0

    short, long = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

    try:
        aligned = edlib.align(short, long, task="path")
        matches = sum(1 for char in edlib.getNiceAlignment(aligned, short, long)["matched_aligned"] if char == "|")
        return matches / len(short)
    except Exception as e:
        logging.error(f"Error in alignment: {e}")
        return 0

# 3. Function to extract BLAST hits into a DataFrame
def extract_hits(xml_res)->pd.DataFrame:
    """
    Parse BLAST XML results into a DataFrame.
    Args:
        xml_res: XML result handle from BLAST.

    Returns:
        pd.DataFrame: DataFrame with hit details.
    """
    try:
        blast_records = NCBIXML.parse(xml_res)
    except Exception as e:
        logging.error(f"Failed to parse BLAST XML: {e}")
        return pd.DataFrame()  # Return empty DataFrame to handle gracefully

    hits = []
    for record in blast_records:
        for aln in record.alignments:
            for hsp in aln.hsps:
                hits.append({
                    "id": aln.accession,
                    "seq": hsp.sbjct,
                    "len": len(hsp.sbjct),
                    "evalue": hsp.expect,
                    "identity": hsp.identities / hsp.align_length * 100
                })
    return pd.DataFrame(hits)

# 4. Function to cluster sequences and select representatives
def cluster_sequences(df, threshold=0.9)->pd.DataFrame:
    """
    Cluster sequences and pick representatives.
    Args:
        df (pd.DataFrame): DataFrame with sequences.
        threshold (float): Similarity threshold for clustering.

    Returns:
        pd.DataFrame: DataFrame with representative sequences.
    """
    """
    Clusters sequences based on similarity and selects representatives.
    """
    seqs = df["seq"].tolist()
    n = len(seqs)

    # === Step 1: Edge Case Handling ===
    if n == 0:
        logging.warning("No sequences provided for clustering.")
        return pd.DataFrame()  # Return empty DataFrame
    if n == 1:
        logging.info("Only one sequence present. No clustering necessary.")
        return df.copy()

    # === Step 2: Compute Pairwise Identity Matrix ===
    identity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):  # Only upper triangle
            identity_matrix[i, j] = seq_identity(seqs[i], seqs[j])
            identity_matrix[j, i] = identity_matrix[i, j]  # Mirror to lower triangle

    # === Step 3: Convert to Condensed Distance Matrix ===
    dst_matrix = 1 - identity_matrix  # Convert similarity to distance
    condensed = squareform(dst_matrix, checks=False)
    link_matrix = linkage(condensed, method="average")
    clusters = fcluster(link_matrix, t=1 - threshold, criterion="distance")

    # Select representative hits
    reps = []
    for cluster in np.unique(clusters):
        indices = np.where(clusters == cluster)[0]
        cluster_df = df.iloc[indices]
        top_hit = cluster_df.sort_values(by=["identity", "evalue", "len"], ascending=[False, True, True]).iloc[0]
        reps.append(top_hit)

    return pd.DataFrame(reps)

# --------------------------------------------------
# API BLAST pipeline
# --------------------------------------------------
def run_blast_api(query_seq, k_max=10, max_eval=1e-5, min_ident=30.0, cluster_thresh=0.9, db="swissprot", prog="blastp"):
    """
    Execute BLAST and process results.
    Args:
        query (str): Query sequence.
        k_max (int): Max hits to keep.
        max_eval (float): Max e-value.
        min_identity (float): Min identity.
        cluster_thresh (float): Clustering threshold.
        db (str): BLAST database.
        prog (str): BLAST program.

    Returns:
        list: Representative hit IDs.
    """
    logging.info(f"Running BLAST on database: {db}...")
    handle = NCBIWWW.qblast(prog, db, query_seq)
    hits_df = extract_hits(handle)

    if hits_df.empty:
        logging.warning("No hits found in BLAST search.")
        return []

    hits_df = hits_df.sort_values("evalue", ascending=True).reset_index(drop=True)
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
        return []

    logging.info(f"{len(hits_df)} hits passed filtering. Proceeding to clustering.")
    clustered_df = cluster_sequences(hits_df, threshold=cluster_thresh)

    logging.info(f"Clustering reduced the hits to {len(clustered_df)} unique representatives.")
    return clustered_df["id"].tolist()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST to get a list of similar genes.")

    parser.add_argument("--fasta", required=True, help="Path to FASTA file containing the query sequence.")
    parser.add_argument("--dbs",nargs='+',default=["swissprot"],help="One or more BLAST databases to query (default: swissprot).")
    parser.add_argument("--prog", default="blastp", help="BLAST program (default: blastp).")
    parser.add_argument("--kMax", type=int, default=10, help="Maximum rank of hits to consider (default: 10).")
    parser.add_argument("--max_eval", type=float, default=1e-5, help="Maximum allowed e-value (default: 1e-5).")
    parser.add_argument("--min_ident", type=float, default=40.0, help="Minimum identity percentage (default: 40.0).")
    parser.add_argument("--cluster", type=float, default=0.9, help="Clustering threshold (default: 0.9).")
    parser.add_argument("--max_len", type=int, default=400, help="Maximum sequence length (default: 400).")
    parser.add_argument("--email", required=False, help="User email address (required for BLAST). If not provided, you will be prompted.")
    parser.add_argument("--output_file", required=True, help="Path to file where results will be saved" )
    args = parser.parse_args()

    try:
        sequence = get_seq(args.fasta, max_len=args.max_len)

    except Exception as e:
        logging.error(f"Error reading FASTA file: {e}")
        sys.exit(1)

    if not args.email:
        print(">>> Please enter your email address (required for NCBI BLAST): ", end="")
        user_input = input().strip()
        if not user_input:
            logging.error("Email address is required for BLAST. Exiting.")
            sys.exit(1)
        args.email = user_input

    if not re.match(r"^[a-zA-Z0-9._%+\-]+@[a-zA-Z0-9.\-]+\.[a-zA-Z]{2,}$", args.email):
        logging.error(f"Invalid email format: {args.email}")
        sys.exit(1)

    # Set email for Entrez
    Entrez.email = args.email

    # Validate the selected program
    if args.prog not in PROGS:
        logging.error(f"Invalid program '{args.prog}'. Available options are: {', '.join(PROGS)}.")
        logging.info("Using default program 'blastp'.")
        args.prog = "blastp"

    valid_dbs = [db for db in args.dbs if db in DBS]

    if not valid_dbs:
        logging.error(f"Invalid databases. Available options are: {', '.join(DBS)}.")
        logging.info("Using default database 'swissprot'.")
        valid_dbs = ["swissprot"]

    timestamp = int(time.time())
    all_results = []
    for db in valid_dbs:
        try:
            results = run_blast_api(
                query_seq=sequence,
                k_max=args.kMax,
                max_eval=args.max_eval,
                min_ident=args.min_ident,
                cluster_thresh=args.cluster,
                db=db,
                prog=args.prog
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
        #logging.info(f"Blast API results saved => {args.output_file}")

    except IOError as e:
        logging.error(f"Error writing results to {args.output_file}: {e}")
        sys.exit(1)
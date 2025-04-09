#!/usr/bin/env python

import os
import obonet
import time
import sys
import argparse
import logging
import networkx as nx
from collections import Counter, defaultdict
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))

sys.path.extend([
    os.path.join(BASE_DIR, 'generate_graph'),
    os.path.join(BASE_DIR, 'generate_graph/graph_edit'),
    os.path.join(BASE_DIR, 'retrieval')
])
from generate_graphviz import generate_graphviz
from retrieval.get_go_terms import load_go_terms_dict
from go_graph_cleaning import process_cleaning
from process_svg import process_svg_file

GO_PIVOTS = {"GO:0003674": "molecular_function",
             "GO:0008150": "biological_process",
             "GO:0005575": "cellular_component"
             }

def format_time(seconds):
    mins, sec = divmod(seconds, 60)
    return f"{mins}m:{sec:02d}s"


def summarize_to_text(uniprot_to_go: dict, go_graph, pivots: list):
    """
    Generate a textual summary of UniProt/ELM - GO annotations.

    Returns:
        summary_text (str): Text summary to print or save.
        summary_dict (dict): Parsed values if JSON output is desired.
    """
    all_go_terms = [term for terms in uniprot_to_go.values() for term in terms]
    go_terms_cnt = Counter(all_go_terms)
    protein_per_go = defaultdict(set)
    reversed_go_graph = go_graph.reverse(copy=True)

    for prot, terms in uniprot_to_go.items():
        for t in terms:
            protein_per_go[t].add(prot)

    total_proteins = len(uniprot_to_go)
    total_terms = len(go_terms_cnt)

    ontology_cnt = {pivot: 0 for pivot in pivots}
    for term in go_terms_cnt:
        if term in go_graph:
            ancestors = nx.ancestors(reversed_go_graph, term)
            for pivot in pivots:
                if pivot in ancestors or pivot == term:
                    ontology_cnt[pivot] += 1

    summary_text = (
        f"GO Annotation Summary:\n"
        f"  - Total proteins: {total_proteins}\n"
        f"  - Total unique GO terms: {total_terms}\n"
        f"\nGO Domain Breakdown (based on pivots):\n"
    )
    for pivot, cnt in ontology_cnt.items():
        label = go_graph.nodes[pivot].get('name', 'Unknown')
        summary_text += f"  - {pivot} ({label}): {cnt} terms\n"

    top_terms = go_terms_cnt.most_common(5)
    summary_text += "\nTop 5 GO terms by frequency:\n"
    for go_id, cnt in top_terms:
        label = go_graph.nodes.get(go_id, {}).get("name", "Unknown")
        summary_text += f"  - {go_id} ({label}): {cnt} annotations\n"

    return summary_text

def process_golizard(uniprot_to_go_path, go_graph_path, pivots, output_dir):
    """
    Main function to process UniProt-to-GO mappings and generate visual GO graphs.

    Args:
        uniprot_to_go_file (str): Path to the file with UniProt-to-GO mappings (JSON format).
        go_graph_path (str): Path to the GO graph file (e.g., .obo file).
        pivots (list of str): List of pivot GO categories (e.g., molecular function, biological process).
        output_dir (str): Directory to save output DOT and SVG files.
    """
    timestamp = int(time.time())

    go_graph = obonet.read_obo(go_graph_path)
    uniprot_to_go = load_go_terms_dict(uniprot_to_go_path)


    # Process GO terms and generate graphs
    logging.info("Graph traversal started")
    cleaned_graphs, pivots, max_dst, selected = process_cleaning(
        uniprot_to_go, go_graph, pivots, verbose=False
    )
    output_graphs = {}

    for ontology in pivots:
        # SUMMARY TEXT FOR EACH SUBGRAPH (GO_TERMS_CNT == NODE_CNT)
        print(f"Nodes {ontology}\n{cleaned_graphs[ontology].nodes}\n")
        output_graphs[ontology] = generate_graphviz(cleaned_graphs[ontology])

        # Define output file paths
        dot_file_path = os.path.join(output_dir, f"{GO_PIVOTS[ontology]}_graphviz.dot")
        svg_file_path = os.path.join(output_dir, f"{GO_PIVOTS[ontology]}_graphviz.svg")

        # Save files
        output_graphs[ontology].write(dot_file_path)

        output_graphs[ontology].draw(path=svg_file_path, format="svg", prog="dot")

        # Process the SVG file for enhancements
        process_svg_file(svg_file_path)

    elapsed = int(time.time()) - timestamp
    logging.info(f"Elapsed: {format_time(elapsed)}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process UniProt-to-GO mappings and generate GO graphs.")
    parser.add_argument(
        "--uniprot_to_go", required=True, help="Path to UniProt-to-GO mapping file (JSON format)."
    )
    parser.add_argument(
        "--go_graph", required=True, help="Path to GO graph file (e.g., .obo file)."
    )
    parser.add_argument(
        "--pivot", nargs="+", required=True,
        help="List of GO pivot categories (e.g., GO:0003674, GO:0008150, GO:0005575).",
        default=["GO:0003674", "GO:0008150", "GO:0005575"]
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to save output files (DOT and SVG)."
    )

    args = parser.parse_args()

    uniprot_path = Path(args.uniprot_to_go)
    if not uniprot_path.is_file():
        logging.error(f"Invalid path to uniprot-go mapping: {uniprot_path}")
        sys.exit(1)

    go_path = Path(args.go_graph)
    if not go_path.is_file():
        logging.error(f"Invalid path to go graph: {go_path}")
        sys.exit(1)
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception as e:
        print(f"[ERROR] Failed to create output directory: {args.output_dir}\n{e}", file=sys.stderr)
        sys.exit(1)

    # Call the main processing function
    process_golizard(args.uniprot_to_go, args.go_graph, args.pivot, args.output_dir)
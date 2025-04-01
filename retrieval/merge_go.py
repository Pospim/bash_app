#!/usr/bin/env python
"""_summary_

    Returns:
        _type_: _description_
"""
import argparse
import json
import sys
import pandas as pd
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def load_json(file_path):
    path = Path(file_path)
    if not path.is_file():
        logging.error(f"File does not exist: {file_path}")
        sys.exit(1)

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        logging.error(f"Failed to parse JSON: {file_path}\n{e}")
        sys.exit(1)

    if not isinstance(data,dict):
        logging.error(f"JSON file must contain a dictionary: {file_path}")
        sys.exit(1)

    return data

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

def merge_dicts(fs_data, elm_data):
    merged = elm_data.copy()
    merged.update(fs_data) # Add/replace with FoldSeek
    return merged

def merge_csv(fs_data, elm_data):
    if fs_data.empty and elm_data.empty:
        logging.error(f"Both .csv files are empty: {fs_data}--{elm_data}")
        return []
    if fs_data.empty:
        return elm_data
    if elm_data.empty:
        return fs_data

    try:
        fs_data = fs_data.reset_index(drop=True)
        elm_data = elm_data.reset_index(drop=True)

        fs_data.columns = ['GO_term', 'score']
        elm_data.columns = ['GO_term', 'score']

        merged_df = pd.concat([fs_data, elm_data], ignore_index=True)

        merged_df = merged_df.drop_duplicates(subset='GO_term', keep='first')
        return merged_df

    except Exception as e:
        logging.error(f"Failed to merge csv files: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Merge GO term JSON files (ELM + FoldSeek)")
    parser.add_argument("--elm_json", required=True, help="ELM GO terms JSON file (elm_id: [go_terms])")
    parser.add_argument("--elm_csv", required=True, help="ELM GO terms with score csv file (go_term, score)")
    parser.add_argument("--fs_json", required=True, help="FoldSeek/BLAST GO terms JSON file (protein_id: [go_terms])")
    parser.add_argument("--fs_csv", required=True, help="FoldSeek/BLAST GO terms with score csv file (go_term, score)")
    parser.add_argument("--output_json", default="merged_go_score.json", help="Path to output merged JSON (default: merged_go_score.csv)")
    parser.add_argument("--output_csv", default="merged_go_score.csv", help="Path to output merged csv (default: merged_go_score.csv)")
    args = parser.parse_args()

    elm_dict = load_json(args.elm_json)
    fs_dict = load_json(args.fs_json)
    merged_dict = merge_dicts(fs_dict, elm_dict)

    elm_df = load_csv(args.elm_csv)
    fs_df = load_csv(args.fs_csv)
    merged_df = merge_csv(elm_df, fs_df)

    try:
        with open(args.output_json, 'w') as out:
            json.dump(merged_dict, out, indent=2)
        merged_df.to_csv(args.output_csv, sep='\t', index=False)
        print("[INFO] Merge completed successfully.")
    except Exception as e:
        print(f"[ERROR] GO mapping merge failed: {e}")

if __name__ == "__main__":
    main()

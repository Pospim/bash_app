#!/usr/bin/env python
"""
merge_go_terms.py

Merge two GO term mapping files in JSON format:

Priority:
- BLAST GO terms take precedence if there are overlapping IDs.

Usage:
    python merge_go_terms.py --elm <elm_go_terms.json> --fs <foldseek_go_terms.json> --output <merged_go_terms.json>
"""

import argparse
import json
import sys
import os
from pathlib import Path

def load_json(file_path):
    path = Path(file_path)
    if not path.is_file():
        print(f"[ERROR] File does not exist: {file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"[ERROR] Failed to parse JSON: {file_path}\n{e}", file=sys.stderr)
        sys.exit(1)

    if not isinstance(data,dict):
        print(f"[ERROR] JSON file must contain a dictionary: {file_path}", file=sys.stderr)
        sys.exit(1)

    return data

def merge_dicts(fs_data, elm_data):
    merged = elm_data.copy()
    merged.update(fs_data) # Add/replace with FoldSeek
    return merged

def main():
    parser = argparse.ArgumentParser(description="Merge GO term JSON files (ELM + FoldSeek)")
    parser.add_argument("--elm", required=True, help="ELM GO terms JSON file (elm_id: [go_terms])")
    parser.add_argument("--fs", required=True, help="FoldSeek/BLAST GO terms JSON file (protein_id: [go_terms])")
    parser.add_argument("--output", required=True, help="Path to output merged JSON")

    args = parser.parse_args()

    elm_data = load_json(args.elm)
    fs_data = load_json(args.fs)
    merged_data = merge_dicts(fs_data, elm_data)

    print(f"[INFO] Saving merged GO terms to: {args.output}")
    output_dir = os.path.dirname(os.path.abspath(args.output))

    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"[ERROR] Failed to create output directory: {output_dir}\n{e}", file=sys.stderr)
        sys.exit(1)
    try:
        with open(args.output, 'w') as out:
            json.dump(merged_data, out, indent=2)
        print("[INFO] Merge completed successfully.")
    except Exception as e:
        print(f"[ERROR] GO mapping merge failed: {e}")

if __name__ == "__main__":
    main()

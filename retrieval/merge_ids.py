import argparse
import sys
import os
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def file_nonempty(file_path: str) -> bool:
    path = Path(file_path)
    return path.is_file() and os.path.getsize(path) > 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Unify 2 results into a single file."
    )
    parser.add_argument("--local_file", required=True,
                        help="Path to local results file (one ID per line).")
    parser.add_argument("--api_file", required=True,
                        help="Path to remote API results file (one ID per line).")
    parser.add_argument("--output_file", required=True,
                        help="Path to the final merged results file.")
    args = parser.parse_args()

    local_exists = file_nonempty(args.local_file)
    api_exists = file_nonempty(args.api_file)

    if not local_exists and not api_exists:
        logging.error("No results found")
        sys.exit(0)

    unique_ids = set()

    if local_exists:
        with open(args.local_file, 'r') as local:
            for line in local:
                line = line.strip()
                if line:
                    unique_ids.add(line)

    if api_exists:
        with open(args.api_file, 'r') as api:
            for line in api:
                line = line.strip()
                if line:
                    unique_ids.add(line)

    if not unique_ids:
        logging.error("No IDs found in results files")
        sys.exit(1)

    with open(args.output_file, 'w') as outf:
        for id in sorted(unique_ids):
            outf.write(id + '\n')

    if not file_nonempty(args.output_file):
        logging.error(f"[ERROR] Output file {args.output_file} was not created or is empty.")
        sys.exit(1)

    print(f"[INFO] Merged {len(unique_ids)} unique IDs => {args.output_file}")
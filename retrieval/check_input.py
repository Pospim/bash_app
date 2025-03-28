import argparse
import sys
import pyfastaq
from pathlib import Path

# Valid amino acid characters
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # Standard single-letter amino acid codes
MAX_SEQ_LEN = 400 # max allowed length for FoldSeek

def is_valid_seq(seq: str, max_len: int=400) -> bool:
    """
    Check if a sequence:
    - Contains only valid amino acid characters.
    - Has length between 1 and max_len.
    Args:
        seq (str): Amino acid sequence.
        max_len (int): Maximum sequence length.

    Returns:
        bool: True if valid, False otherwise.
    """
    return 0 < len(seq) <= max_len and all(char in AMINO_ACIDS for char in seq.upper())

def parse_fasta(file_path: Path, max_len: int) -> str:
    """
    Extract the first valid amino acid sequence from a FASTA file.

    Args:
        file_path (path): Path to the FASTA file.
        max_len (int): Maximum allowed sequence length.

    Returns:
        str: The first valid sequence found.

    Raises:
        ValueError: If no valid sequences are found.
    """
    valid_seq = None

    for record in pyfastaq.sequences.file_reader(str(file_path)):
        if is_valid_seq(record.seq, max_len):
            valid_seq = record.seq
            break

    if not valid_seq:
        raise ValueError(f"No valid sequences found in {file_path} (max_len={max_len}).")

    return valid_seq


def get_seq(input_path: Path = None, max_len: int = MAX_SEQ_LEN) -> str:
    """
    Validate and retrieve an amino acid sequence from a file.

    Args:
        input_path (str): Path to the input file.
        max_len (int): Maximum allowed sequence length.

    Returns:
        str: Validated amino acid sequence.

    Raises:
        ValueError: If input is invalid or no sequences are found.
    """
    if max_len <= 0:
        raise ValueError("max_len must be greater than 0.")

    if input_path:
        return parse_fasta(input_path, max_len)
    else:
        raise ValueError("No input provided. Provide a file path or a pasted sequence.")

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Extract an amino acid sequence from a FASTA file.")
        parser.add_argument("-f", "--file", type=str, required=True, help="Path to the FASTA file.")
        parser.add_argument("-m", "--max_len", type=int, default=MAX_SEQ_LEN, help=f"Maximum sequence length (default: 400).")
        args = parser.parse_args()

        file_path = Path(args.file)
        if not file_path.is_file():
            print(f"Error: File '{args.file}' does not exist.")
            sys.exit(1)

        seq = get_seq(input_path=file_path, max_len=args.max_len)
        print(seq)
    except ValueError as err:
        print(f"Error reading fasta file: {err}", file=sys.stderr)
        sys.exit(1)
    except pyfastaq.sequences.Error as err:
        print(f"[ERROR] Invalid FASTA format: {err}", file=sys.stderr)
        sys.exit(1)
    except Exception as err:
        print(f"[ERROR] Unexpected error: {err}", file=sys.stderr)
        sys.exit(1)

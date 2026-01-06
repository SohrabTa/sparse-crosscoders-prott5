
import gzip
import argparse
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm

def count_proteins(
    input_path: str,
    max_length: int = 512
) -> None:
    """
    Counts the number of sequences in a gzipped FASTA file based on length.

    Args:
        input_path: Path to the input gzipped FASTA file.
        max_length: Maximum length of sequences to include.
    """
    input_file = Path(input_path)
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Processing {input_path}...")
    print(f"Target: sequences with length <= {max_length}")

    count = 0
    
    with gzip.open(input_file, "rt") as handle_in:

        records = SeqIO.parse(handle_in, "fasta")
        # total extracted from uniref50.release_note
        pbar = tqdm(total=32474829, desc="Sequences processed", unit="seq")
        
        for record in records:
            pbar.update(1)
            if len(record.seq) <= max_length:
                count += 1
                
        pbar.close()

    print(f"Found {count} sequences with length <= {max_length}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Counts the number of sequences in a gzipped FASTA file based on length.")
    parser.add_argument("--input", default="data/uniref50/uniref50.fasta.gz", help="Path to input .fasta.gz file")
    parser.add_argument("--max_length", type=int, default=512, help="Maximum sequence length")
    
    args = parser.parse_args()
    
    count_proteins(args.input, args.max_length)

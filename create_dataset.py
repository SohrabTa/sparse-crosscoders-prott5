import argparse
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm

def create_dataset(
    input_path: str,
    output_path: str,
    max_sequences: int = 3_000_000,
    max_length: int = 512
) -> None:
    """
    Extracts a subset of sequences from a FASTA file based on length.

    Args:
        input_path: Path to the input FASTA file.
        output_path: Path to save the extracted sequences in FASTA format.
        max_sequences: Maximum number of sequences to extract.
        max_length: Maximum length of sequences to include.
    """
    input_file = Path(input_path)
    output_file = Path(output_path)
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Processing {input_path}...")
    print(f"Target: {max_sequences} sequences with length <= {max_length}")
    print(f"Output: {output_path}")

    count = 0
    
    # ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "rt") as handle_in, \
         open(output_file, "w") as handle_out:
        
        seqs_processed_pbar = tqdm(total=32474829, desc="Sequences processed", unit="seq")
        seqs_collected_pbar = tqdm(total=max_sequences, desc="Sequences collected", unit="seq")

        for record in SeqIO.parse(handle_in, "fasta"):
            seqs_processed_pbar.update(1)
            if len(record.seq) <= max_length:
                SeqIO.write(record, handle_out, "fasta")
                count += 1
                seqs_collected_pbar.update(1)
                
                if count >= max_sequences:
                    break
        
        seqs_processed_pbar.close()
        seqs_collected_pbar.close()

    if count < max_sequences:
        print(f"Warning: Only found {count} sequences matching criteria (requested {max_sequences}).")
    else:
        print(f"Successfully extracted {count} sequences.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a dataset of protein sequences from UniRef50.")
    parser.add_argument("--input", default="data/uniref50/uniref50.fasta", help="Path to input .fasta file")
    parser.add_argument("--output", default="data/uniref50/uniref50_3M_length_512.fasta", help="Path to output .fasta file")
    parser.add_argument("--limit", type=int, default=3_000_000, help="Number of sequences to extract")
    parser.add_argument("--max_length", type=int, default=512, help="Maximum sequence length")
    
    args = parser.parse_args()
    
    create_dataset(args.input, args.output, args.limit, args.max_length)

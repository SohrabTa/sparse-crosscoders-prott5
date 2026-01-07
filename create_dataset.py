import argparse
import random
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm

def filter_sequences(
    input_path: Path,
    output_path: Path,
    max_length: int
) -> int:
    """
    Filters sequences by length and saves them to an intermediate file.
    If the output file already exists, the entries will be counted and returned.
    Returns the count of valid sequences found.
    """
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    count = 0
    if output_path.exists():
        print(f"Found existing output file at {output_path}. Counting sequences...")
        with open(output_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                count += 1
        return count
    
    print(f"Filtering sequences (length <= {max_length}) from {input_path}...")
    with open(input_path, "rt") as handle_in, \
         open(output_path, "w") as handle_out:
        
        # We don't know total sequences in input easily without reading, so maybe just show processed count
        pbar = tqdm(desc="Filtering sequences", unit="seq")
        
        for record in SeqIO.parse(handle_in, "fasta"):
            pbar.update(1)
            if len(record.seq) <= max_length:
                SeqIO.write(record, handle_out, "fasta")
                count += 1
        
        pbar.close()
        
    print(f"Found {count} valid sequences.")
    return count

def sample_sequences(
    input_path: Path,
    output_path: Path,
    total_count: int,
    max_sequences: int,
    seed: int | None
) -> None:
    """
    Randomly samples sequences from the filtered dataset.
    """
    if seed is not None:
        random.seed(seed)
        
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if total_count <= max_sequences:
        print(f"Total valid sequences ({total_count}) is not greater than limit ({max_sequences}).")
        print("Copying all valid sequences to output...")
        indices_to_keep = set(range(total_count))
    else:
        print(f"Sampling {max_sequences} sequences from {total_count} valid candidates...")
        indices_to_keep = set(random.sample(range(total_count), max_sequences))
    
    print(f"Writing result to {output_path}")
    
    with open(input_path, "rt") as handle_in, \
         open(output_path, "w") as handle_out:
        
        pbar = tqdm(total=max_sequences, desc="Sampling sequences", unit="seq")
        
        written_count = 0
        for i, record in enumerate(SeqIO.parse(handle_in, "fasta")):
            if written_count >= max_sequences:
                break   
            if i in indices_to_keep:
                SeqIO.write(record, handle_out, "fasta")
                pbar.update(1)
                written_count += 1
            
        pbar.close()

    print(f"Successfully created dataset with {max_sequences} sequences.")

def create_dataset(
    input_path: str,
    output_path: str,
    max_sequences: int,
    max_length: int,
    seed: int | None
) -> None:
    """
    Extracts a subset of randomly sampled sequences from a FASTA file based on length.
    
    Process:
    1. Filter all sequences <= max_length to an intermediate file.
    2. Randomly sample from the intermediate file.
    """
    input_file = Path(input_path)
    output_file = Path(output_path)
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Determine intermediate filename
    # e.g., data/uniref50/uniref50.fasta -> data/uniref50/uniref50_length_512.fasta
    stem = input_file.stem
    start_dir = input_file.parent
    suffix = input_file.suffix
    intermediate_name = f"{stem}_length_{max_length}{suffix}"
    intermediate_file = start_dir / intermediate_name
    
    # Step 1: Filter
    total_valid = filter_sequences(input_file, intermediate_file, max_length)
    
    # Step 2: Sample
    sample_sequences(intermediate_file, output_file, total_valid, max_sequences, seed)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a dataset of protein sequences from UniRef50 with random sampling.")
    parser.add_argument("--input", default="data/uniref50/uniref50.fasta", help="Path to input .fasta file")
    parser.add_argument("--output", default="data/uniref50/uniref50_3M_length_512.fasta", help="Path to output .fasta file")
    parser.add_argument("--limit", type=int, default=3_000_000, help="Number of sequences to extract")
    parser.add_argument("--max_length", type=int, default=512, help="Maximum sequence length")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    create_dataset(args.input, args.output, args.limit, args.max_length, args.seed)

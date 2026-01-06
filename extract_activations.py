import torch
from transformers import T5EncoderModel, T5Tokenizer
import re
import argparse
from typing import List, Tuple
from Bio import SeqIO
from pathlib import Path

def get_device() -> torch.device:
    """Determines the best available device.

    Returns:
        torch.device: The best available device (CUDA or MPS).
    """
    if torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
        return torch.device("cuda")
    elif torch.backends.mps.is_available():
        print("Using Apple Silicon (MPS)")
        return torch.device("mps")
    else:
        raise Exception("No GPU (CUDA or MPS) found!")

def load_model(model_name: str, device: torch.device, dtype: torch.dtype = torch.float16) -> Tuple[T5EncoderModel, T5Tokenizer]:
    """Loads the T5 model and tokenizer.

    Args:
        model_name: The name of the model to load from HuggingFace.
        device: The device to load the model onto.
        dtype: The data type to use for the model.

    Returns:
        A tuple containing the loaded model and tokenizer.
    """
    print(f"Loading model: {model_name}...")
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name, dtype=dtype)
    model = model.to(device)
    model.eval()
    print("Model loaded successfully.")
    return model, tokenizer

def load_sequences_from_fasta(file_path: str) -> List[str]:
    """Loads sequences from a FASTA file.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        list of protein sequences as strings.
        
    Raises:
        FileNotFoundError: If the file does not exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    
    print(f"Loading sequences from {file_path}...")
    records = list(SeqIO.parse(path, "fasta"))
    sequences = [str(record.seq) for record in records]
    print(f"Loaded {len(sequences)} sequences.")
    return sequences

def extract_activations(
    model: T5EncoderModel, 
    tokenizer: T5Tokenizer, 
    sequences: List[str], 
    device: torch.device
) -> torch.Tensor:
    """Extracts activations from the model for the given sequences.

    Args:
        model: The T5 encoder model.
        tokenizer: The T5 tokenizer.
        sequences: A list of protein sequences.
        device: The device to run inference on.

    Returns:
        A tensor containing the extracted activations.
        Shape: (num_layers, batch_size, seq_len, hidden_dim).
        Note: Layers includes the initial embedding layer.
    """
    # Pre-processing (Regex replace UZOB -> X, add spaces)
    processed_seqs = [" ".join(list(re.sub(r"[UZOB]", "X", seq))) for seq in sequences]

    # Tokenize
    ids = tokenizer.batch_encode_plus(processed_seqs, add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)

    print(f"Input IDs shape: {input_ids.shape}")

    # Forward Pass
    print("Running inference...")
    with torch.no_grad():
        output = model(input_ids=input_ids, attention_mask=attention_mask, output_hidden_states=True)

    # Extract Hidden States (excluding initial embedding layer)
    # output.hidden_states is a tuple of (embedding_output, layer_1, ..., layer_N)
    # We want all except the very first one (which is the embedding layer)
    all_hidden_states = output.hidden_states[1:25]

    # Stack them: (num_layers, batch_size, seq_len, hidden_dim)
    stacked_activations = torch.stack(all_hidden_states)

    # Verify shape
    print(f"Extracted activations shape: {stacked_activations.shape}")
    print("(Layers, Batch_Size, Sequence_Length, Hidden_Dim)")
    
    return stacked_activations

def main() -> None:
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Extract activations from ProtT5 model.")
    parser.add_argument("--batch_size", type=int, default=5, help="Number of sequences to process from the file.")
    parser.add_argument("--input", type=str, default="data/uniref50/uniref50_3M_length_512.fasta", help="Input FASTA file path.")
    parser.add_argument("--output", type=str, default="random_protein_activations_batch.pt", help="Output file path.")
    args = parser.parse_args()

    device = get_device()
    model_name = 'Rostlab/prot_t5_xl_half_uniref50-enc'
    model, tokenizer = load_model(model_name, device)
    
    all_sequences = load_sequences_from_fasta(args.input)
    if not all_sequences:
        print("No sequences found in the file.")
        return
        
    # Process only the first batch_size sequences
    sequences = all_sequences[:args.batch_size]
    print(f"Processing first {len(sequences)} sequences...")

    activations = extract_activations(model, tokenizer, sequences, device)
    
    torch.save(activations, args.output)
    print(f"Saved activations to {args.output}")

if __name__ == "__main__":
    main()



import torch
from transformers import T5EncoderModel, T5Tokenizer
import requests
import random
import re
import argparse
from typing import List, Tuple

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
    model = T5EncoderModel.from_pretrained(model_name, torch_dtype=dtype)
    model = model.to(device)
    model.eval()
    print("Model loaded successfully.")
    return model, tokenizer

def get_random_uniref50_sequences(batch_size: int = 5) -> List[str]:
    """Fetches random sequences from UniRef50.

    Args:
        batch_size: The number of sequences to fetch. Defaults to 5.

    Returns:
        list of protein sequences as strings.
        
    Raises:
        ValueError: If no results are found from the UniProt API.
    """
    # TODO replace this with loading from preprocessed fasta file 
    # Filter: Modified before 2019, length <= 512
    # UniRef search API
    url = "https://rest.uniprot.org/uniref/search"
    query = "date_modified:[* TO 2019-01-01] AND length:[1 TO 512] identity:0.5"
    
    params = {
        'query': query,
        'format': 'json',
        'size': 50  # Fetch a larger pool to sample from
    }
    
    print(f"Fetching random sequences (batch_size={batch_size}) from UniRef50...")
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()
    
    if 'results' not in data or not data['results']:
        raise ValueError("No results found for the query.")
    
    # Pick random entries
    entries = random.sample(data['results'], min(batch_size, len(data['results'])))
    
    sequences = []
    for entry in entries:
        # Extract sequence
        try:
            seq = entry['representativeMember']['sequence']['value']
            accession = entry['id']
            print(f"Selected UniRef50 Entry: {accession}, Length: {len(seq)}")
            sequences.append(seq)
        except KeyError:
            # Fallback if structure is different
            print("Could not parse sequence from entry, dumping keys:", entry.keys())
            continue
            
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
    print("Note: Layers does not include the initial embedding layer.")
    
    return stacked_activations

def main() -> None:
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Extract activations from ProtT5 model.")
    parser.add_argument("--batch_size", type=int, default=5, help="Number of sequences to process.")
    parser.add_argument("--output", type=str, default="random_protein_activations_batch.pt", help="Output file path.")
    args = parser.parse_args()

    device = get_device()
    model_name = 'Rostlab/prot_t5_xl_half_uniref50-enc'
    model, tokenizer = load_model(model_name, device)
    
    sequences = get_random_uniref50_sequences(batch_size=args.batch_size)
    if not sequences:
        print("No sequences found/parsed.")
        return

    activations = extract_activations(model, tokenizer, sequences, device)
    
    torch.save(activations, args.output)
    print(f"Saved activations to {args.output}")

if __name__ == "__main__":
    main()

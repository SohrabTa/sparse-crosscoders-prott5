import urllib.request
import urllib.parse
import os
import ssl
from pathlib import Path

def get_uniprot_url(annotation_scores=[5], max_length=512):
    """
    Constructs the UniProt REST API URL for downloading a tailored Swiss-Prot dataset.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    
    # Define the fields to retrieve (matching InterPLM's extract_annotations.py expectations)
    fields = [
        "accession", "reviewed", "protein_name", "length", "sequence", "ec", 
        "ft_act_site", "ft_binding", "cc_cofactor", "ft_disulfid", "ft_carbohyd", 
        "ft_lipid", "ft_mod_res", "ft_signal", "ft_transit", "ft_helix", "ft_turn", 
        "ft_strand", "ft_coiled", "cc_domain", "ft_compbias", "ft_domain", "ft_motif", 
        "ft_region", "ft_zn_fing", "xref_alphafolddb"
    ]
    
    # Define query parameters
    query_parts = [
        "(reviewed:true)",          # Swiss-Prot only
        "(database:alphafolddb)",   # Has AlphaFold DB cross-reference
        f"(length:[1 TO {max_length}])"
    ]
    
    if annotation_scores:
        scores_query = " OR ".join([f"(annotation_score:{s})" for s in annotation_scores])
        query_parts.append(f"({scores_query})")

    query = " AND ".join(query_parts)
    
    # Construct final URL
    params = {
        "compressed": "true",
        "fields": ",".join(fields),
        "format": "tsv",
        "query": query
    }
    
    url = f"{base_url}?{urllib.parse.urlencode(params)}"
    return url

def download_dataset(output_path, annotation_scores=[5], max_length=512):
    url = get_uniprot_url(annotation_scores, max_length)
    print(f"Constructed Query URL:\n{url}\n")
    
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Downloading dataset to {out_file}...")
    
    # Disable SSL verification issues if running locally on some macs, but normally not needed
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    
    try:
        with urllib.request.urlopen(url, context=ctx) as response, open(out_file, 'wb') as out:
            data = response.read()
            out.write(data)
        print("\nDownload complete!")
        print(f"File size: {os.path.getsize(out_file) / (1024*1024):.2f} MB")
        
    except Exception as e:
        print(f"Error downloading dataset: {e}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Download UniProt dataset for InterPLM annotation processing.")
    parser.add_argument("--output", type=str, default="data/uniprotkb_modern/proteins.tsv.gz",
                        help="Path to save the compressed TSV dataset.")
    parser.add_argument("--score", type=int, nargs='+', default=[5],
                        help="Annotation score(s) to fetch (1-5). e.g. --score 4 5. Use 0 to fetch all scores.")
    parser.add_argument("--max-length", type=int, default=512,
                        help="Maximum sequence length.")
    
    args = parser.parse_args()
    scores = None if 0 in args.score else args.score
    
    download_dataset(args.output, scores, args.max_length)

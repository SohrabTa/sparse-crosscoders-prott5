#!/usr/bin/env python3
"""
Build the full-UniRef50 (<=512 AA) training corpus, split into walltime-sized
chunks for chunked/resumed crosscoder training.

Why chunks
----------
The crosscoder trainer is single-GPU and (with the resume patch) can only resume
at fasta-file boundaries, not mid-stream. The LRZ walltime cap is 34h and a
~1.29B-residue corpus took ~31h, so we split the ~5.63B-residue <=512 corpus into
K chunks each sized to run comfortably under walltime. Each chunk is trained in
one Slurm job that resumes from the previous chunk's checkpoint.

Why random assignment
---------------------
Each <=512 sequence is assigned to a chunk by a seeded RNG in a single streaming
pass. This (a) balances chunks by expectation and (b) globally shuffles cluster
membership across chunks, which fixes the "first 3M" file-order bias of the
original corpus (data/.../uniref50_length_512_first_3M.fasta) at the same time.
Within a chunk the trainer's 400k shuffle buffer handles local ordering.

Reproducibility
---------------
- Input: --input uniref50.fasta (UniRef50 2019_01 release).
- Outputs: <out_dir>/chunk_{i:02d}.fasta for i in [0,K), plus manifest.json with
  per-chunk sequence/residue/step counts and the global step budget.
- Determinism: seeded (--seed, default 42); logged in the manifest. Re-running
  with the same seed + input reproduces identical chunks.
- No dedup: UniRef50 representatives are already unique per 50%-identity cluster.

Usage
-----
    python chunk_uniref_corpus.py \
      --input data/external/uniprot/release-2019_01/uniref/uniref50/uniref50.fasta \
      --out_dir data/external/uniprot/release-2019_01/uniref/uniref50/chunks_512 \
      --max_length 512 --n_chunks 7 --batch_size 512 --seed 42

Disposition: COMMIT (defines the training corpus for the preprint retrain).
"""
import argparse
import json
import random
from pathlib import Path


def iter_fasta(path):
    """Stream (header_line, seq_len, [seq_lines]) records without buffering the file."""
    header, seq_lines, seq_len = None, [], 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_len, seq_lines
                header, seq_lines, seq_len = line, [], 0
            else:
                seq_lines.append(line)
                seq_len += len(line) - 1  # drop newline
    if header is not None:
        yield header, seq_len, seq_lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--max_length", type=int, default=512)
    ap.add_argument("--n_chunks", type=int, default=7)
    ap.add_argument("--batch_size", type=int, default=512,
                    help="residues per training step (for the step-budget estimate)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    handles = [open(args.out_dir / f"chunk_{i:02d}.fasta", "w") for i in range(args.n_chunks)]
    n_seq = [0] * args.n_chunks
    n_res = [0] * args.n_chunks
    total_seen = 0

    for header, seq_len, seq_lines in iter_fasta(args.input):
        total_seen += 1
        if seq_len > args.max_length or seq_len == 0:
            continue
        c = rng.randrange(args.n_chunks)
        h = handles[c]
        h.write(header)
        h.writelines(seq_lines)
        n_seq[c] += 1
        n_res[c] += seq_len

    for h in handles:
        h.close()

    total_seq = sum(n_seq)
    total_res = sum(n_res)
    chunks = []
    for i in range(args.n_chunks):
        steps = n_res[i] // args.batch_size
        chunks.append({
            "chunk": i,
            "file": f"chunk_{i:02d}.fasta",
            "n_sequences": n_seq[i],
            "n_residues": n_res[i],
            "steps_one_epoch": steps,
        })

    manifest = {
        "input": str(args.input),
        "max_length": args.max_length,
        "n_chunks": args.n_chunks,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "total_sequences_seen": total_seen,
        "total_sequences_kept": total_seq,
        "total_residues_kept": total_res,
        "global_num_steps_one_epoch": total_res // args.batch_size,
        "chunks": chunks,
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(json.dumps(manifest, indent=2))
    print(f"\nglobal_num_steps_one_epoch = {manifest['global_num_steps_one_epoch']:,} "
          f"(set train.num_steps to this)")
    print(f"per-chunk steps ~ {chunks[0]['steps_one_epoch']:,} "
          f"(~{chunks[0]['steps_one_epoch']/81300:.1f} h at 81.3k steps/h)")


if __name__ == "__main__":
    main()

"""
embeddings.py

Traverses subfolders in `--input_dir` for .fna/.fasta + .gff/.gtf pairs.
For each pair (same base name):
  1) Extract protein sequences via get_protein_sequences().
  2) Compute ESM-2 embeddings for proteins below `--length_threshold`.
  3) Generate metadata (normalized length + 20 AA proportions).
  4) Save:
       *_embeddings.txt  (2D array)
       *_metadata.txt    (2D array)
       *_genename.txt    (1D list of gene names)
in the `--output_dir`, preserving the subfolder structure.
"""

import os
import sys
import gc
import re
import math
import glob
import argparse
import numpy as np
from tqdm import tqdm

import torch
import esm
from Bio import SeqIO
from collections import defaultdict, Counter, OrderedDict


from utils import get_protein_sequences


def main():
    parser = argparse.ArgumentParser(description="Compute ESM-2 embeddings + metadata from FASTA/GFF pairs.")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to input directory with subfolders.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to output directory.")
    parser.add_argument("-l", "--length_threshold", type=int, default=1500, help="Maximum protein length to include.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    length_thr = args.length_threshold

    # -------------------------------------------------------------------------
    # 1) Load ESM-2 model
    # -------------------------------------------------------------------------
    print("Loading ESM-2 model...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # no dropout

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)
    model = model.to(device)

    # -------------------------------------------------------------------------
    # 2) Recursively scan subfolders for FASTA + GFF pairs
    # -------------------------------------------------------------------------
    print(f"Searching subfolders of '{input_dir}' for .fasta/.fna + .gff/.gtf pairs...")

    # For convenience, we'll find all possible .fna/.fasta and .gff/.gtf files:
    fasta_candidates = []
    gff_candidates = []
    # Use os.walk to explore subfolders
    for root, dirs, files in os.walk(input_dir):
        for f in files:
            # Check if it's .fna or .fasta
            if f.lower().endswith(".fna") or f.lower().endswith(".fasta"):
                fasta_candidates.append(os.path.join(root, f))
            # Check if it's .gff or .gtf
            elif f.lower().endswith(".gff") or f.lower().endswith(".gtf"):
                gff_candidates.append(os.path.join(root, f))

    # Create a dict: base_name -> (paths to all fna/fasta)
    # We'll match base_name ignoring extension
    # e.g. "something.fna" => base_name "something"
    # We'll then pair up if there's a matching GFF
    def extract_base(filepath):
        return os.path.splitext(os.path.basename(filepath))[0]

    fasta_map = {}
    for ff in fasta_candidates:
        base = extract_base(ff)
        fasta_map.setdefault(base, []).append(ff)

    gff_map = {}
    for gg in gff_candidates:
        base = extract_base(gg)
        gff_map.setdefault(base, []).append(gg)

    # Now find pairs
    paired = {}
    for base in fasta_map:
        if base in gff_map:
            # We'll just pick the first match if there's more than one
            paired[base] = (fasta_map[base][0], gff_map[base][0])

    print(f"Found {len(paired)} paired FASTA/GFF sets.")

    # -------------------------------------------------------------------------
    # 3) Process each paired set, compute embeddings + metadata (+ gene names)
    # -------------------------------------------------------------------------
    for base, (fasta_file, gff_file) in paired.items():
        print(f"\n[Pair: {base}]\n  FASTA: {fasta_file}\n  GFF: {gff_file}")

        # Derive subfolder path for output
        rel_dir = os.path.relpath(os.path.dirname(fasta_file), input_dir)
        out_subdir = os.path.join(output_dir, rel_dir)
        os.makedirs(out_subdir, exist_ok=True)

        out_embeddings = os.path.join(out_subdir, f"{base}_embeddings.txt")
        out_metadata = os.path.join(out_subdir, f"{base}_metadata.txt")
        out_genenames = os.path.join(out_subdir, f"{base}_genename.txt")

        # Skip if embeddings file exists
        if os.path.exists(out_embeddings):
            print("  Embedding file already exists. Skipping.")
            continue

        try:
            # Extract all protein sequences + gene info
            protein_seqs, gene_names, gene_lengths, aa_props = get_protein_sequences(fasta_file, gff_file)
            if not protein_seqs:
                print(f"  No proteins found for {base} or zero-length. Skipping.")
                continue

            print(f"  Extracted {len(protein_seqs)} total protein sequences.")

            # Collect final embeddings / metadata / gene names for sequences < length_thr
            final_embeddings = []
            final_meta = []
            final_gene_names = []

            max_len = max(gene_lengths) if gene_lengths else 1

            for idx, seq in tqdm(enumerate(protein_seqs), total=len(protein_seqs), desc="Proteins"):
                seq_len = len(seq)
                if seq_len >= length_thr:
                    continue  # skip

                # ESM embedding
                protein_name = f"{base}_prot_{idx+1}"
                data = [(protein_name, str(seq))]
                batch_labels, batch_strs, batch_tokens = batch_converter(data)
                batch_tokens = batch_tokens.to(device)

                with torch.no_grad():
                    results = model(batch_tokens, repr_layers=[33], return_contacts=False)
                token_reprs = results["representations"][33]
                # average across tokens except special start/end
                actual_seq_len = (batch_tokens != alphabet.padding_idx).sum(1)[0]
                seq_repr = token_reprs[0, 1 : actual_seq_len - 1].mean(0).cpu().numpy()
                final_embeddings.append(seq_repr)

                # Build row of metadata: [normalized length + 20 AA proportions]
                length_norm = seq_len / max_len if max_len else 0
                row_meta = [length_norm] + aa_props[idx]
                final_meta.append(row_meta)

                # Store the gene name for reference
                final_gene_names.append(gene_names[idx])

                # Clean up GPU mem
                del batch_tokens, token_reprs, results
                torch.cuda.empty_cache()
                gc.collect()

            if len(final_embeddings) == 0:
                print(f"  No proteins under length {length_thr}. Skipping.")
                continue

            # Convert to arrays
            embeddings_matrix = np.vstack(final_embeddings)
            meta_matrix = np.vstack(final_meta)

            # Save
            np.savetxt(out_embeddings, embeddings_matrix)
            np.savetxt(out_metadata, meta_matrix)
            with open(out_genenames, "w") as f:
                for gname in final_gene_names:
                    f.write(gname + "\n")

            print(f"  Saved embeddings to:  {out_embeddings}")
            print(f"  Saved metadata to:    {out_metadata}")
            print(f"  Saved gene names to:  {out_genenames}")

        except Exception as e:
            print(f"  Error processing {base}: {e}")

    print("\nDone! All paired FASTA/GFF sets processed.")


if __name__ == "__main__":
    main()
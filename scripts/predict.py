"""
predict.py

Usage:
    python predict.py \
      --input_dir /path/to/input \
      --model_dir /path/to/model \
      --output_dir /path/to/output

Description:
  1) Looks in `--input_dir` for triples of files:
      - {NAME}_embeddings.txt
      - {NAME}_metadata.txt
      - {NAME}_genename.txt
  2) Loads each triple, concatenates embeddings + metadata (horizontal).
  3) Loads a single saved model from `--model_dir` (model.pth).
  4) Makes a prediction for each row in the input matrix.
  5) Saves results to `--output_dir/{NAME}_predictions.csv`, with two columns:
      Gene Name, Prediction
"""

import os
import sys
import argparse
import glob
import numpy as np
import pandas as pd
import torch
from torch.utils.data import TensorDataset, DataLoader

from model import (
    PositionalEncoding,
    TransformerModel,
    CombinedModel
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir",  required=True, help="Folder with embeddings, metadata, and genename files.")
    parser.add_argument("-m","--model_dir",  required=True, help="Folder containing the saved model file (model.pth).")
    parser.add_argument("-o","--output_dir", required=True, help="Folder to save the prediction CSV files.")
    args = parser.parse_args()

    input_dir  = args.input_dir
    model_dir  = args.model_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1) Construct or load your model
    #    Adjust model1_params, combined_dim, output_dim if needed
    # -------------------------------------------------------------------------
    model1_params = {
        'input_dim': 1301,   # e.g. (embedding_dim + meta_dim)
        'num_heads': 4,
        'num_layers': 1,
        'dim_feedforward': 128,
        'dropout': 0.0
    }
    combined_dim = 256
    output_dim = 1

    # Initialize model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)
    model = CombinedModel(model1_params, combined_dim, output_dim).to(device)
    model.eval()

    # -------------------------------------------------------------------------
    # 2) Load the saved model weights (assumes "model.pth" in model_dir)
    # -------------------------------------------------------------------------
    model_path = model_dir
    if not os.path.exists(model_path):
        print(f"Error: Could not find saved model at {model_path}")
        sys.exit(1)

    print(f"Loading model weights from: {model_path}")
    model.load_state_dict(torch.load(model_path, map_location=device))

    # -------------------------------------------------------------------------
    # 3) Find matching file sets in input_dir
    #    We'll assume a pattern: {BASE}_embeddings.txt, {BASE}_metadata.txt, {BASE}_genename.txt
    # -------------------------------------------------------------------------
    embeddings_files = sorted(glob.glob(os.path.join(input_dir, "*_embeddings.txt")))

    # For each embeddings file, check if we have matching metadata and genename
    for emb_file in embeddings_files:
        base_name = emb_file.rsplit("_embeddings.txt", 1)[0]
        meta_file = base_name + "_metadata.txt"
        gname_file = base_name + "_genename.txt"

        # Check existence
        if not all([os.path.exists(meta_file), os.path.exists(gname_file)]):
            print(f"Skipping {emb_file}, missing {meta_file} or {gname_file}.")
            continue

        # Extract just the filename (or subfolder) for final output
        out_name = os.path.basename(base_name)  # e.g. "GCF_000011065.1"

        print(f"\nProcessing set: {out_name}")
        print(f"  Embeddings: {emb_file}")
        print(f"  Metadata:   {meta_file}")
        print(f"  GeneNames:  {gname_file}")

        # ---------------------------------------------------------------------
        # 4) Load data
        # ---------------------------------------------------------------------
        # Embeddings
        try:
            X_embedding = np.loadtxt(emb_file)
            # Metadata
            X_meta = np.loadtxt(meta_file)
            # Gene names
            with open(gname_file, "r") as f:
                gene_names = [line.strip() for line in f]
        except Exception as e:
            print(f"  Error loading files for {out_name}: {e}")
            continue

        # Check dimensions
        if X_embedding.shape[0] != X_meta.shape[0]:
            print(f"  Mismatch: Embeddings rows={X_embedding.shape[0]}, Meta rows={X_meta.shape[0]}. Skipping.")
            continue
        if len(gene_names) != X_embedding.shape[0]:
            print(f"  Mismatch: {len(gene_names)} gene names, but {X_embedding.shape[0]} data rows. Skipping.")
            continue

        # Combine horizontally
        X_combined = np.hstack([X_embedding, X_meta])
        print(f"  Combined shape = {X_combined.shape}")

        # ---------------------------------------------------------------------
        # 5) Model Inference
        # ---------------------------------------------------------------------
        X_tensor = torch.FloatTensor(X_combined)
        test_dataset = TensorDataset(X_tensor)
        test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=32, shuffle=False)

        predictions = []
        model.eval()
        with torch.no_grad():
            for (batch_x,) in test_loader:
                batch_x = batch_x.to(device)
                batch_out = model(batch_x)
                predictions.append(batch_out.cpu().numpy())

        predictions = np.concatenate(predictions, axis=0).flatten()

        # ---------------------------------------------------------------------
        # 6) Save output as a CSV: gene names + predictions
        # ---------------------------------------------------------------------
        result_df = pd.DataFrame({
            "Gene_Name": gene_names,
            "Prediction": predictions
        })
        out_csv = os.path.join(output_dir, f"{out_name}_predictions.csv")
        result_df.to_csv(out_csv, index=False)
        print(f"  Predictions saved to: {out_csv}")

    print("\nAll data processed.")


if __name__ == "__main__":
    main()

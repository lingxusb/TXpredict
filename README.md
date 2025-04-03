# TXpredict：predicting microbial transcriptome using genome sequence
![github_v2](https://github.com/user-attachments/assets/89c81779-19f0-4184-a526-36ca36188abf)


We present TXpredict, a transcriptome prediction tool that generalizes to novel microbial genomes. By leveraging information learned from a large protein language model (ESM2), TXpredict achieves an average Spearman correlation of 0.53 and 0.62 in predicting gene expression for new bacterial and fungal genomes. We further extend this framework to predict transcriptomes for 2,685 additional microbial genomes spanning 1,744 genera, a large proportion of which remain uncharacterized at the transcriptional level. Additionally, TXpredict enables the prediction of condition-specific gene expression, providing a powerful tool for understanding microbial adaptation and facilitating rational design of gene regulatory sequences.

## Table of Contents

- [Installation](#Installation)
- [Command lines](#Command-lines)
- [Jupyter notebooks](#Juypter-notebooks)
- [Colab notebooks](#Colab-notebooks)
- [Trained models](#Trained-models)
- [Acknowledgement](#Acknowledgement)

## Installation
Python package dependencies:
- torch 2.0.1
- pandas 2.2.0
- seaborn 0.13.2
- biopython 1.81

We recommend using [Conda](https://docs.conda.io/en/latest/index.html) to install our packages. For convenience, we have provided a conda environment file with package versions that are compatiable with the current version of the program. The conda environment can be setup with the following comments:

1. Clone this repository:
   ```bash
     git clone https://github.com/lingxusb/TXpredict.git
     cd TXpredict
   ```

2. Create and activate the Conda environment:
   ```bash
   conda env create -f env.yml
   conda activate TXpredict
   ```
## Command lines
### Embedding generation
The `embeddings.py` script:
1. Traverses subfolders looking for `.fna`/`.fasta` + `.gff`/`.gtf` pairs with matching base names
2. Extracts protein sequences from these pairs
3. Computes ESM-2 embeddings for proteins below the specified length threshold
4. Generates metadata (normalized length + amino acid proportions)
5. Saves the results while preserving the input folder structure

#### Usage

```bash
python embeddings.py -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY [-l LENGTH_THRESHOLD]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input_dir` | Path to input directory containing subfolders with FASTA/GFF pairs (required) |
| `-o`, `--output_dir` | Path where output files will be saved (required) |
| `-l`, `--length_threshold` | Maximum protein length to include (default: 1500) |

#### Output Files

For each processed FASTA/GFF pair, the script generates:

- `*_embeddings.txt`: 2D array of protein embeddings
- `*_metadata.txt`: 2D array of metadata (normalized length + 20 AA proportions)
- `*_genename.txt`: 1D list of gene names

### Model prediction

The `predict.py` script:
1. Finds matching sets of embedding, metadata, and gene name files in the input directory
2. Loads and combines the data (embeddings + metadata)
3. Applies a pre-trained model to make predictions
4. Saves the results as CSV files with gene names and their corresponding predictions

#### Usage

```bash
python predict.py --input_dir INPUT_DIRECTORY --model_dir MODEL_DIRECTORY --output_dir OUTPUT_DIRECTORY
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input_dir` | Directory containing the input files (required) |
| `-m`, `--model_dir` | Directory containing the trained model file (required) |
| `-o`, `--output_dir` | Directory where prediction results will be saved (required) |

#### Input Files

The script expects sets of three files with matching names in the input directory:
- `{NAME}_embeddings.txt`: File containing protein embeddings
- `{NAME}_metadata.txt`: File containing metadata features
- `{NAME}_genename.txt`: File containing gene names

#### Output

For each processed set of input files, the script generates:
- `{NAME}_predictions.csv`: A CSV file with two columns:
  - `Gene_Name`: The name of the gene
  - `Prediction`: The model's prediction value
 
## Jupyter notebooks
### Data preprocessing
`01_data_preprocessing.ipynb` handles the preprocessing steps needed to prepare data for the TXpredict model.
The notebook performs three main tasks:

1. **Calculating normalized gene expression** - Processes RNA-seq count data to compute TPM (Transcripts Per Million) values with z-score normalization.
2. **Generating ESM embeddings** - Uses the ESM-2 model to create protein embeddings from genome annotation files.
3. **Preparing training data** - Combines embedding data with expression data to create the final training dataset.

The notebook produces several output files including:
- `{strain}_filtered_log_tpm.csv` - Normalized expression data
- `{strain}_filtered_embeddings.txt` - Protein embeddings for model input
- `{strain}_filtered_meta.txt` - Metadata for each protein

### Model training
`02_model_training.ipynb` demonstrates how to train the TXpredict model for gene expression prediction. The notebook covers the complete workflow:

1. **Data loading** - Loads preprocessed embeddings and metadata from the previous preprocessing step
2. **Model definition** - Implements the model to learn from sequence embeddings
3. **Training process** - Trains the model with evaluation metrics
4. **Model saving** - Saves the trained model for later use in prediction




## Colab notebooks
We have provided [Colab notebooks](https://colab.research.google.com/drive/1Kd-QIwTgESIg_62b4rstuT1KO-NMqtPL?usp=sharing) for transcriptome prediction in the web browser. Please also check our [Colab instruction](https://github.com/lingxusb/TXpredict/blob/main/Colab_instruction.md). We also provided a [Colab notebook](https://colab.research.google.com/drive/1xvgQlRsz8vUMW_R7MZgNlHvIh5hEyqhN?usp=sharing) for fungal transcriptome prediction.
- The only required inputs are genome sequence file (.fna or .fasta) and the annotation file (.gtf, .gff or .gff3). Please check our [example data](https://github.com/lingxusb/TXpredict/tree/main/example_data)
- Please connect to a GPU instance (e.g. T4, Runtime -> Change runtime type -> T4 GPU).
- It takes ~20min to predict transcriptome for a genome with 4k genes.

## Trained models
Our transcriptome prediciton models are available from [Huggingface](https://huggingface.co/lingxusb/TXpredict/tree/main).

## Acknowledgement
We deeply appreciate the [experimental works and datasets](https://github.com/lingxusb/TXpredict/blob/main/Acknowledgement.md) that make our work possible.

## References
- ESM2: https://github.com/facebookresearch/esm
- [Predicting microbial transcriptome using genome sequence](https://www.biorxiv.org/content/10.1101/2024.12.30.630741v1)

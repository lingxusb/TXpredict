import re
from Bio import SeqIO
from collections import Counter, OrderedDict

def get_protein_sequences(genome_file, annotation_file):
    """
    Reads a FASTA genome_file and a GFF annotation_file, extracts all CDS regions,
    and translates them to protein sequences. Returns:
      - protein_sequences         : list of Bio.Seq objects (proteins)
      - gene_names               : list of gene identifiers
      - gene_lengths             : list of integer lengths for each protein
      - amino_acid_proportions   : list of 20-element vectors with AA proportions
    """
    # Load genome sequences
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    # Define the 20 standard amino acids
    standard_aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_template = OrderedDict.fromkeys(standard_aas, 0)

    protein_sequences = []
    gene_names = []
    gene_lengths = []
    amino_acid_proportions = []

    with open(annotation_file, "r") as gff:
        for line in gff:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Skip lines that don't have enough fields

            feature_type = parts[2].lower()
            if feature_type == "cds":
                seqid, source, feature, start, end, score, strand, phase, attributes = parts

                # Example: capturing a gene name from 'gene-...' in attributes
                match = re.search(r"gene-([^;]+)", attributes)
                if match:
                    gene = match.group(1)
                else:
                    gene = ""

                # Convert to 0-based indexing
                start, end = int(start) - 1, int(end)
                if seqid in genome:
                    # Extract DNA sequence
                    dna_seq = genome[seqid].seq[start:end]
                    if strand == '-':
                        dna_seq = dna_seq.reverse_complement()

                    # Translate to protein sequence (stop at first stop codon)
                    protein_seq = dna_seq.translate(to_stop=True)

                    # Filter out invalid or very large proteins
                    if 0 < len(protein_seq) < 1500:
                        protein_sequences.append(protein_seq)
                        gene_names.append(gene)
                        gene_lengths.append(len(protein_seq))

                        # Calculate amino acid proportions
                        aa_counts = Counter(protein_seq)
                        total_aa = sum(aa_counts.values())
                        aa_props = aa_template.copy()
                        for aa in aa_props:
                            aa_props[aa] = aa_counts.get(aa, 0) / total_aa
                        amino_acid_proportions.append(list(aa_props.values()))

    return protein_sequences, gene_names, gene_lengths, amino_acid_proportions

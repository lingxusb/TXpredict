# TXpredict：predicting microbial transcriptome using genome sequence
![github](https://github.com/user-attachments/assets/697aeda2-d6d4-421d-8240-2368c4570c65)

We present TXpredict, a transcriptome prediction tool that generalizes to novel microbial genomes. By leveraging information learned from a large protein language model (ESM2), TXpredict achieves an average Spearman correlation of 0.53 in predicting gene expressions for new bacterial genomes. We further extend this framework to predict transcriptomes for 900 additional microbial genomes spanning 280 genera, a large proportion of which remain uncharacterized at the transcriptional level. Additionally, TXpredict enables the prediction of condition-specific gene expression, providing a powerful tool for understanding microbial adaptation and facilitating rational design of gene regulatory sequences.

### Models
Our transcriptome prediciton models are available from [Huggingface](https://huggingface.co/lingxusb/TXpredict/tree/main).

### Colab notebooks
We have provided [Colab notebooks](https://colab.research.google.com/drive/1Kd-QIwTgESIg_62b4rstuT1KO-NMqtPL?usp=sharing) for transcriptome prediction in the web browser. Please also check our [Colab instruction](https://github.com/lingxusb/TXpredict/blob/main/Colab_instruction.md)
- The only required inputs are genome sequence file (.fna or .fasta) and the annotation file (.gtf, .gff or .gff3). Please check our [example data](https://github.com/lingxusb/TXpredict/tree/main/example_data)
- Please connect to a GPU instance (e.g. T4, Runtime -> Change runtime type -> T4 GPU).
- It takes ~20min to predict transcriptome for a genome with 4k genes.

### References
- ESM2: https://github.com/facebookresearch/esm
- Predicting microbial transcriptome using genome sequence

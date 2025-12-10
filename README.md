# BIOL4300 Project: Phylogenetics of Human Herpesvirus 6B

Scripts and data necessary to complete gene-specific analysis of Human Herpesvirus 6B (HHV-6B).

## Folder organization
1. `datasets/code`: All scripts used (all R scripts ran in R version 4.5.1, all python scripts ran in Python 3.11.9)
2. `datasets/input_data`: All aligned, trimmed seqeunce data input into IQ-TREE v3.0.1 (`https://iqtree.github.io/`) and BEAST v2.7.7 (`https://www.beast2.org/`)
3. `datasets/input_data`: Raw data from IQ-TREE and BEAST
4. `datasets/trees`: Final trees visualizations (in pdf files) and tree files (.treefile from IQ-TREE and .tree from BEAST)

## Steps
1. Get raw FASTA gene-specific data from Genbank (`datasets/code/01.get_gene_fastas.R`).
2. Get whole genome data from GenBank (`datasets/code/02.get_complete_genomes.R`).
3. Extract specific genes from whole genome data, add gene sequences to original gene FASTA files (`datasets/code/03.extract_genes.py`).
4. Get alignments by inputting gene FASTA files into MAFFT version 7 (`https://mafft.cbrc.jp/alignment/software/`)
5. Trim Alignments (`datasets/code/05.trim.R`).
6. Construct Trees:
   - Maximum Likelihood trees (IQ-TREE v3.0.1)
   - Bayesian trees (BEAST v2.7.7)
   - tree files are output to `datasets/trees`
7. Tree visualization (R [`datasets/code/07.plot_Bayesian.R`, `datasets/code/07.plot_ML.R`] and FigTree v1.4.4 (`https://tree.bio.ed.ac.uk/software/figtree/`)

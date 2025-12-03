# BIOL4300 Project: Phylogenetics of Human Herpesvirus 6B

Scripts and data necessary to complete gene-specific analysis of Human Herpesvirus 6B (HHV-6B).

## Folder organization
1. code: All scripts used (all R scripts ran in R version 4.5.1, all python scripts ran in Python 3.11.9)
2. input_data: All aligned, trimmed seqeunce data input into IQ-TREE v. 3.0.1 and BEAST v. 2.7.7
3. output_data: Raw data from IQ-TREE and BEAST
4. trees: Final trees visualizations (in pdf files) and tree files (.treefile from IQ-TREE and .tree from BEAST)

## Steps
1. Get raw FASTA gene-specific data from Genbank.
2. Get whole genome data from GenBank.
3. Extract specific genes from whole genome data, add gene sequences to original gene FASTA files.
4. Get alignments by input gene FASTA files into MAFFT version 7 (`https://mafft.cbrc.jp/alignment/software/`)
5. Trimming Alignments
7. Trees:
   - ML trees
   - Bayesian trees in BEAST
   - tree files are output to `./trees`
8. Tree visualization (R and FigTree)

## Currently under construction :)

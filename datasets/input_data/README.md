# Data Folder Organization

| Folder                    | Description                                                                                                                                                           |
| :------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`01.raw_data/`**        | Contains the raw FASTA sequence files and GenBank files downloaded from NCBI using the scripts `code/01.get_gene_fastas.R` and `code/02.get_complete_genomes.R`                    |
| **`02.mafft_input/`**     | Raw FASTA files input into MAFFT to create Multiple sequence alignments (MSAs). Created by script `code/03.extract_genes.py`                                |
| **`alignments/`**         | Cleaned and trimmed versions of the MAFFT alignments. Gaps and poorly aligned regions were filtered out using script `code/05.trim.R`. These files were inputs for both maximum likelihood and Bayesian tree inference.                                           |

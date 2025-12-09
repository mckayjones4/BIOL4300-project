# Import required libraries
# Bio.SeqIO is part of Biopython, a toolkit for working with biological sequence data.
# os is a standard Python library for interacting with files and directories.
from Bio import SeqIO
import os

# 1. Set working directory
# This is the folder where both GenBank (.gb) files and FASTA outputs are stored.
# Change this path as needed.
gb_dir = "C:/Users/jones/OneDrive/Documents/BIOL4300-project/datasets/input_data/1.raw_data"

# 2. Define the genes of interest
# These are the HHV-6B genes you want to extract from each genome file.
# Add or remove gene names depending on project.
genes_of_interest = ['U39', 'U47', 'U69', 'U83', 'U90']

# 3. Loop through every file in the directory
# os.listdir(gb_dir) gives you a list of all files in the folder.
# We check each file to see if it ends with ".gb" (GenBank format).
for gb_file in os.listdir(gb_dir):
    if gb_file.endswith(".gb"):
        # Build the full file path
        gb_path = os.path.join(gb_dir, gb_file)

        # Read the GenBank file into a SeqRecord object
        # "genbank" tells Biopython the file format.
        gb_obj = SeqIO.read(gb_path, "genbank")

        # Extract the accession number
        accession = gb_obj.id

        # 4. Look through all features in the GenBank file
        # A GenBank file contains "features" such as genes, CDS, etc.
        for feature in gb_obj.features:
            # Check if the feature is a gene and has a "gene" qualifier
            if feature.type == "gene" and "gene" in feature.qualifiers:
                # Get the gene name from the qualifiers
                gene_name = feature.qualifiers["gene"][0]

                # Only process if the gene is in list of interest
                if gene_name in genes_of_interest:
                    # Extract the sequence corresponding to this gene
                    extract = feature.extract(gb_obj)

                    # 5. Clean up the FASTA header
                    # Format: HHV6_<Gene>_<Accession>
                    extract.id = f"HHV6_{gene_name}_{accession}"
                    extract.name = ""          # remove extra name info
                    extract.description = ""   # remove description text

                    # 6. Append to the gene-specific FASTA file
                    # Each gene gets its own FASTA file (e.g., hhv_6b_U39.fasta).
                    # We open the file in "a" (append) mode, so new sequences are added without overwriting existing ones.
                    out_fasta = os.path.join(gb_dir, f"hhv_6b_{gene_name}.fasta")
                    with open(out_fasta, "a") as handle:
                        SeqIO.write(extract, handle, "fasta")

                    print(f"Appended {gene_name} from {gb_file} to {out_fasta}")

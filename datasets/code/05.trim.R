############################################################
# HHV-6B Gene Sequence Trimming
# ----------------------------------------------------------
# This script:
# 1. Reads in a FASTA file of gene sequences
# 2. Removes duplicate sequences
# 3. Trims alignment columns with too many gaps
# 4. Writes out a cleaned FASTA file
# 5. Ensures only unique sequences remain
############################################################

# Load required packages
library(evobiR)
library(ape)
library(msaR)
library(phangorn)
library(Biostrings) 

# Step 1: Read in FASTA file
# Replace "hhv6b_U90.fasta" with your input file as needed
u90_sequences <- read.FASTA("hhv6b_U90.fasta")

# Inspect how many sequences were read
length(u90_sequences)

# Check how many unique sequence names (headers) are present
length(unique(names(u90_sequences)))

# Visualize multiple sequence alignment
msaR(u90_sequences)

# Step 2: Remove duplicate sequences by name
u90_unique <- u90_sequences[!duplicated(names(u90_sequences))]

# Convert to phyDat object (needed for trimming)
u90_phyDat <- phyDat(u90_unique, type = "DNA")

# Step 3: Trim alignment columns with >50% gaps
u90_trimmed <- u90_phyDat[, colMeans(as.character(u90_phyDat) == "-") < 0.5]

# Convert back to DNAbin format for writing/inspection
u90_trimmed_bin <- as.DNAbin(u90_trimmed)

# Visualize trimmed alignment
msaR(u90_trimmed_bin)

# Write trimmed FASTA file
write.FASTA(u90_trimmed_bin, "U90_trimmed.fasta")

# Step 4: Ensure only unique sequences remain
# Read the trimmed FASTA back in as a DNAStringSet
u90_trimmed_set <- readDNAStringSet("U90_trimmed.fasta")

# Count how many sequences
length(u90_trimmed_set)

# Remove duplicates (based on sequence content) (289 for gene U90)
u90_final_unique <- unique(u90_trimmed_set)

# Count how many unique sequences remain (27 for gene U90)
length(u90_final_unique)

# Visualize final alignment
msaR(u90_final_unique)

# Write final unique FASTA file
writeXStringSet(u90_final_unique, filepath = "U90_unique.fasta") #these files are input into BEAST and IQ-TREE

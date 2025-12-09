############################################################
# 1. HHV-6B Gene Sequence Downloader and Header Cleaner
# ----------------------------------------------------------
# This script:
# 1. Searches NCBI GenBank for HHV-6B gene sequences
# 2. Fetches sequences in FASTA format
# 3. Cleans FASTA headers to a standard format
# 4. Writes cleaned FASTA files for each gene of interest
############################################################

# Load required package
library(rentrez)

# Define a reusable header cleaning function
# This function takes a FASTA header line and returns a cleaned version:
# Format: >HHV6_<Gene>_<Accession>
clean_header <- function(header_line, gene_prefix = "HHV6", gene = "U90") {
  # Extract accession (first token after ">")
  acc <- sub("^>(\\S+).*", "\\1", header_line)
  paste0(">", gene_prefix, "_", gene, "_", acc)
}

# ----------------------------------------------------------
# HHV-6B U90 gene (with date filtering)
u90_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',   # Organism: HHV-6B
                              'U90[ALL]',          # Gene of interest
                              '2019:2025[PDAT]',   # Publication date range
                              'partial cds[Title]',# Restrict to partial CDS
                              sep = " AND "
                            ),
                            retmax = 1000)

cat("Number of U90 candidates found:", u90_search$count, "\n")

# Inspect search results
u90_summary <- entrez_summary(db = "nuccore", id = u90_search$ids)
titles <- sapply(u90_summary, `[[`, "title")
tail(titles, 10)

# Fetch sequences in FASTA format
u90_output <- entrez_fetch(db = "nuccore", id = u90_search$ids, rettype = "fasta")
u90_fasta_raw <- strsplit(u90_output, "\n")[[1]]

# Clean headers
u90_out_lines <- ifelse(startsWith(u90_fasta_raw, ">"),
                        vapply(u90_fasta_raw, clean_header, character(1),
                               gene_prefix = "HHV6", gene = "U90"),
                        u90_fasta_raw)

# Write cleaned FASTA
writeLines(u90_out_lines, "hhv6b_U90.fasta")

# ----------------------------------------------------------
# HHV-6B U83 gene (no date filtering)
u83_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U83[ALL]',
                              'cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 1000)

cat("Number of U83 candidates found:", u83_search$count, "\n")

u83_summary <- entrez_summary(db = "nuccore", id = u83_search$ids)
titles <- sapply(u83_summary, `[[`, "title")
tail(titles, 10)

u83_output <- entrez_fetch(db = "nuccore", id = u83_search$ids, rettype = "fasta")
u83_fasta_raw <- strsplit(u83_output, "\n")[[1]]

u83_out_lines <- ifelse(startsWith(u83_fasta_raw, ">"),
                        vapply(u83_fasta_raw, clean_header, character(1),
                               gene_prefix = "HHV6", gene = "U83"),
                        u83_fasta_raw)

writeLines(u83_out_lines, "hhv6b_U83.fasta")

# ----------------------------------------------------------
# HHV-6B U69 gene (no date filtering)
u69_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U69[ALL]',
                              'cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 50)

cat("Number of U69 candidates found:", u69_search$count, "\n")

u69_summary <- entrez_summary(db = "nuccore", id = u69_search$ids)
titles <- sapply(u69_summary, `[[`, "title")
tail(titles, 10)

u69_output <- entrez_fetch(db = "nuccore", id = u69_search$ids, rettype = "fasta")
u69_fasta_raw <- strsplit(u69_output, "\n")[[1]]

u69_out_lines <- ifelse(startsWith(u69_fasta_raw, ">"),
                        vapply(u69_fasta_raw, clean_header, character(1),
                               gene_prefix = "HHV6", gene = "U69"),
                        u69_fasta_raw)

writeLines(u69_out_lines, "hhv6b_U69.fasta")

# ----------------------------------------------------------
# HHV-6B U47 gene (no date filtering)
u47_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U47[Title]',
                              'cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 50)

cat("Number of U47 candidates found:", u47_search$count, "\n")

u47_summary <- entrez_summary(db = "nuccore", id = u47_search$ids)
titles <- sapply(u47_summary, `[[`, "title")
head(titles, 10)

u47_output <- entrez_fetch(db = "nuccore", id = u47_search$ids, rettype = "fasta")
u47_fasta_raw <- strsplit(u47_output, "\n")[[1]]

u47_out_lines <- ifelse(startsWith(u47_fasta_raw, ">"),
                        vapply(u47_fasta_raw, clean_header, character(1),
                               gene_prefix = "HHV6", gene = "U47"),
                        u47_fasta_raw)

writeLines(u47_out_lines, "hhv6b_U47.fasta")

# ----------------------------------------------------------
# HHV-6B U39 gene (no date filtering)
u39_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U39[ALL]',
                              'complete cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 50)

cat("Number of U39 candidates found:", u39_search$count, "\n")

u39_summary <- entrez_summary(db = "nuccore", id = u39_search$ids)
titles <- sapply(u39_summary, `[[`, "title")
head(titles, 10)

u39_output <- entrez_fetch(db = "nuccore", id = u39_search$ids, rettype = "fasta")
u39_fasta_raw <- strsplit(u39_output, "\n")[[1]]

u39_out_lines <- ifelse(startsWith(u39_fasta_raw, ">"),
                        vapply(u39_fasta_raw, clean_header, character(1),
                               gene_prefix = "HHV6", gene = "U39"),
                        u39_fasta_raw)

writeLines(u39_out_lines, "hhv6b_U39.fasta")


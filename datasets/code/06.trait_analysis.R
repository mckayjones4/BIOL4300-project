############################################################
# HHV-6B U39 Metadata Extraction for Trait Analysis
# ----------------------------------------------------------
#   This script reads a FASTA file of HHV-6B U39 sequences,
#   extracts metadata from GenBank records (collection date,
#   geographic location, host, etc.), computes GC content,
#   and merges everything into a single CSV file for use in
#   downstream visualization.
#
# Inputs:
#   - FASTA file containing unique sequences (e.g., "Gene_u39_unique.fasta")
#
# Outputs:
#   - CSV file ("genbank_u39_metadata.csv") with metadata and GC content
############################################################

# Load required packages
library(Biostrings)  # for reading FASTA and computing GC content
library(rentrez)     # for fetching GenBank records
library(stringr)     # for regex-based text extraction

# Step 1: Read FASTA file
# This FASTA file should contain the sequences used in IQ-TREE/BEAST2
fasta <- readDNAStringSet("Gene_u39_unique.fasta")

# Extract sequence names (headers)
full_names <- names(fasta)  

# Extract GenBank accession IDs from headers
# Assumes headers are formatted like HHV6_U39_<Accession>
ids <- sub(".*_", "", full_names)

# Step 2: Define metadata extraction function
# Given an accession ID, fetch the GenBank record and parse relevant metadata fields (using regex).
get_metadata <- function(acc, full_name) {
  gb_record <- entrez_fetch(db = "nucleotide", id = acc,
                            rettype = "gb", retmode = "text")
  
  # Extract collection date
  date <- str_match(gb_record, "/collection_date=\"([^\"]+)\"")[,2]
  
  # Extract geographic location (country or geo_loc_name)
  geo <- str_match(gb_record, "/country=\"([^\"]+)\"")[,2]
  if (is.na(geo)) {
    geo <- str_match(gb_record, "/geo_loc_name=\"([^\"]+)\"")[,2]
  }
  
  # Extract isolation source
  source <- str_match(gb_record, "/isolation_source=\"([^\"]+)\"")[,2]
  
  # Extract sequence length (bp)
  bp <- str_match(gb_record, "LOCUS\\s+\\S+\\s+(\\d+) bp")[,2]
  
  # Extract host
  host <- str_match(gb_record, "/host=\"([^\"]+)\"")[,2]
  
  # Extract protein ID
  protein <- str_match(gb_record, "/protein_id=\"([^\"]+)\"")[,2]
  
  # Return metadata as a data frame row
  return(data.frame(
    TaxaName = full_name,
    Accession = acc,
    CollectionDate = date,
    GeoLocation = geo,
    IsolationSource = source,
    Host = host,
    LengthBP = bp,
    Protein = protein,
    stringsAsFactors = FALSE
  ))
}

# Step 3: Apply metadata function to all sequences
# lapply() runs get_metadata() for each accession ID
meta_list <- lapply(seq_along(ids), function(i) get_metadata(ids[i], full_names[i]))

# Combine list of data frames into one table
meta_table <- do.call(rbind, meta_list)

# Step 4: Compute GC content from FASTA sequences
# letterFrequency() counts bases; as.prob=TRUE gives proportions
gc_content <- letterFrequency(fasta, letters = c("G", "C"), as.prob = TRUE)
gc_percent <- rowSums(gc_content) * 100  # convert to percentage

# Create GC content data frame
gc_df <- data.frame(TaxaName = full_names, GCpercent = gc_percent)

# Step 5: Merge GC content into metadata table
meta_table <- merge(meta_table, gc_df, by = "TaxaName")

# Step 6: Save metadata table to CSV
write.csv(meta_table, "genbank_u39_metadata.csv", row.names = FALSE) # change name for each gene


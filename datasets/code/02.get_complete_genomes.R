############################################################
# HHV-6B Complete Genome Downloader
# ----------------------------------------------------------

#   This script uses the NCBI Entrez utilities (via the rentrez
#   package) to search for and download complete HHV-6B genomes
#   from GenBank. It filters out iciHHV-6B genomes and saves each
#   genome as an individual GenBank (.gb) file.
############################################################

#load rentrez package
library(rentrez)

# Define output directory for downloaded GenBank files
out_dir <- "C:/Users/jones/OneDrive/Documents/BIOL4300-project/datasets/input_data/1.raw_data" #change as needed
# Create directory
dir.create(out_dir, showWarnings = FALSE)

# ----------------------------------------------------------
# Define Search parameters
# - TaxID for HHV-6B: 32604
# - Restrict to complete genomes
# ----------------------------------------------------------
search_term <- paste(
  "txid32604[ORGN]",        # Organism: HHV-6B
  "complete genome[Title]", # Only complete genomes
  sep = " AND "
)

# Perform search using entrez_search() function
genome_search <- entrez_search(db = "nuccore",
                               term = search_term,
                               retmax = 500)

#we found 7 complete genomes
cat("Number of complete HHV-6B genomes found:", genome_search$count, "\n")

# Summarize search results
genome_summary <- entrez_summary(db = "nuccore", id = genome_search$ids)
titles <- sapply(genome_summary, `[[`, "title")

#Visually inspect titles (optional)
head(titles, 10)

# Remove iciHHV-6B genomes (confounds our data, this is optional)
keep <- !grepl("ici", titles, ignore.case = TRUE)
filtered_ids <- genome_search$ids[keep]

data <- data.frame(
  id = filtered_ids,
  accession_num = sapply(genome_summary[keep], function(x) x$accessionversion),
  stringsAsFactors = FALSE
)

# Visually inspect output
print(data)

# Download each genome as a GenBank (.gb) file
for (i in seq_len(nrow(data))) {
  acc <- data$accession_num[i]
  gb_text <- entrez_fetch(db = "nuccore", id = data$id[i],
                          rettype = "gb", retmode = "text")
  gb_file <- file.path(out_dir, paste0(acc, ".gb"))
  writeLines(gb_text, gb_file)
  cat("Saved:", gb_file, "\n")
}

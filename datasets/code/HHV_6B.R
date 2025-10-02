#load releveant packages
library(rentrez)
library(XML)
library(dplyr)
library(httr)
library(seqinr)

#We first search for gene isolates, excluding partial and complete genome results
#taxid for HHV-6B: 32604

#HHV-6B U90 gene
u90_search <- entrez_search(db = "nuccore",
                            term = paste(
                            'txid32604[ORGN]',
                            'U90[ALL]',
                            '2019:2025[PDAT]',
                            'partial cds[Title]',
                            sep = " AND "
                            ),
                            retmax = 1000)

#How many candidates we found: 281
u90_search$count

sum_u90_search <- entrez_summary(db = "nuccore", id=u90_search$ids)
sum_u90_search$ids

#manually checking the titles
titles <- sapply(sum_u90_search, `[[`, "title")
tail(titles, 100)

#writing to a FASTA file
output <- entrez_fetch(db="nuccore", id=u90_search$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="hhv_6b_u90.fasta")

#HHV-6B U86 gene
u86_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U86[ALL]',
                              '2019:2025[PDAT]',
                              'partial cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 1000)

#How many candidates we found: 49
u86_search$count

sum_u86_search <- entrez_summary(db = "nuccore", id=u86_search$ids)
sum_u86_search$ids

#manually checking the titles
titles <- sapply(sum_u86_search, `[[`, "title")
tail(titles, 100)

#writing to a FASTA file
output <- entrez_fetch(db="nuccore", id=u86_search$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="hhv_6b_u86.fasta")


#HHV-6B U39 gene, no date filtering because Granger's paper didn't conduct a U39 (glycoprotein B) gene analysis 
u39_search <- entrez_search(db = "nuccore",
                            term = paste(
                              'txid32604[ORGN]',
                              'U39[ALL]',
                              'complete cds[Title]',
                              sep = " AND "
                            ),
                            retmax = 50)

#How many candidates we found: 281
u39_search$count

sum_u39_search <- entrez_summary(db = "nuccore", id=u39_search$ids)
sum_u39_search$ids

#manually checking the titles
titles <- sapply(sum_u39_search, `[[`, "title")
tail(titles, 100)

#writing to a FASTA file
output <- entrez_fetch(db="nuccore", id=u39_search$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="hhv_6b_u39.fasta")


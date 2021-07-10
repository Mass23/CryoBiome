#!/usr/bin/Rscript

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

# Running gRodon on MAGs
suppressMessages(library(Biostrings))
suppressMessages(library(coRdon))
suppressMessages(library(matrixStats))
library(gRodon)

# Load your *.ffn file into R
genes <- readSet(file=snakemake@input[["FNA"]])

# Storing codon counts for each sequence
gene_table <- codonTable(genes)

# Estimating codon counts
cc <- codonCounts(gene_table)

# Writing the output to file
write.table(cc, file = snakemake@output[["COUNT"]], sep="\t", row.names=FALSE, quote=FALSE)

#!/usr/bin/Rscript

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

# Running gRodon on MAGs
suppressMessages(library(Biostrings))
suppressMessages(library(coRdon))
suppressMessages(library(matrixStats))
library(gRodon)

# Load your *.ffn file into R
genes <- readDNAStringSet(snakemake@input[["FFN"]])

# Subset your sequences to those that code for proteins
CDS_IDs <- readLines(snakemake@input[["CDS"]])
gene_IDs <- gsub(" .*","",names(genes)) #Just look at first part of name before the space
genes <- genes[gene_IDs %in% CDS_IDs]

#Search for genes annotated as ribosomal proteins
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)

# Since some MAGs are not very complete the Growth Prediction is run using "tryCatch"
# Example usage: https://statisticsglobe.com/using-trycatch-function-to-handle-errors-and-warnings-in-r
# Running growth prediction
pred_growth <- tryCatch({
    print("Running growth prediction")
    result <- predictGrowth(genes, highly_expressed, mode="partial")
}, error = function(e) {
    print("Creating empty file if errors are thrown")
    result <- NULL
})

# pred_growth <- predictGrowth(genes, highly_expressed, mode="partial")

# Writing the output to file
if(!is.null(pred_growth)){ write.table(pred_growth, file = snakemake@output[["PRED"]], sep="\t", row.names=FALSE, quote=FALSE) } else { write.table(pred_growth, file = snakemake@output[["PRED"]], sep="\t", row.names=FALSE, quote=FALSE) }

#!/usr/bin/Rscript

# Script to merge outputs from gRodon for all MAGs

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

all_paths <-
  list.files(path = dirname(snakemake@input[["PRED"]]),
             pattern = "*.txt",
             full.names = TRUE)

# Removing EMPTY files from the list             
all_paths <- all_paths[which(file.info(all_paths)$size>1)]

# Gathering data from all files
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t")

all_filenames <- all_paths %>%
  basename() %>%
  as.list()

all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# Combining all files based on paths
all_result <- rbindlist(all_lists, fill = T) %>% 
                rename(., "Sample"="V1") 
                
# Creating a dataframe and cleaning up Sample names
merged_df <- all_result %>% mutate(Sample=gsub("_growth_prediction.txt", "", Sample)) 

# Writing dataframe to file
write.table(merged_df, file = snakemake@output[["DF"]], sep="\t", row.names=FALSE, quote=FALSE)

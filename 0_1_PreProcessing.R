# Clean the taxonomy string and put it to the genus level
c

# Merge Tab and Taxonomy, remove spaces in taxonomy and Unassigned/Archael/Chloroplast and mitochonria taxonomies, adapt 
MergeTabTax <- function(table,taxonomy){
out_table <- as.data.frame(table)
out_table$Taxonomy = taxonomy$Taxonomy

# remove spaces
out_table$Taxonomy <- gsub("\\s", "", out_table$Taxonomy)

# filter based on taxonomy
out_table<-out_table[!grepl("Unassigned", out_table$Taxonomy),]
out_table<-out_table[!grepl("k__Archaea", out_table$Taxonomy),]
out_table<-out_table[!grepl("Chloroplast", out_table$Taxonomy),]
out_table<-out_table[!grepl("k__mitochondria", out_table$Taxonomy),]

return(out_table)
}

# Filter samples by sequence number, remove SVs with no more sequences
FilterSeqNumber <- function(table, seq_number, sample_names){
table_filtered = table[,colSums(table[,sample_names]) >= seq_number]
table_filtered$Taxonomy = table$Taxonomy
table_filtered$SV_ID = table$OTU.ID
sample_names = sample_names[sample_names %in% colnames(table_filtered)]
table_filtered = table_filtered[rowSums(table_filtered[,sample_names]) > 0,]
return(table_filtered)}

# Update the metadata table based on the samples that were kept
UpdateMetadata <- function(metadata,table){
  filtered_metadata = metadata[metadata$sample_name %in% colnames(table),]
  return(filtered_metadata)
}

# Create function to sub-sample EUCI data
SubSampleTable <- function(table, samples_to_keep){
  table_filtered = table[rowSums(table[,samples_to_keep]) > 0,c('Taxonomy','SV_ID',samples_to_keep)]
  return(table_filtered)
}

# Create function to create presence/absence matrix for the N most abundant SVs
MostAbundantPresAbs <- function(table, number_to_keep, samples){
  most_abundant = table
  most_abundant[samples] = sweep(most_abundant[samples],2,colSums(most_abundant[samples]),'/')
  most_abundant$row_sum <- rowSums(most_abundant[,samples])
  most_abundant = most_abundant[order(most_abundant$row_sum, decreasing = T),]
  most_abundant = most_abundant[1:number_to_keep,]
  most_abundant[most_abundant > 0] = 1
  most_abundant = most_abundant[,c(samples, 'Taxonomy', 'SV_ID')]
  return(most_abundant)}

#########################################################################################################################
# 0. Data pre-processing
setwd('/Users/admin/Documents/Academia/PhD/Chapter I/')

# Load PP1
PP1_raw_tab = read.table('Datasets/Data/PP1_raw_tab.tsv', sep = '\t', header = T)
PP1_raw_tax = read.table('Datasets/Data/PP1_raw_tax.tsv', sep = '\t', header = T)
PP1_raw_tax$Taxonomy = as.character(PP1_raw_tax$Taxonomy)
PP1_metadata_raw = read.table('Datasets/Metadata/PP1_metadata_raw.csv', sep = ',', header = T)
PP1_EUCI_samples = as.vector(PP1_metadata_raw$sample_name[PP1_metadata_raw$EUCI == 'Yes'][PP1_metadata_raw$sample_name[PP1_metadata_raw$EUCI == 'Yes'] %in% colnames(PP1_raw_tab)])
PP1_COMP_samples = as.vector(PP1_metadata_raw$sample_name[PP1_metadata_raw$EUCI == 'No'][PP1_metadata_raw$sample_name[PP1_metadata_raw$EUCI == 'No'] %in% colnames(PP1_raw_tab)])

# Load PP2
PP2_raw_tab = read.table('Datasets/Data/PP2_raw_tab.tsv', sep = '\t', header = T)
PP2_raw_tax = read.table('Datasets/Data/PP2_raw_tax.tsv', sep = '\t', header = T)
PP2_raw_tax$Taxonomy = as.character(PP2_raw_tax$Taxonomy)
PP2_metadata_raw = read.table('Datasets/Metadata/PP2_metadata_raw.csv', sep = ',', header = T)
PP2_EUCI_samples = as.vector(PP2_metadata_raw$sample_name[PP2_metadata_rawt$EUCI == 'Yes'][PP2_metadata_raw$sample_name[PP2_metadata_raw$EUCI == 'Yes'] %in% colnames(PP2_raw_tab)])
PP2_COMP_samples = as.vector(PP2_metadata_raw$sample_name[PP2_metadata_raw$EUCI == 'No'][PP2_metadata_raw$sample_name[PP2_metadata_raw$EUCI == 'No'] %in% colnames(PP2_raw_tab)])

# Clean the taxonomy string
PP1_raw_tax$Taxonomy = CleanTaxString(PP1_raw_tax$Taxonomy)
PP2_raw_tax$Taxonomy = CleanTaxString(PP2_raw_tax$Taxonomy)

# Merge taxonomy and tables
PP1_tab = MergeTabTax(PP1_raw_tab, PP1_raw_tax)
PP2_tab = MergeTabTax(PP2_raw_tab, PP2_raw_tax)

# Remove samples with less than 2000 sequences
PP1_tab = FilterSeqNumber(PP1_tab, 2000, c(PP1_EUCI_samples, PP1_COMP_samples))
PP2_tab = FilterSeqNumber(PP2_tab, 2000, c(PP2_EUCI_samples, PP2_COMP_samples))
write.table(PP1_tab, file = 'Datasets/Data/PP1_tab.tsv', sep = "\t", row.names = FALSE)
write.table(PP2_tab, file = 'Datasets/Data/PP2_tab.tsv', sep = "\t", row.names = FALSE)

PP1_metadata = UpdateMetadata(PP1_metadata_raw, PP1_tab)
PP2_metadata = UpdateMetadata(PP2_metadata_raw, PP2_tab)
write.table(PP1_metadata, file = 'Datasets/Metadata/PP1_metadata.tsv', sep = "\t", row.names = FALSE)
write.table(PP2_metadata, file = 'Datasets/Metadata/PP2_metadata.tsv', sep = "\t", row.names = FALSE)

# Create EUCI sub sample
PP1_EUCI_samples = as.vector(PP1_metadata$sample_name[PP1_metadata$EUCI == 'Yes'][PP1_metadata$sample_name[PP1_metadata$EUCI == 'Yes'] %in% colnames(PP1_tab)])
PP2_EUCI_samples = as.vector(PP2_metadata$sample_name[PP2_metadata$EUCI == 'Yes'][PP2_metadata$sample_name[PP2_metadata$EUCI == 'Yes'] %in% colnames(PP2_tab)])
PP1_EUCI_tab = SubSampleTable(PP1_tab, PP1_EUCI_samples)
PP2_EUCI_tab = SubSampleTable(PP2_tab, PP2_EUCI_samples)
write.table(PP1_EUCI_tab, file = 'Datasets/Data/PP1_EUCI_tab.tsv', sep = "\t", row.names = FALSE)
write.table(PP2_EUCI_tab, file = 'Datasets/Data/PP2_EUCI_tab.tsv', sep = "\t", row.names = FALSE)

# Create Random-forest data (10'000 most abundant SVs, presence/absence matrix)
PP1_COMP_samples = as.vector(PP1_metadata$sample_name[PP1_metadata$EUCI == 'No'][PP1_metadata$sample_name[PP1_metadata$EUCI == 'No'] %in% colnames(PP1_tab)])
PP2_COMP_samples = as.vector(PP2_metadata$sample_name[PP2_metadata$EUCI == 'No'][PP2_metadata$sample_name[PP2_metadata$EUCI == 'No'] %in% colnames(PP2_tab)])

PP1_RF_tab = MostAbundantPresAbs(PP1_tab, 10000, c(PP1_EUCI_samples, PP1_COMP_samples))
PP2_RF_tab = MostAbundantPresAbs(PP2_tab, 10000, c(PP2_EUCI_samples, PP2_COMP_samples))
write.table(PP1_RF_tab, file = 'Datasets/Data/PP1_RF_tab.tsv', sep = "\t", row.names = FALSE)
write.table(PP2_RF_tab, file = 'Datasets/Data/PP2_RF_tab.tsv', sep = "\t", row.names = FALSE)



library(data.table)
library(tibble)
library(stringr)
library(Biostrings)

MergeTabTax <- function(ASV_table, Tax_table){
  merged_table <- ASV_table
  merged_table$Taxonomy = as.character(Tax_table$Taxon)
  return(merged_table)
}

CleanTaxString <- function(Tax_column){
  print('Removing species level taxonomies...')
  Clean_Tax_column = vapply(Tax_column, function(x) {strsplit(x, '; s__')[[1]][1]}, FUN.VALUE = character(1))
  return(Clean_Tax_column)}

FilterTabTax <- function(merged_table){
  print('Removing sequences with taxonomy classification containing: Unassigned, Eukaryota, Chloroplast, Mitochondria and Archaea...')
  merged_table<-merged_table[!grepl("Unassigned", merged_table$Taxonomy),]
  merged_table<-merged_table[!grepl("Archaea", merged_table$Taxonomy),]
  merged_table<-merged_table[!grepl("Eukaryota", merged_table$Taxonomy),]
  merged_table<-merged_table[!grepl("Chloroplast", merged_table$Taxonomy),]
  merged_table<-merged_table[!grepl("Mitochondria", merged_table$Taxonomy),]
  return(merged_table)}

FilterSeqNumber <- function(merged_table, seq_number){
  # remove columns with less than seq_number sequences after taxonomy filtering
  print('Initial filtering: removing samples with n or less sequences...')
  sample_names = colnames(merged_table)[!(colnames(merged_table) %in% c('ASV','Taxonomy'))]
  samples_sums = colSums(merged_table[,..sample_names])
  columns_to_keep = c('ASV', 'Taxonomy', sample_names[samples_sums > seq_number])
  merged_table = merged_table[,..columns_to_keep]
  
  # remove ASVs that now have 0 counts due to filtering or are singletons
  print('Removing ASVs with less than 2 counts...')
  sample_names = colnames(merged_table)[!(colnames(merged_table) %in% c('ASV','Taxonomy'))]
  merged_table = merged_table[rowSums(merged_table[,..sample_names])>1,]

  # filter again the samples, in case singleton filtering removed samples
  print('Filtering again the samples with n or less sequences after singletons removal...')
  updated_sample_names = colnames(merged_table)[!(colnames(merged_table) %in% c('ASV','Taxonomy'))]
  updated_samples_sums = colSums(merged_table[,..updated_sample_names])
  updated_columns_to_keep = c('ASV', 'Taxonomy', updated_sample_names[updated_samples_sums > seq_number])
  filtered_table = merged_table[,..updated_columns_to_keep]
  
  # remove again ASVs that now have 0 counts due to filtering
  print('Removing ASVs with 0 counts...')
  filtered_table[, `:=`(ASV_sum = rowSums(.SD)), .SDcols=updated_sample_names[updated_samples_sums > seq_number]]
  filtered_table = filtered_table[ASV_sum > 0,]
  filtered_table = filtered_table[,ASV_sum:=NULL]

  return(filtered_table)}

UpdateMetadata <- function(metadata,table){
  print('Updating metadata...')
  filtered_metadata = metadata[metadata$Sample %in% colnames(table),]
  filtered_metadata$Longitude = UpdateCoord(strsplit(as.character(filtered_metadata$Longitude),''))
  filtered_metadata$Latitude = UpdateCoord(strsplit(as.character(filtered_metadata$Latitude),''))
  return(filtered_metadata)
}

UpdateMetadataMTG <- function(metadata,table,kegg){
  print('Updating metadata...')
  keep_samples = union(colnames(table), colnames(kegg))
  print(keep_samples)
  filtered_metadata = metadata[metadata$Sample %in% keep_samples,]
  filtered_metadata$Longitude = UpdateCoord(strsplit(as.character(filtered_metadata$Longitude),''))
  filtered_metadata$Latitude = UpdateCoord(strsplit(as.character(filtered_metadata$Latitude),''))
  return(filtered_metadata)
}

UpdateCoord <- function(gps_list){
  new_gps_list = c()
  for(gps_coord in gps_list){
    if (is.na(gps_coord)){
      new_gps_list = c(new_gps_list, paste(gps_coord, collapse = ''))}
    else{
      if (gps_coord[length(gps_coord)] == 'N'){
        new_coord = paste(c(gps_coord[1:length(gps_coord)-1]), collapse = '')
        new_gps_list = c(new_gps_list,new_coord)}
      else if (gps_coord[length(gps_coord)] == 'S'){
        new_coord = paste(c('-',gps_coord[1:length(gps_coord)-1]), collapse = '')
        new_gps_list = c(new_gps_list,new_coord)}
      else if (gps_coord[length(gps_coord)] == 'E'){
        new_coord = paste(c(gps_coord[1:length(gps_coord)-1]), collapse = '')
        new_gps_list = c(new_gps_list,new_coord)}
      else if (gps_coord[length(gps_coord)] == 'W'){
        new_coord = paste(c('-',gps_coord[1:length(gps_coord)-1]), collapse = '')
        new_gps_list = c(new_gps_list,new_coord)}
      else{new_gps_list = c(new_gps_list, paste(gps_coord, collapse = ''))}}
  }
  return(as.vector(new_gps_list))}

SamplesBadTax <- function(data, metadata){
            samples_to_remove = c()
            samples = unique(metadata$Sample)[unique(metadata$Sample) %in% colnames(data)]
            for (sample in samples){
              cols_to_keep = c('ASV', 'Taxonomy', sample)
              sample_data = as.data.frame(data[, ..cols_to_keep])
              tax_vec = as.vector(sample_data$Taxonomy[sample_data[sample] > 0])
              if (length(tax_vec) > 0){
              prop_phylum = length(tax_vec[tax_vec %like% 'p__'])/length(tax_vec)
              prop_genus = length(tax_vec[tax_vec %like% 'g__'])/length(tax_vec)
              if ((prop_phylum < 0.75) & (prop_genus < 0.5)){
                print('-----------------------------')
                print(paste0('Sample: ', sample))
                print(paste0('  - ', length(tax_vec), ' ASVs'))
                print(paste0('  - Proportion of phylum taxonomy:', prop_phylum))
                print(paste0('  - Proportion of genus taxonomy:', prop_genus))
                samples_to_remove = c(samples_to_remove, sample)}}
              else{samples_to_remove = c(samples_to_remove, sample)}}
            return(samples_to_remove)}

setwd('SET_THE_WORKING_DIRECTORY_HERE')

###############################################################################################
# PP1
PP1_raw_tab = fread('Data/raw/PP1_all_merged_raw_table.tsv')
PP1_raw_tax = fread('Data/raw/PP1_all_merged_raw_tax.tsv')
PP1_raw_metadata = as_tibble(fread('Metadata/PP1_metadata_raw.csv', header = TRUE))
# Add taxonomy to the table, keep only taxonomy up to genus, filter based on taxonomy
PP1_merged_table = MergeTabTax(PP1_raw_tab, PP1_raw_tax)
PP1_merged_table$Taxonomy = as.vector(CleanTaxString(PP1_merged_table$Taxonomy))
PP1_merged_table = FilterTabTax(PP1_merged_table)
# Remove samples with bad quality taxonomy
PP1_samples_to_remove = SamplesBadTax(PP1_merged_table, PP1_raw_metadata)
intersect(PP1_raw_metadata$Sample[PP1_raw_metadata$Cryosphere == 'Yes'], PP1_samples_to_remove)
PP1_merged_table = PP1_merged_table[,!PP1_samples_to_remove, with = FALSE]
# Filter samples with less thank 5k sequences, remove singletons, filter again
PP1_merged_table = FilterSeqNumber(PP1_merged_table, 4999)
# Write Table
write.table(PP1_merged_table, file='Data/PP1_table.tsv', sep = "\t", row.names = FALSE)
# Update metadata
PP1_metadata = UpdateMetadata(PP1_raw_metadata,PP1_merged_table)
PP1_metadata$Latitude = UpdateCoord(PP1_metadata$Latitude)
PP1_metadata$Longitude = UpdateCoord(PP1_metadata$Longitude)
write.table(PP1_metadata, file='Metadata/PP1_metadata.tsv', sep = "\t", row.names = FALSE)

###############################################################################################
# PP2
PP2_raw_tab = fread('Data/raw/PP2_all_merged_raw_table.tsv')
PP2_raw_tax = fread('Data/raw/PP2_all_merged_raw_tax.tsv')
PP2_raw_metadata = as_tibble(fread('Metadata/PP2_metadata_raw.csv', header = TRUE))
# Add taxonomy to the table, keep only taxonomy up to genus, filter based on taxonomy
PP2_merged_table = MergeTabTax(PP2_raw_tab, PP2_raw_tax)
PP2_merged_table$Taxonomy = as.vector(CleanTaxString(PP2_merged_table$Taxonomy))
PP2_merged_table = FilterTabTax(PP2_merged_table)
# Remove samples with bad quality taxonomy
PP2_samples_to_remove = SamplesBadTax(PP2_merged_table, PP2_raw_metadata)
intersect(PP2_raw_metadata$Sample[PP2_raw_metadata$Cryosphere == 'Yes'], PP2_samples_to_remove)
PP2_merged_table = PP2_merged_table[,!PP2_samples_to_remove, with = FALSE]
# Filter samples with less thank 5k sequences
PP2_merged_table = FilterSeqNumber(PP2_merged_table, 4999)
# Write Table
write.table(PP2_merged_table, file='Data/PP2_table.tsv', sep = "\t", row.names = FALSE)
# Update metadata
PP2_metadata = UpdateMetadata(PP2_raw_metadata, PP2_merged_table)
PP2_metadata$Latitude = UpdateCoord(PP2_metadata$Latitude)
PP2_metadata$Longitude = UpdateCoord(PP2_metadata$Longitude)
write.table(PP2_metadata, file='Metadata/PP2_metadata.tsv', sep = "\t", row.names = FALSE)

# Sequences
library(seqinr)
PP1_raw_fasta = read.fasta(file = "Data/raw/PP1_all_merged_raw_seqs.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
PP2_raw_fasta = read.fasta(file = "Data/raw/PP2_all_merged_raw_seqs.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

PP1_fasta = PP1_raw_fasta[names(PP1_raw_fasta) %in% PP1_merged_table$ASV]
PP2_fasta = PP2_raw_fasta[names(PP2_raw_fasta) %in% PP2_merged_table$ASV]

write.fasta(sequences = PP1_fasta, names = names(PP1_fasta), file.out = 'Data/trees/PP1_ASV_seqs.fasta')
write.fasta(sequences = PP2_fasta, names = names(PP2_fasta), file.out = 'Data/trees/PP2_ASV_seqs.fasta')

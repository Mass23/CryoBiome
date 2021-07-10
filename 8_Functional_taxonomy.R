library(metacoder)
library(taxize)
library(rentrez)
library(foreach)
library(dplyr) 
library(ggplot2)
library(ggrepel)
# Colors light:
#---------------
# PP1: #1F9DDB
# PP2: #DB0600
# MTG tax: #6BC90C
# MTG func: #29C275  

# Colors dark:
#--------------
# PP1: #1A87BD
# PP2: #BA0600
# MTG tax: #5CAD0A
# MTG func: #23A671

setwd('/Users/admin/Documents/Academia/PhD/Chapter I/')
############################################################################################################
KEGG_taxonomy_table = read.csv('7_MTG_functional/KEGG_sign_tax.tsv', sep = '\t')
GetGenus <- function(id){
  genera_vec = c()
  tax = extract_tax_data(id,key="taxon_id",regex = "(.+)")
  if ('genus' %in% tax$data$tax_data$ncbi_rank){
    genus = tax$data$tax_data$ncbi_name[tax$data$tax_data$ncbi_rank == 'genus']
    return(genus)}
  else{return('Unassigned')}}

set_entrez_key("1b37cec2fdc64d05c5cc5337b17b09b4c609")
all_genera <- foreach(i=1:nrow(KEGG_taxonomy_table), .export = c('GetGenus'), .combine=rbind, .verbose = TRUE) %do% {
  taxids = KEGG_taxonomy_table$Taxid[i]
  taxids = as.vector(strsplit(taxids, ','))[[1]]
  df_tax = data.frame('Genus'=c(), 'KEGG'=c(), 'Sample'=c(), 'n_match'=c())
  for (id in taxids){
    genus = GetGenus(id)
    df_tax = rbind(df_tax, data.frame('Genus'=genus, 'KEGG'=KEGG_taxonomy_table$KEGG[i], 'Sample'=KEGG_taxonomy_table$Sample[i], 'n_match'=1))}
  return(df_tax)
}

all_genera = all_genera %>% group_by(Genus) %>% summarise(sum_match = sum(n_match), n_samp = n_distinct(Sample), n_KO = n_distinct(KEGG))


KEGG_taxonomy = read.table('8_Functional_taxonomy/KEGG_sign_tax_genera.csv')
KEGG_taxonomy = KEGG_taxonomy[KEGG_taxonomy$Genus != 'Unassigned',]

strict_file = read.csv('4_Taxonomic_analysis/strict_genera.csv')
strict_genera = strict_file$NCBI
relaxed_file = read.csv('4_Taxonomic_analysis/relaxed_genera.csv')
relaxed_genera = relaxed_file$NCBI

KEGG_taxonomy$Group = 'Others'
KEGG_taxonomy$Group[KEGG_taxonomy$Genus %in% relaxed_genera] = 'Ancillary'
KEGG_taxonomy$Group[KEGG_taxonomy$Genus %in% strict_genera] = 'Core'
KEGG_taxonomy$Group  <- factor(KEGG_taxonomy$Group , levels = c("Core","Ancillary", "Others"))

ggplot(KEGG_taxonomy, aes(x=n_KO, y=n_samp/35, size=sum_match, colour=Group)) +  theme_bw() + geom_point() +
  geom_label_repel(show.legend = F,aes(label = ifelse(sum_match>2000,as.character(Genus),'')), box.padding = 0.7) + 
  scale_colour_manual(values=c('#1D8C78', '#29C275', 'grey')) + labs(y='Prevalence in EUCI', x='EUCI enriched KOs', size='Contigs n.')
ggsave('8_Functional_taxonomy/KEGG_taxonomy.pdf', width = 8, height = 6)

ggplot(KEGG_taxonomy) + geom_boxplot(aes(y=sum_match, x = Group, colour = Group)) + scale_y_log10() + theme_bw() +
  scale_colour_manual(values=c('#1D8C78', '#29C275', 'grey')) + ylab('Contigs n.') + xlab('')
ggsave('8_Functional_taxonomy/KEGG_n_matches.pdf', width = 4, height = 3)






library(taxize)
library(rentrez)
library(foreach)
library(dplyr) 
library(ggplot2)
library(taxa)
library(parallel)
library(ggpubr)
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

setwd('/Users/mabourqu/Desktop/cryobiome_revisions/')
############################################################################################################
# Gather Genus info, first run the KEGG_
KEGG_taxonomy_table = read.csv('Data/KEGG_sign_tax.tsv', sep = '\t')
GetGenus <- function(id){
  genera_vec = c()
  tax = extract_tax_data(id,key="taxon_id",regex = "(.+)")
  if ('genus' %in% tax$data$tax_data$ncbi_rank){
    genus = tax$data$tax_data$ncbi_name[tax$data$tax_data$ncbi_rank == 'genus']
    return(genus)}
  else{return('Unassigned')}}

set_entrez_key("1b37cec2fdc64d05c5cc5337b17b09b4c609")

cl <- parallel::makeCluster(10)

all_genera <- foreach(i=1:nrow(KEGG_taxonomy_table), .export = c('GetGenus'), .combine=rbind, .verbose = TRUE) %dopar% {
  taxids = as.character(KEGG_taxonomy_table$Taxid[i])
  if (length(taxids) > 0){
  taxids = as.vector(strsplit(taxids, ','))[[1]]
  df_tax = data.frame('Genus'=c(), 'KEGG'=c(), 'Sample'=c(), 'n_match'=c())
  for (id in taxids){
    genus = GetGenus(id)
    df_tax = rbind(df_tax, data.frame('Genus'=genus, 'KEGG'=KEGG_taxonomy_table$KEGG[i], 'Sample'=KEGG_taxonomy_table$Sample[i], 'n_match'=1))}
  return(df_tax)
}}

all_genera = all_genera %>% group_by(Genus) %>% summarise(sum_match = sum(n_match), n_samp = n_distinct(Sample), n_KO = n_distinct(KEGG))
write.csv(all_genera, file='Data/KEGG_sign_tax_genera.csv',sep=',')

############################################################################################################
KEGG_taxonomy = read.table('Data/KEGG_sign_tax_genera.csv')
KEGG_taxonomy = KEGG_taxonomy[KEGG_taxonomy$Genus != 'Unassigned',]

# Assign the genera groups by partial string match
over_res = read.csv('Data/Amplicon_ancom_res.csv')
over_res$DA = as.character(over_res$DA)
over_res$DA[over_res$DA == 'Overrepresented'] = 'Cryo.'
over_res$SilvaGenus = vapply(over_res$taxa_id, function(x){gsub('_',' ',strsplit(as.character(x), split='g__')[[1]][2])}, FUN.VALUE = character(1))
over_genera = over_res$SilvaGenus[over_res$DA == 'Cryo.']

KEGG_taxonomy$DA = vapply(KEGG_taxonomy$Genus, function(x){ifelse(sum(grep(x, over_genera)) > 0, 'Cryo.','Others')}, FUN.VALUE = character(1))

# Plots
ggplot(KEGG_taxonomy, aes(x=n_KO, y=n_samp/35, size=sum_match, colour=DA)) +  theme_linedraw() + geom_point() +
  geom_label_repel(show.legend = F,aes(label = ifelse(sum_match>2000,as.character(Genus),'')), box.padding = 0.6)  + theme(panel.grid = element_line(colour = 'darkgrey')) +
  scale_colour_manual(values=c('#00A087FF','darkgrey')) + labs(y='Prevalence in cryospheric samples', x='Cryosphere enriched KOs', size='Contigs n.')
ggsave('3_Functional_analysis/3_2_Taxonomy/KEGG_taxonomy.pdf', width = 6, height = 5)

ggplot(KEGG_taxonomy) + geom_density(aes(x=sum_match, fill = DA), alpha=0.5)  + scale_x_log10() + theme_linedraw() + guides(fill=FALSE) + 
  scale_fill_manual(values=c('#00A087FF','darkgrey')) + xlab('Contigs n.') + ylab('Density') +
  theme(panel.grid = element_line(colour = 'darkgrey')) 
ggsave('3_Functional_analysis/3_2_Taxonomy/KEGG_tax_n_matches.pdf', width = 3, height = 3)
wilcox.test(KEGG_taxonomy$sum_match[KEGG_taxonomy$DA =='Cryo.'],KEGG_taxonomy$sum_match[KEGG_taxonomy$DA == 'Others'])
# data:  KEGG_taxonomy$sum_match[KEGG_taxonomy$DA == "Cryo."] and KEGG_taxonomy$sum_match[KEGG_taxonomy$DA == "Others"]
# W = 413103, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0



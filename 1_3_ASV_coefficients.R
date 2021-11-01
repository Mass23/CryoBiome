library(ggtree)
library(treeio)
library(ape)
library(seqinr)
library(ggplot2)
library(phangorn)
library(ggstance)
library(RColorBrewer)
library(ggtreeExtra)
library(ggnewscale)
library(phytools)
library(reshape2)

# colors light:
# PP1: #1F9DDB
# PP2: #DB0600
# MTG tax: #6BC90C
# MTG func: #29C275

# colors dark:
# PP1: #1A87BD
# PP2: #BA0600
# MTG tax: #5CAD0A
# MTG func: #23A671
setwd('/Users/mabourqu/Documents/PhD/C1/')

####################################################################################################
# 1. Phylogenetics
# Coefficients of SVM across the tree
PP1_tab_abund = read.csv('Data/trees/PP1_005_ASVs_table.csv')
PP1_ASV_coefs = read.csv('Data/PP1_Logistics_coefs.csv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep ='\t')
PP1_metadata$Ecosystem[PP1_metadata$Ecosystem %in% c('Ice','Snow')] = 'Ice/Snow'
rownames(PP1_metadata) = PP1_metadata$Sample
PP1_tree = read.newick(file='Data/trees/PP1_005_trim.fasta.treefile')
PP1_tree = midpoint(PP1_tree)

PP2_tab_abund = read.csv('Data/trees/PP2_005_ASVs_table.csv')
PP2_ASV_coefs = read.csv('Data/PP2_Logistics_coefs.csv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep ='\t')
PP2_metadata$Ecosystem[PP2_metadata$Ecosystem %in% c('Ice','Snow')] = 'Ice/Snow'
PP2_tree = read.newick(file='Data/trees/PP2_005_trim.fasta.treefile')
PP2_tree = midpoint(PP2_tree)

# Coef matrix
PP1_coef_tab = data.frame(ASV = PP1_ASV_coefs$ASV, Coefficient = PP1_ASV_coefs$Importance)
PP1_coef_tab$Phylum = vapply(1:nrow(PP1_ASV_coefs), function(x) {ifelse(grepl('; p__', PP1_ASV_coefs$Taxonomy[x], fixed=TRUE), strsplit(strsplit(PP1_ASV_coefs$Taxonomy[x], split = '; c__')[[1]][1], split = '; p__')[[1]][2], 'Unassigned')}, FUN.VALUE = character(1))
PP1_coef_tab$Phylum[!(PP1_coef_tab$Phylum %in% c('Acidobacteriota','Actinobacteriota','Armatimonadata','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota'))] = 'Others'
# Taking the log as we have the exp of the coefs
PP1_coef_tab$Coefficient = log(PP1_coef_tab$Coefficient)
PP1_coef_tab_pos_only = PP1_coef_tab
PP1_coef_tab_pos_only[PP1_coef_tab_pos_only < 0] = 0

PP2_coef_tab = data.frame(ASV = PP2_ASV_coefs$ASV, Coefficient = PP2_ASV_coefs$Importance)
PP2_coef_tab$Phylum = vapply(1:nrow(PP2_ASV_coefs), function(x) {ifelse(grepl('; p__', PP2_ASV_coefs$Taxonomy[x], fixed=TRUE), strsplit(strsplit(PP2_ASV_coefs$Taxonomy[x], split = '; c__')[[1]][1], split = '; p__')[[1]][2], 'Unassigned')}, FUN.VALUE = character(1))
PP2_coef_tab$Phylum[!(PP2_coef_tab$Phylum %in% c('Acidobacteriota','Actinobacteriota','Armatimonadata','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota'))] = 'Others'
# Taking the log as we have the exp of the coefs
PP2_coef_tab$Coefficient = log(PP2_coef_tab$Coefficient)
PP2_coef_tab_pos_only = PP2_coef_tab
PP2_coef_tab_pos_only[PP2_coef_tab_pos_only < 0] = 0

# Ecosystem matrix
PP1_eco_tab = expand.grid(ASV = PP1_tab_abund$ASV, Ecosystem = c('Ice/Snow', 'Terrestrial', 'Marine', 'Freshwater'))
PP1_eco_tab$Phylum = vapply(PP1_eco_tab$ASV, function(x){PP1_coef_tab$Phylum[PP1_coef_tab$ASV == x]}, character(1))
PP1_eco_tab$Presence = vapply(1:nrow(PP1_eco_tab), function(x){ifelse(sum(PP1_tab_abund[PP1_eco_tab$ASV[x], PP1_metadata$Sample[(PP1_metadata$Ecosystem == PP1_eco_tab$Ecosystem[x]) & (PP1_metadata$Cryosphere == 'Yes')]]) > 0,'Yes','No')},character(1))
PP1_eco_tab$Phylum[!(PP1_eco_tab$Phylum %in% c('Acidobacteriota','Actinobacteriota','Armatimonadata','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota'))] = 'Others'

PP2_eco_tab = expand.grid(ASV = PP2_tab_abund$ASV, Ecosystem = c('Ice/Snow', 'Terrestrial', 'Marine', 'Freshwater'))
PP2_eco_tab$Phylum = vapply(PP2_eco_tab$ASV, function(x){PP2_coef_tab$Phylum[PP2_coef_tab$ASV == x]}, character(1))
PP2_eco_tab$Presence = vapply(1:nrow(PP2_eco_tab), function(x){ifelse(sum(PP2_tab_abund[PP2_eco_tab$ASV[x], PP2_metadata$Sample[(PP2_metadata$Ecosystem == PP2_eco_tab$Ecosystem[x]) & (PP2_metadata$Cryosphere == 'Yes')]]) > 0,'Yes','No')},character(1))
PP2_eco_tab$Phylum[!(PP2_eco_tab$Phylum %in% c('Acidobacteriota','Actinobacteriota','Armatimonadata','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota'))] = 'Others'

write.csv(PP1_eco_tab,file='Data/trees/PP1_eco_tab_tree.csv')
write.csv(PP2_eco_tab,file='Data/trees/PP2_eco_tab_tree.csv')
# Trees
# phylum: c('Acidobacteriota','Actinobacteriota','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Others,'Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota')
# colors: c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')
# The circular layout tree.
PP1_p <- ggtree(PP1_tree, layout="fan", size=0.15, open.angle=10)

PP1_p1 <- PP1_p + geom_fruit(data=PP1_eco_tab, geom=geom_tile,
                             mapping=aes(y=ASV, x=Ecosystem, fill=Phylum, alpha=Presence),
                             offset = 0.04, size = 0.02, axis.params = list(axis = "x", text.angle = -90, text.size = 2.5)) + scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')) + 
                             scale_alpha_discrete() + new_scale_fill() +
                  geom_fruit(data=PP1_coef_tab_pos_only, geom=geom_bar,
                             mapping=aes(y=ASV, x=Coefficient, fill=Phylum, colour=Phylum),
                             pwidth=0.38, orientation="y", stat="identity", axis.params=list(
                             axis = "x",text.size = 1.8, hjust = 1,vjust = 0.5, nbreak = 3), grid.params=list()) + 
                             scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')) +
                             scale_color_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF'))

ggsave('1_ASV_analysis/1_3_ASV_coefficients/PP1_tree.png', plot = PP1_p1, width = 12, height = 10)



PP2_p <- ggtree(PP2_tree, layout="fan", size=0.15, open.angle=10)

PP2_p1 <- PP2_p + geom_fruit(data=PP2_eco_tab, geom=geom_tile,
                             mapping=aes(y=ASV, x=Ecosystem, fill=Phylum, alpha=Presence),
                             offset = 0.04, size = 0.02, axis.params = list(axis = "x", text.angle = -90, text.size = 2.5)) + scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')) + 
                             scale_alpha_discrete() + new_scale_fill() +
                  geom_fruit(data=PP2_coef_tab_pos_only, geom=geom_bar,
                             mapping=aes(y=ASV, x=Coefficient, fill=Phylum, colour=Phylum),
                             pwidth=0.38, orientation="y", stat="identity", axis.params=list(
                             axis = "x",text.size = 1.8, hjust = 1,vjust = 0.5, nbreak = 3), grid.params=list()) + 
  scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')) +
  scale_color_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF'))

ggsave('1_ASV_analysis/1_3_ASV_coefficients/PP2_tree.png', plot = PP2_p1, width = 12, height = 10)


###########################################################################################################################
# 2. Taxonomy
PP1_ASV_coefs = read.csv('Data/PP1_Logistics_coefs.csv')
PP1_ASV_coefs$Phylum = vapply(PP1_ASV_coefs$Taxonomy, function(x)strsplit(strsplit(x,split='; c__')[[1]][1], split='; p__')[[1]][2], FUN.VALUE = character(1)) 
PP1_ASV_coefs$Genus = vapply(PP1_ASV_coefs$Taxonomy, function(x)strsplit(x,split='; g__')[[1]][2], FUN.VALUE = character(1)) 
PP1_ASV_coefs$Dataset = 'PP1'

PP2_ASV_coefs = read.csv('Data/PP2_Logistics_coefs.csv')
PP2_ASV_coefs$Phylum = vapply(PP2_ASV_coefs$Taxonomy, function(x)strsplit(strsplit(x,split='; c__')[[1]][1], split='; p__')[[1]][2], FUN.VALUE = character(1)) 
PP2_ASV_coefs$Genus = vapply(PP2_ASV_coefs$Taxonomy, function(x)strsplit(x,split='; g__')[[1]][2], FUN.VALUE = character(1)) 
PP2_ASV_coefs$Dataset = 'PP2'

ASV_coefs = rbind(PP1_ASV_coefs, PP2_ASV_coefs)

library(dplyr)
highest_coefs = ASV_coefs[ASV_coefs$Importance > 1,]
phyla_to_keep = names(table(highest_coefs$Phylum))[rowSums(table(highest_coefs$Phylum, highest_coefs$Dataset) > 20) == 2]
genera_to_keep = names(table(highest_coefs$Genus))[rowSums(table(highest_coefs$Genus, highest_coefs$Dataset) > 10) == 2]
genera_to_keep = genera_to_keep[genera_to_keep != 'uncultured']
highest_coefs$Phylum[!(highest_coefs$Phylum %in% c('Acidobacteriota','Actinobacteriota','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota'))] = 'Others'
highest_coefs[highest_coefs$Phylum %in% phyla_to_keep,] %>%  ggplot(aes(y=Phylum, fill=Phylum)) + xlab('N [ASVs with OR > 1]') + ylab('') + xlim(c(0,1300)) +
  geom_bar(stat='count') + facet_grid(Phylum~Dataset, scales = 'free_y', space = 'free_y') + theme_linedraw() + theme(strip.text.y = element_blank(), legend.position = 'none') + 
  scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')) 
ggsave('1_ASV_analysis/1_3_ASV_coefficients/1_3_Phylum_coefs.pdf', width=4.5, height = 4.5)

highest_coefs[highest_coefs$Genus %in% genera_to_keep,] %>% ggplot(aes(x=..count..,y=Genus, fill=Phylum)) + xlab('N [ASVs with OR > 1]') + ylab('') +
  geom_bar(stat='count') + facet_grid(Phylum~Dataset, scales = 'free_y', space = 'free_y') + theme_linedraw() + theme(strip.text.y = element_blank(), legend.position = 'none') + 
  scale_fill_manual(values = c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#7E6148FF','#3C5488FF','#8491B4FF')) 
ggsave('1_ASV_analysis/1_3_ASV_coefficients/1_3_Genera_coefs.pdf', width=4.5, height = 5)


library(data.table)
library(dplyr) 
library(ggplot2)
library(ape)
library(seqinr)
library(foreach)
library(doMC)
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

aln_files = list.files(path = "Data/Filtered/EUCI_genus_aln/.")
PP1_aln_files = aln_files[grep('PP1_', aln_files)]
PP2_aln_files = aln_files[grep('PP2_', aln_files)]

tree_files = list.files(path = 'Data/Filtered/EUCI_genus_trees/.')
PP1_tree_files = tree_files[grep('PP1_', tree_files)]
PP2_tree_files = tree_files[grep('PP2_', tree_files)]

PP1_genera = as.vector(vapply(PP1_tree_files, function(x) strsplit(strsplit(x,'_trim.fasta.treefile')[[1]][1], 'PP1_')[[1]][2], FUN.VALUE = character(1)))
PP2_genera = as.vector(vapply(PP2_tree_files, function(x) strsplit(strsplit(x,'_trim.fasta.treefile')[[1]][1], 'PP2_')[[1]][2], FUN.VALUE = character(1)))

########################################################################################
PP1_tab = fread('Data/PP1_table.tsv')
PP1_ASV = PP1_tab$ASV
PP1_metadata = read.table('Metadata/PP1_metadata.tsv', header = T, sep = '\t')
PP1_EUCI_samples = as.vector(PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes'][PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes'] %in% colnames(PP1_tab)])
PP1_COMP_samples = as.vector(PP1_metadata$Sample[PP1_metadata$EUCI == 'No'][PP1_metadata$Sample[PP1_metadata$EUCI == 'No'] %in% colnames(PP1_tab)])
PP1_samples = c(PP1_EUCI_samples, PP1_COMP_samples)

PP2_tab = fread('Data/PP2_table.tsv')
PP2_ASV = PP2_tab$ASV
PP2_metadata = read.table('Metadata/PP2_metadata.tsv', header = T, sep = '\t')
PP2_EUCI_samples = as.vector(PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes'][PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes'] %in% colnames(PP2_tab)])
PP2_COMP_samples = as.vector(PP2_metadata$Sample[PP2_metadata$EUCI == 'No'][PP2_metadata$Sample[PP2_metadata$EUCI == 'No'] %in% colnames(PP2_tab)])
PP2_samples = c(PP2_EUCI_samples, PP2_COMP_samples)

PP1_tab = PP1_tab[, (PP1_samples) := lapply(.SD, function(col) col / sum(col)), .SDcols = PP1_samples]
PP2_tab = PP2_tab[, (PP2_samples) := lapply(.SD, function(col) col / sum(col)), .SDcols = PP2_samples]

# To pres/abs: 1/2000 abundance threshold: as we kept only samples with more than 2000 sequences
PP1_tab = PP1_tab[,..PP1_samples]
PP1_tab[PP1_tab > 0] = 1
PP2_tab = PP2_tab[,..PP2_samples]
PP2_tab[PP2_tab > 0] = 1

PP1_tab$EUCI_prevalence = rowSums(PP1_tab[,..PP1_EUCI_samples]) / length(PP1_EUCI_samples)
PP1_tab$COMP_prevalence = rowSums(PP1_tab[,..PP1_COMP_samples]) / length(PP1_COMP_samples)
PP1_tab$ASV = PP1_ASV
PP2_tab$EUCI_prevalence = rowSums(PP2_tab[,..PP2_EUCI_samples]) / length(PP2_EUCI_samples)
PP2_tab$COMP_prevalence = rowSums(PP2_tab[,..PP2_COMP_samples]) / length(PP2_COMP_samples)
PP2_tab$ASV = PP2_ASV

to_keep = c('ASV','EUCI_prevalence','COMP_prevalence')
PP1_preval = PP1_tab[,..to_keep]
PP2_preval = PP2_tab[,..to_keep]
PP1_EUCI_ASVs = PP1_preval$ASV[PP1_preval$EUCI_prevalence > 0]
PP1_COMP_ASVs = PP1_preval$ASV[PP1_preval$COMP_prevalence > 0]
PP2_EUCI_ASVs = PP2_preval$ASV[PP2_preval$EUCI_prevalence > 0]
PP2_COMP_ASVs = PP2_preval$ASV[PP2_preval$COMP_prevalence > 0]

########################################################################################
ProcessGenus <- function(aln_file, tree_file, Group1_ASVs, Group2_ASVs,max_iter=20){
    pairwise_mat = as.matrix(dist.alignment(read.alignment(aln_file, format='fasta'), matrix = "identity"))
    mean_pairwise_dist = mean(pairwise_mat, na.rm=T)**2
    median_pairwise_dist = median(pairwise_mat, na.rm=T)**2
    
    p_tree = read.tree(file = tree_file)
    n_G1_ASV = length(p_tree$tip.label[p_tree$tip.label %in% Group1_ASVs])
    n_G2_ASV = length(p_tree$tip.label[p_tree$tip.label %in% Group2_ASVs])
    
    if ((n_G1_ASV > 2) && (n_G2_ASV > 2)){
      p_tree = keep.tip(p_tree, p_tree$tip.label[p_tree$tip.label %in% c(Group1_ASVs, Group2_ASVs)])
      coph_tree = cophenetic(p_tree)
      coph_tree = coph_tree / max(coph_tree)
      diag(coph_tree) = NA
      
      nearest_intra_group = c()
      nearest_inter_group = c()
      
      for (curr_ASV in p_tree$tip.label[p_tree$tip.label %in% Group1_ASVs]){
        curr_intra = c()
        curr_inter = c()
        
        curr_G1_n = length(colnames(coph_tree)[colnames(coph_tree) %in% Group1_ASVs])
        curr_G2_n = length(colnames(coph_tree)[colnames(coph_tree) %in% Group2_ASVs])
        min_n = min(c(curr_G1_n,curr_G2_n))
        
        for (iteration in range(0,max_iter)){
          curr_G1_data = coph_tree[rownames(coph_tree)[rownames(coph_tree) == curr_ASV],sample(colnames(coph_tree)[colnames(coph_tree) %in% Group1_ASVs], min_n)]
          curr_G2_data = coph_tree[rownames(coph_tree)[rownames(coph_tree) == curr_ASV],sample(colnames(coph_tree)[colnames(coph_tree) %in% Group2_ASVs], min_n)]
          curr_intra = c(curr_intra, min(curr_G1_data, na.rm = T))
          curr_inter = c(curr_inter, min(curr_G2_data, na.rm = T))}
        nearest_intra_group = c(nearest_intra_group, mean(curr_intra, na.rm = T))
        nearest_inter_group = c(nearest_inter_group, mean(curr_inter, na.rm = T))}}
    else{nearest_intra_group = NA
         nearest_inter_group = NA}
    return(data.frame('Genus'=genus, 'Mean_pairwise_dist'=mean_pairwise_dist, 'Median_pairwise_dist'=median_pairwise_dist,
                      'G1_Mean_NN_EUCI'= mean(nearest_intra_group),'G2_Mean_NN_COMP'= mean(nearest_inter_group), 
                      'G1_Median_NN_EUCI'= median(nearest_intra_group),'G2_Median_NN_COMP'= median(nearest_inter_group), 
                      'n_G1_ASV'=n_G1_ASV, 'n_G2_ASV'=n_G2_ASV))}

PP1_MNPD_df <- foreach(genus=PP1_genera, .export = c('ProcessGenus'), .combine=rbind, .verbose = TRUE) %do% {
  aln_file = paste0('Data/Filtered/EUCI_genus_aln/PP1_',genus,'_aln.fasta')
  tree_file = paste0('Data/Filtered/EUCI_genus_trees/PP1_',genus,'_trim.fasta.treefile')
  ProcessGenus(aln_file, tree_file, PP1_EUCI_ASVs, PP1_COMP_ASVs, max_iter = 10)
}
write.csv(PP1_MNPD_df, file='5_Phylogenetics/PP1_MNNPD.csv')

PP2_MNPD_df <- foreach(genus=PP2_genera, .export = c('ProcessGenus'), .combine=rbind, .verbose = TRUE) %do% {
  aln_file = paste0('Data/Filtered/EUCI_genus_aln/PP2_',genus,'_aln.fasta')
  tree_file = paste0('Data/Filtered/EUCI_genus_trees/PP2_',genus,'_trim.fasta.treefile')
  ProcessGenus(aln_file, tree_file, PP2_EUCI_ASVs, PP2_COMP_ASVs, max_iter = 10)
}
write.csv(PP2_MNPD_df, file='5_Phylogenetics/PP2_MNNPD.csv')


PP1_MNPD_df = read.csv('5_Phylogenetics/PP1_MNNPD.csv')
PP2_MNPD_df = read.csv('5_Phylogenetics/PP2_MNNPD.csv')

strict_genera = read.csv('4_Taxonomic_analysis/strict_genera.csv')$x
relaxed_genera = read.csv('4_Taxonomic_analysis/relaxed_genera.csv')$x

PP1_MNPD_df$Group = 'Others'
PP1_MNPD_df$Group[PP1_MNPD_df$Genus %in% relaxed_genera] = 'Ancillary'
PP1_MNPD_df$Group[PP1_MNPD_df$Genus %in% strict_genera] = 'Core'
PP2_MNPD_df$Group = 'Others'
PP2_MNPD_df$Group[PP2_MNPD_df$Genus %in% relaxed_genera] = 'Ancillary'
PP2_MNPD_df$Group[PP2_MNPD_df$Genus %in% strict_genera] = 'Core'

PP1_MNPD_df$Group  <- factor(PP1_MNPD_df$Group , levels = c("Core","Ancillary", "Others"))
PP2_MNPD_df$Group  <- factor(PP2_MNPD_df$Group , levels = c("Core","Ancillary", "Others"))

# Mean nearest Phylogenetic distance to EUCI ASVs or Other ASVs
ggplot(PP1_MNPD_df[(PP1_MNPD_df$n_G1_ASV > 9) & (PP1_MNPD_df$n_G2_ASV > 9),], aes(x=G1_Mean_NN_EUCI, y=G2_Mean_NN_COMP, size=n_G1_ASV, colour=Group)) + 
  geom_point() + geom_abline(slope = 1,colour='grey20',lty=2) + 
  guides(size=guide_legend(title='EUCI ASV n.')) + xlim(c(0,0.4)) + ylim(c(0,0.4)) + 
  #geom_smooth(show.legend = F,colour = 'grey10', lty=5) +
  #geom_abline(intercept = 0.027335, slope = 1.223699 ,colour = 'grey10', lty=5, show.legend = F) +
  theme_bw() + xlab('MNNPD within EUCI') + ylab('MNNPD with non-EUCI') + scale_color_manual(values=c('#154194','#1F9DDB','grey'))
nrow(PP1_MNPD_df[(PP1_MNPD_df$n_G1_ASV > 9) & (PP1_MNPD_df$n_G2_ASV > 9),])
# n=558
ggsave('5_Phylogenetics/PP1_EUCI_MNPD.pdf', width = 6, height = 5)

ggplot(PP2_MNPD_df[(PP2_MNPD_df$n_G1_ASV > 9) & (PP2_MNPD_df$n_G2_ASV > 9),], aes(x=G1_Mean_NN_EUCI, y=G2_Mean_NN_COMP, size=n_G1_ASV, colour=Group)) +
  geom_point() + geom_abline(slope = 1,colour='grey20',lty=2) + 
  guides(size=guide_legend(title='EUCI ASV n.')) + xlim(0,0.4) + ylim(0,0.4) + 
  #geom_smooth(show.legend = F,colour = 'grey10', lty=5) +
  #geom_abline(intercept = 0.023815, slope = 1.038753 ,colour = 'grey10', lty=5, show.legend = F) +
  theme_bw() + xlab('MNNPD within EUCI') + ylab('MNNPD with non-EUCI') + scale_color_manual(values=c('#860600','#DB0600','grey'))
nrow(PP2_MNPD_df[(PP2_MNPD_df$n_G1_ASV > 9) & (PP2_MNPD_df$n_G2_ASV > 9),])
# 721
ggsave('5_Phylogenetics/PP2_EUCI_MNPD.pdf', width = 6, height = 5)

# Compare Pairwise Sequence Similarity
ggplot(PP1_MNPD_df) +
  geom_point(aes(y=1-Mean_pairwise_dist,x=n_G1_ASV,colour=Group)) + 
  xlab('EUCI ASV n.') + ylab('Mean pairwise similarity') + 
  xlim(0,1500) + ylim(0.75,1) +
  scale_color_manual(values=c('#154194','#1F9DDB','grey')) + theme_bw()
ggsave('5_Phylogenetics/PP1_pairwise_sim.pdf', width = 6, height = 5)

ggplot(PP1_MNPD_df) +
  geom_violin(aes(x=n_G1_ASV,y=Group,colour=Group)) + 
  ylab('') + xlab('EUCI ASV n.') + xlim(0,1500) + scale_x_log10() +
  scale_colour_manual(values=c('#154194','#1F9DDB','grey')) + theme_bw()
ggsave('5_Phylogenetics/PP1_EUCI_ASV_n.pdf', width = 3.5, height = 2)
wilcox.test(PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Core'], PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Others'])
# W = 77521, p-value < 2.2e-16
wilcox.test(PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Ancillary'], PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Others'])
# W = 487933, p-value < 2.2e-16
wilcox.test(PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Core'], PP1_MNPD_df$n_G1_ASV[PP1_MNPD_df$Group == 'Ancillary'])
# W = 25802, p-value = 4.042e-11

ggplot(PP2_MNPD_df) +
  geom_point(aes(y=1-Mean_pairwise_dist,x=n_G1_ASV,colour=Group)) + 
  xlab('EUCI ASV n.') + ylab('Mean pairwise similarity') + scale_x_log10() +
  xlim(0,1500) + ylim(0.75,1) +
  scale_color_manual(values=c('#860600','#DB0600','grey')) + theme_bw()
ggsave('5_Phylogenetics/PP2_pairwise_sim.pdf', width = 6, height = 5)

ggplot(PP2_MNPD_df) +
  geom_violin(aes(x=n_G1_ASV,y=Group,colour=Group)) + 
  ylab('') + xlab('EUCI ASV n.') + xlim(0,1500) + scale_x_log10() +
 scale_colour_manual(values=c('#860600','#DB0600','grey')) + theme_bw()
ggsave('5_Phylogenetics/PP2_EUCI_ASV_n.pdf', width = 3.5, height = 2)
wilcox.test(PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Core'], PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Others'])
# W = 73668, p-value < 2.2e-16
wilcox.test(PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Ancillary'], PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Others'])
# W = 393274, p-value < 2.2e-16
wilcox.test(PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Core'], PP2_MNPD_df$n_G1_ASV[PP2_MNPD_df$Group == 'Ancillary'])
# W = 23174, p-value = 1.008e-10



# Pairwise dist comparison
common_genera = intersect(PP1_MNPD_df$Genus, PP2_MNPD_df$Genus)
common_genera_df = data.frame('Genus'=c(), 'PP1_dist'=c(), 'PP2_dist'=c(), 'PP1_median_dist'=c(), 'PP2_median_dist'=c())
for (genus in common_genera){
  common_genera_df = rbind(common_genera_df, data.frame('Genus'=genus, 
                                                        'PP1_dist'=PP1_MNPD_df$Mean_pairwise_dist[PP1_MNPD_df$Genus == genus], 
                                                        'PP2_dist'=PP2_MNPD_df$Mean_pairwise_dist[PP2_MNPD_df$Genus == genus],
                                                        'PP1_median_dist'=PP1_MNPD_df$Median_pairwise_dist[PP1_MNPD_df$Genus == genus], 
                                                        'PP2_median_dist'=PP2_MNPD_df$Median_pairwise_dist[PP1_MNPD_df$Genus == genus]))
}

common_genera_df$Group = 'Others'
common_genera_df$Group[common_genera_df$Genus %in% relaxed_genera] = 'Ancillary'
common_genera_df$Group[common_genera_df$Genus %in% strict_genera] = 'Core'

ggplot(common_genera_df) + theme_bw() + geom_abline(intercept = 0, slope = 1, color='black', lty=2) + 
  geom_point(aes(x=PP1_dist,y=PP2_dist), colour='black', alpha = 0.2) + xlim(0,0.25) + ylim(0,0.25) +
  xlab('Mean pairwise nuc. distance (PP1)') + ylab('Mean pairwise nuc. distance (PP2)')
ggsave('5_Phylogenetics/Pairwise_distance_datasets_comparison.pdf', width = 5, height = 4) 

cor.test(common_genera_df$PP1_dist, common_genera_df$PP2_dist)
# t = 61.49, df = 1411, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.8385362 0.8669366
# sample estimates:
#   cor 
#   0.8533684 


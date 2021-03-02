library(ggplot2)
library(dplyr)
library(stringr)
library(circlize)

# Data loading
setwd('/Users/admin/Documents/Academia/PhD/Chapter I/')
stats_file = read.csv('Data/Unassigned_clusters_stats.tsv',sep='\t')
clusters_file = read.csv('9_Unassigned_genes/clusters_min_9_seq_2_samp.tsv',sep='\t')
metadata = read.csv('Metadata/MTG_metadata.tsv',sep='\t')

EUCI_samples = metadata$Sample[metadata$EUCI == 'Yes']

clusters_trna = read.csv('9_Unassigned_genes/cluster_tRNA_collection.txt',sep='\t', header = F)
clusters_trna$V1 = gsub('_tRNA.txt','',clusters_trna$V1)
clusters_trna$V1 = gsub('_','',clusters_trna$V1)

clusters_file$KEGG = str_split_fixed(clusters_file$Sequence, "_", 3)[,1]
clusters_file$Accession = str_split_fixed(clusters_file$Sequence, "_", 3)[,2]
clusters_file$SeqID =  str_split_fixed(clusters_file$Sequence, "_", 3)[,3]
stats_file$ClusterID = gsub('_','',stats_file$ClusterID)

stats_file$tRNA = 'No'
stats_file$EUCI = NA
stats_file$specificity = NA
stats_file$Unassigned = NA

for (cluster in unique(stats_file$ClusterID)){
  cluster = gsub("_", "", cluster)
  print(cluster)
  # tRNA
  if (cluster %in% clusters_trna$V1){stats_file$tRNA[stats_file$ClusterID == cluster] <- clusters_trna$V2[clusters_trna$V1 == cluster]}
  
  # Proportion of EUCI sequences
  sample_list = clusters_file$Accession[clusters_file$ClusterID == cluster]
  stats_file$specificity[stats_file$ClusterID == cluster] = sum(sample_list %in% EUCI_samples) / length(sample_list)
  if (stats_file$specificity[stats_file$ClusterID == cluster] == 1){stats_file$EUCI[stats_file$ClusterID == cluster] = 'EUCI'}
  else if (stats_file$specificity[stats_file$ClusterID == cluster] == 0){stats_file$EUCI[stats_file$ClusterID == cluster] = 'non-EUCI'}
  else {stats_file$EUCI[stats_file$ClusterID == cluster] <- 'Shared'}
  
  # Prop of assigned KEGG in cluster
  kegg_list = clusters_file$KEGG[clusters_file$ClusterID == cluster]
  if (sum(kegg_list %in% 'Unassigned') / length(kegg_list) < 1){stats_file$Unassigned[stats_file$ClusterID == cluster] <- 'KEGG'}
  else {stats_file$Unassigned[stats_file$ClusterID == cluster] <- 'Unassigned'}}

# A) tRNA analysis
library(ComplexHeatmap)
euci_trna = data.frame(as.list(table(stats_file$tRNA[(stats_file$EUCI == 'EUCI') & (stats_file$tRNA != 'No')])))
noneuci_trna = data.frame(as.list(table(stats_file$tRNA[(stats_file$EUCI == 'non-EUCI') & (stats_file$tRNA != 'No')])))
shared_trna = table(stats_file$tRNA[(stats_file$EUCI == 'Shared') & (stats_file$tRNA != 'No')])
shared_trna['Undet'] = 0
shared_trna = data.frame(as.list(shared_trna))
full_trna = rbind('EUCI' = euci_trna, 'non-EUCI' = noneuci_trna, 'Shared' = shared_trna)
full_trna = as.data.frame(t(as.matrix(full_trna / rowSums(full_trna))))
full_trna$diff_euci = full_trna$EUCI - full_trna$`non-EUCI`
full_trna$row_split = 'Positive'
full_trna$row_split[full_trna$diff_euci < 0] = 'Negative'
full_trna$row_split[rownames(full_trna) == 'Undet'] = 'Undet'

pdf('9_Unassigned_genes/tRNA_summary.pdf', width = 4, height = 6.6)
f1 = colorRamp2(seq(min(as.matrix(full_trna[,c('EUCI','non-EUCI','Shared')])), max(as.matrix(full_trna[,c('EUCI','non-EUCI','Shared')])), length = 3), c('#23A671','white','#A673D1'))
Heatmap(as.matrix(full_trna[,c('EUCI','non-EUCI','Shared')]), col = f1, show_column_dend = FALSE, cluster_rows = FALSE, name = 'Frequency', right_annotation =  rowAnnotation('Diff' = anno_barplot(full_trna$diff_euci)))
dev.off()

trna_data = stats_file[stats_file$tRNA != 'No',]
ggplot(trna_data, aes(x=EUCI,y=meanGC,colour=EUCI)) + geom_boxplot() + theme_minimal() + 
  scale_color_manual(values = c('#23A671','dimgrey','#A673D1')) + ylab('GC content') + xlab('') + theme(legend.position = "none")
ggsave('9_Unassigned_genes/tRNA_gc_content.pdf', width = 3, height = 2.2)

ggplot(trna_data, aes(x=EUCI,y=100-meanPID,colour=EUCI)) + geom_boxplot() + theme_minimal() + scale_y_log10() +
  scale_color_manual(values = c('#23A671','dimgrey','#A673D1')) + ylab('Divergence') + xlab('') + theme(legend.position = "none")
ggsave('9_Unassigned_genes/tRNA_divergence.pdf', width = 3, height = 2.2)

ggplot(trna_data, aes(x=EUCI,y=meanSeqLen,colour=EUCI)) + geom_boxplot() + theme_minimal() + scale_y_log10() + ylim(c(60,100)) +
  scale_color_manual(values = c('#23A671','dimgrey','#A673D1')) + ylab('Length(bp)') + xlab('') + theme(legend.position = "none")
ggsave('9_Unassigned_genes/tRNA_length.pdf', width = 3, height = 2.2)

# B) EUCI only analysis
euci_data = stats_file[(stats_file$tRNA == 'No') & (stats_file$EUCI == 'EUCI'),]
ggplot(euci_data, aes(size=SeqNum, x=meanGC, y=100-meanPID, colour=Unassigned)) + 
  geom_point(alpha=0.5) + scale_y_log10() + theme_minimal() + ylab('Divergence') + xlab('GC content')
ggsave('9_Unassigned_genes/EUCI_clusters.pdf', width = 6, height = 5)

# C) Stats comparison with others
comp_df = data.frame('Ecosystem' = c(), 'Category' = c(), 'Cluster_number' = c())
for (i in c('EUCI','non-EUCI','Shared')){
    trna_n = nrow(stats_file[(stats_file$EUCI == i) & (stats_file$tRNA != 'No'),])
    kegg_n = table(stats_file$Unassigned[(stats_file$EUCI == i) & (stats_file$tRNA == 'No')])['KEGG']
    unass_n = table(stats_file$Unassigned[(stats_file$EUCI == i) & (stats_file$tRNA == 'No')])['Unassigned']
    comp_df = rbind(comp_df, data.frame('Ecosystem' = c(i,i,i), 'Category' = c('tRNA','KEGG','Unassigned'), 'Cluster_number' = c(trna_n, kegg_n, unass_n)))
    }
ggplot(comp_df, aes(x=Category, y=Cluster_number, fill=Ecosystem)) + xlab('') + scale_fill_manual(values = c('#23A671','grey','#A673D1')) +
  geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + scale_y_log10() + ylab('Cluster number')
ggsave('9_Unassigned_genes/Unassigned_proportions.pdf', width = 4, height = 3)
# n EUCI = 9492
# n non-EUCI = 10129
# n Shared = 1341

# '#23A671','#14B8CE','#91B9CF'
#f2 = colorRamp2(seq(min(comp_df), max(comp_df), length = 2), c('white','#23A671'))
#pdf('9_Unassigned_genes/Unassigned_proportions.pdf', width = 4, height = 3)
#Heatmap(comp_df, name = 'Proportion', col = f2, show_column_dend = FALSE, cluster_rows = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
#  grid.text(sprintf("%.2f", comp_df[i, j]), x, y, gp = gpar(fontsize = 10))
#})
#dev.off()

summary(stats_file$meanGC[(stats_file$EUCI == 'EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanGC[(stats_file$EUCI == 'non-EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanGC[(stats_file$EUCI == 'Shared') & (stats_file$tRNA == 'No')])
ggplot(stats_file, aes(x=Unassigned, colour=EUCI, y=meanGC)) + geom_boxplot()  + xlab('') + ylab('GC content') +  theme(legend.position = "none") +
  theme_minimal() + scale_color_manual(values = c('#23A671','dimgrey','#A673D1'))
ggsave('9_Unassigned_genes/Unassigned_gc_content.pdf', width = 3, height = 2)

summary(stats_file$meanPID[(stats_file$EUCI == 'EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanPID[(stats_file$EUCI == 'non-EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanPID[(stats_file$EUCI == 'Shared') & (stats_file$tRNA == 'No')])
ggplot(stats_file, aes(x=Unassigned, colour=EUCI, y=100-meanPID)) + geom_boxplot() + xlab('') + ylab('Divergence') +
  theme_minimal() + scale_y_log10() + scale_color_manual(values = c('#23A671','dimgrey','#A673D1'))
ggsave('9_Unassigned_genes/Unassigned_PID.pdf', width = 3, height = 2)

summary(stats_file$meanSeqLen[(stats_file$EUCI == 'EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanSeqLen[(stats_file$EUCI == 'non-EUCI') & (stats_file$tRNA == 'No')])
summary(stats_file$meanSeqLen[(stats_file$EUCI == 'Shared') & (stats_file$tRNA == 'No')])
ggplot(stats_file, aes(x=Unassigned, colour=EUCI, y=meanSeqLen)) + geom_boxplot() + scale_y_log10() + theme_minimal() + scale_color_manual(values = c('#23A671','dimgrey','#A673D1'))


stats_file = na.omit(stats_file)
euci_only = stats_file[stats_file$EUCI == 'EUCI',]
write.csv(stats_file, file='9_Unassigned_genes/Clusters_stats.csv')


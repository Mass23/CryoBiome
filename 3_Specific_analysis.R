library(DESeq2)
library(EnhancedVolcano)

# Colors light:
#---------------
# PP1: #1F9DDB
# PP2: #DB0600
# MTG: #6BC90C
  
# Colors dark:
#--------------
# PP1: #1A87BD
# PP2: #BA0600
# MTG: #5CAD0A

setwd('/Users/admin/Documents/Academia/PhD/Chapter I/')
############################################################################################################
# Loading data

PP1_genus_tab = read.csv('Data/PP1_genus.csv')
PP2_genus_tab = read.csv('Data/PP2_genus.csv')
MTG_genus_tab = read.csv('Data/MTG_genus.csv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')

# for PP2, remove the snow study: bad sequences (not assigned further than bacteria, blast: 87% with 16s seqs.)
PP1_data = PP1_genus_tab[,PP1_metadata$Sample[PP1_metadata$Sample %in% colnames(PP1_genus_tab)]]
rownames(PP1_data) = PP1_genus_tab$Genus

PP2_data = PP2_genus_tab[,PP2_metadata$Sample[PP2_metadata$Sample %in% colnames(PP2_genus_tab)]]
rownames(PP2_data) = PP2_genus_tab$Genus

MTG_data = MTG_genus_tab[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_genus_tab)]]
rownames(MTG_data) = MTG_genus_tab$Genus

############################################################################################################
# DESeq2 analysis

PP1_dds <- DESeqDataSetFromMatrix(countData=PP1_data+1, colData=PP1_metadata, design=~EUCI)
PP1_deseq <- DESeq(PP1_dds)
PP1_res <- results(PP1_deseq)
PP1_res$padj[is.na(PP1_res$padj)] = 1
PP1_sign_genera = rownames(PP1_res)[(PP1_res$padj < 0.05) & (PP1_res$log2FoldChange > 1)]

PP2_dds <- DESeqDataSetFromMatrix(countData=PP2_data+1, colData=PP2_metadata, design=~EUCI)
PP2_deseq <- DESeq(PP2_dds)
PP2_res <- results(PP2_deseq)
PP2_res$padj[is.na(PP2_res$padj)] = 1
PP2_sign_genera = rownames(PP2_res)[(PP2_res$padj < 0.05) & (PP2_res$log2FoldChange > 1)]

MTG_dds <- DESeqDataSetFromMatrix(countData=MTG_data+1, colData=MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_genus_tab),], design=~EUCI)
MTG_deseq <- DESeq(MTG_dds)
MTG_res <- results(MTG_deseq)
MTG_res$padj[is.na(MTG_res$padj)] = 1
MTG_sign_genera = rownames(MTG_res)[(MTG_res$padj < 0.05) & (MTG_res$log2FoldChange > 1)]

EnhancedVolcano(PP1_res,
                lab = rownames(PP1_res),
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey30", "grey30", "#1F9DDB"),
                x = 'log2FoldChange',
                y = 'pvalue')
ggsave('3_Specific_analysis/PP1_specific_genera.pdf', width = 6, height = 7)

EnhancedVolcano(PP2_res,
                lab = rownames(PP2_res),
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey30", "grey30", "#DB0600"),
                x = 'log2FoldChange',
                y = 'pvalue')
ggsave('3_Specific_analysis/PP2_specific_genera.pdf', width = 6, height = 7)

EnhancedVolcano(MTG_res,
                lab = rownames(MTG_res),
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey30", "grey30", "#6BC90C"),
                x = 'log2FoldChange',
                y = 'pvalue')
ggsave('3_Specific_analysis/MTG_specific_genera.pdf', width = 6, height = 7)

write.csv(PP1_res, '3_Specific_analysis/PP1_deseq_res.csv')
write.csv(PP2_res, '3_Specific_analysis/PP2_deseq_res.csv')
write.csv(MTG_res, '3_Specific_analysis/MTG_deseq_res.csv')

length(PP1_sign_genera)
# 452
length(PP2_sign_genera)
# 292
length(MTG_sign_genera)
# 59



both_sign_genera = intersect(PP1_sign_genera, PP2_sign_genera)
# 102
write.csv(PP1_sign_genera, '3_Specific_analysis/PP1_significant_genera.csv')
write.csv(PP2_sign_genera, '3_Specific_analysis/PP2_significant_genera.csv')
write.csv(both_sign_genera, '3_Specific_analysis/Both_significant_genera.csv')
write.csv(MTG_sign_genera, '3_Specific_analysis/MTG_significant_genera.csv')

nrow(PP1_res)
# 3101
nrow(PP2_res)
# 3012
length(union(rownames(PP1_res), rownames(PP2_res)))
# 3441

############################################################################################################
# Comparison of PP1 with PP2
Fold_change_comp = data.frame(Genus = intersect(rownames(PP1_res), rownames(PP2_res)))
Fold_change_comp$PP1_fc = vapply(Fold_change_comp$Genus, function(x) PP1_res$log2FoldChange[rownames(PP1_res) == x], FUN.VALUE = numeric(1))
Fold_change_comp$PP2_fc = vapply(Fold_change_comp$Genus, function(x) PP2_res$log2FoldChange[rownames(PP2_res) == x], FUN.VALUE = numeric(1))
Fold_change_comp$Group = 'Others'
Fold_change_comp$Group[Fold_change_comp$Genus %in% PP1_sign_genera] = 'PP1 sign.'
Fold_change_comp$Group[Fold_change_comp$Genus %in% PP2_sign_genera] = 'PP2 sign.'
Fold_change_comp$Group[Fold_change_comp$Genus %in% both_sign_genera] = 'Both sign.'

cor = cor.test(Fold_change_comp$PP1_fc, Fold_change_comp$PP2_fc)
cor_p = cor$p.value
cor_st = cor$statistic
cor_est = cor$estimate

ggplot(Fold_change_comp) +
  geom_point(aes(x = PP1_fc, y = PP2_fc, color = Group), alpha =0.5) + 
  scale_color_manual(values = c('black', 'grey', '#1A87BD', '#BA0600')) +
  xlab('PP1 log2 Fold change') +
  ylab('PP2 log2 Fold change') + 
  ggtitle(paste0("Pearson's corr: ",round(cor_est,3), ", p: ",round(cor_p,3), ", statistic: ", round(cor_st,3))) +
  xlim(c(-10,10)) +
  ylim(c(-7,7)) +
  geom_hline(yintercept = 0, color = 'grey40', lty = 2) +
  geom_vline(xintercept = 0, color = 'grey40', lty = 2) +
  theme_bw()
ggsave('3_Specific_analysis/Comp_genera.pdf', width = 5, height = 5)



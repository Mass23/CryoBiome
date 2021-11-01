library(data.table)
library(ggplot2)

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

###################################################################
ancom_table = read.csv('Data/Amplicon_ancom_res.csv')
ancom_table$DA = as.character(ancom_table$DA)
ancom_table$DA[(ancom_table$CLR_mean_diff < 0 ) & (ancom_table$W > 1251)] = 'Underrepresented'
ancom_table$Genus = vapply(ancom_table$taxa_id, function(x) strsplit(as.character(x),split = 'g__') [[1]][2],FUN.VALUE = character(1))

refseq_table = read.csv('Data/prokaryotes.txt',sep='\t')
cub_table = read.csv('Data/merged_all_growth_prediction.txt',sep='\t')
refseq_table = refseq_table[refseq_table$Status == 'Complete Genome',]
refseq_table$d = NA
refseq_table$d_lowerCI = NA
refseq_table$d_upperCI = NA
refseq_table$CUBHE = NA
refseq_table$ConsistencyHE = NA
refseq_table$CPB = NA

# Add the CUB information to the refseq table
for (accession in refseq_table$Assembly.Accession){
  cub_data = cub_table[grep(accession,cub_table$Sample),]
  if (nrow(cub_data) == 1){
    refseq_table$d[refseq_table$Assembly.Accession == accession] = cub_data$d
    refseq_table$CUBHE[refseq_table$Assembly.Accession == accession] = cub_data$CUBHE
    refseq_table$ConsistencyHE[refseq_table$Assembly.Accession == accession] = cub_data$ConsistencyHE
    refseq_table$CPB[refseq_table$Assembly.Accession == accession] = cub_data$CPB
  }
}

# Merge deseq and refseq information in another table
ancom_refseq_table = data.frame(Genus=c(),DA=c(),clr_diff=c(),mean_GC=c(),mean_size_mb=c(),mean_gene_number=c(),n_genomes=c(),
                                d=c(),CUBHE=c(),ConsistencyHE=c(),CPB=c())

refseq_table$Genes = as.integer(as.character(refseq_table$Genes))
refseq_table$GC. = as.numeric(as.character(refseq_table$GC.))

CleanName <- function(genus_name){
  genus_split = unlist(strsplit(as.character(genus_name), split = ' '))
  if (genus_split[1] == 'Candidatus'){genus_name = paste(genus_split[1], genus_split[2],collapse = ' ')}
  else{genus_name = genus_split[1]}
  if (unlist(gregexpr("[A-Z]", genus_name)) == 1){
    print(genus_name)
    return(genus_name)}
  else{return('')}}

refseq_table$clean_genus = vapply(refseq_table$X.Organism.Name, function(x)CleanName(x), character(1))

for (genus in ancom_table$Genus){
  genus_ncbi = gsub("_", " ", genus)
  table_subset = refseq_table[refseq_table$clean_genus %like% genus_ncbi,]

  if (nrow(table_subset) > 0){
    ancom_refseq_table = rbind(ancom_refseq_table, data.frame(Genus = genus, 
                                                              Taxonomy = ancom_table$taxa_id[ancom_table$Genus == genus],
                                                              DA = ancom_table$DA[ancom_table$Genus == genus], 
                                                              clr_diff = ancom_table$CLR_mean_diff[ancom_table$Genus == genus], 
                                                              mean_GC = mean(table_subset$GC.),
                                                              mean_size_mb = mean(table_subset$Size..Mb.),
                                                              mean_gene_number = mean(as.numeric(table_subset$Genes)),
                                                              n_genomes = c(nrow(table_subset)),
                                                              d=mean(table_subset$d),
                                                              CUBHE=mean(table_subset$CUBHE),
                                                              ConsistencyHE=mean(table_subset$ConsistencyHE),
                                                              CPB=mean(table_subset$CPB)))}
  else{ancom_refseq_table = rbind(ancom_refseq_table, data.frame(Genus = genus, 
                                                                 Taxonomy = ancom_table$taxa_id[ancom_table$Genus == genus],
                                                                 DA = ancom_table$DA[ancom_table$Genus == genus], 
                                                                 clr_diff = ancom_table$CLR_mean_diff[ancom_table$Genus == genus], 
                                                                 mean_GC = NA,
                                                                 mean_size_mb = NA,
                                                                 mean_gene_number = NA,
                                                                 n_genomes = 0,
                                                                 d=NA,
                                                                 CUBHE=NA,
                                                                 ConsistencyHE=NA,
                                                                 CPB=NA))}}

write.csv(ancom_refseq_table, file = 'Data/ancom_refseq_genomes.csv', row.names = FALSE)
###################################################################
library(ggpubr)
library(performance)
library(lme4)
library(see)
library(dplyr)
ancom_refseq_table = read.csv('Data/ancom_refseq_genomes.csv')

# Mean GC
gc_df = as.data.frame(na.omit(ancom_refseq_table[,c('mean_GC','clr_diff','DA')]))
gc_df$DA = as.character(gc_df$DA)
gc_df$DA[gc_df$DA == 'Overrepresented'] = 'Cryo.'
compare_means(data=gc_df, mean_GC ~ DA, ref.group='Others')
# .y.       group1  group2        p  p.adj p.format p.signif method  
# mean_GC   Others  Cryo.    0.0122  0.012 0.012    *        Wilcoxon
gc_df %>% group_by(DA) %>% summarise(median(mean_GC)) # Others = 52.8, Over = 56.6
ggplot(gc_df,aes(x=DA,y=mean_GC,fill=DA)) + geom_boxplot() + theme_linedraw() + ylab('GC content [%]') + xlab('') + theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", ref.group='Others') + scale_fill_manual(values=c('#3C5488FF','#F39B7FFF'))
ggsave('2_Genus_analysis/2_2_Refseq/2_2_Refseq_GC.pdf', width=2, height=3)

# Genome size
size_df = as.data.frame(na.omit(ancom_refseq_table[,c('mean_size_mb','clr_diff','DA')]))
size_df$DA = as.character(size_df$DA)
size_df$DA[size_df$DA == 'Overrepresented'] = 'Cryo.'
compare_means(data=size_df, mean_size_mb ~ DA, ref.group='Others')
# .y.          group1  group2        p  p.adj  p.format p.signif method  
# mean_size_mb Others  Cryo.   0.00524  0.0052 0.0052   **       Wilcoxon
size_df %>% group_by(DA) %>% summarise(median(mean_size_mb)) # Others = 3.96, Over = 4.28
ggplot(size_df,aes(x=DA,y=mean_size_mb,fill=DA)) + geom_boxplot() + theme_linedraw() + ylab('Genome size [mbp]') + xlab('') + theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", ref.group='Others') + scale_fill_manual(values=c('#3C5488FF','#F39B7FFF')) + scale_y_log10()
ggsave('2_Genus_analysis/2_2_Refseq/2_2_Refseq_size.pdf', width=2, height=3)

# Growth doubling time
d_df = as.data.frame(na.omit(ancom_refseq_table[,c('d','clr_diff','DA')]))
d_df$DA = as.character(d_df$DA)
d_df$DA[d_df$DA == 'Overrepresented'] = 'Cryo.'
d_df %>% group_by(DA) %>% summarise(median(d)) 
compare_means(data=d_df, d ~ DA, ref.group='Others')
# .y.   group1  group2     p p.adj p.format p.signif method  
# d     Others Cryo.   0.145  0.15 0.15     ns       Wilcoxon

# Codon Usage Bias
cub_df = as.data.frame(na.omit(ancom_refseq_table[,c('CUBHE','clr_diff','DA')]))
cub_df$DA = as.character(cub_df$DA)
cub_df$DA[cub_df$DA == 'Overrepresented'] = 'Cryo.'
cub_df %>% group_by(DA) %>% summarise(median(CUBHE)) 
compare_means(data=cub_df, CUBHE ~ DA, ref.group='Others')
# .y.   group1 group2     p  p.adj p.format p.signif method  
# CUBHE Others Cryo.  0.978  0.98  0.98     ns       Wilcoxon


# Consistency HE
con_df = as.data.frame(na.omit(ancom_refseq_table[,c('ConsistencyHE','clr_diff','DA')]))
con_df$DA = as.character(con_df$DA)
con_df$DA[con_df$DA == 'Overrepresented'] = 'Cryo.'
con_df %>% group_by(DA) %>% summarise(median(ConsistencyHE)) 
compare_means(data=con_df, ConsistencyHE ~ DA, ref.group='Others')
# .y.           group1 group2      p p.adj p.format p.signif method  
# ConsistencyHE Others Cryo.  0.0415 0.042 0.042    *        Wilcoxon
ggplot(con_df,aes(x=DA,y=ConsistencyHE,fill=DA)) + geom_boxplot() + theme_linedraw() + ylab('Consistency HE') + xlab('') + theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", ref.group='Others') + scale_fill_manual(values=c('#3C5488FF','#F39B7FFF')) 
ggsave('2_Genus_analysis/2_2_Refseq/2_2_Refseq_ConsistencyHE.pdf', width=2, height=3)

# CPB
cpb_df = as.data.frame(na.omit(ancom_refseq_table[,c('CPB','clr_diff','DA')]))
cpb_df$DA = as.character(cpb_df$DA)
cpb_df$DA[cpb_df$DA == 'Overrepresented'] = 'Over.'
cpb_df %>% group_by(DA) %>% summarise(median(CPB)) 
compare_means(data=cpb_df, CPB ~ DA, ref.group='Others')
# .y.   group1 group2     p p.adj p.format p.signif method  
# CPB   Others Over.  0.410  0.41 0.41     ns       Wilcoxon


###################################################################
# 2. Codon / Amino Acid Composition
library(DESeq2)
library(dplyr)
setwd('/Users/mabourqu/Documents/PhD/C1/')

codon_counts = read.csv('Data/merged_all_codon_counts.txt',sep='\t')
refseq_table = read.csv('Data/prokaryotes.txt',sep='\t')
ancom_table = read.csv('Data/Amplicon_ancom_res.csv')
ancom_table$DA = as.character(ancom_table$DA)
ancom_table$DA[(ancom_table$CLR_mean_diff < 0 ) & (ancom_table$W > 1251)] = 'Others'
ancom_table$DA[ancom_table$DA == 'Overrepresented'] = 'Cryo.'
ancom_table$Genus = vapply(ancom_table$taxa_id, function(x) as.character(strsplit(as.character(x),split = 'g__') [[1]][2]),FUN.VALUE = character(1))

# Aggregate by genome
codon_counts = aggregate(codon_counts[,colnames(codon_counts) != 'Sample'], by=list(Sample=codon_counts$Sample), FUN=sum)
# Add genus information
codon_counts$Sample = vapply(as.character(codon_counts$Sample), function(x){strsplit(x,split = '_ASM')[[1]][1]}, FUN.VALUE = character(1))
refseq_table$Genus = vapply(as.character(refseq_table$X.Organism.Name), function(x) strsplit(x,split = ' ')[[1]][1], FUN.VALUE = character(1))
codon_counts$Genus = vapply(as.character(codon_counts$Sample), function(x){ifelse(x %in% refseq_table$Assembly.Accession, refseq_table$Genus[refseq_table$Assembly.Accession == x], 'NA')}, FUN.VALUE = character(1))
# Calculate frequencies, aggregate by genus
codon_counts[,!names(codon_counts) %in% c('Genus','Sample')] = codon_counts[,!names(codon_counts) %in% c('Genus','Sample')] / rowSums(codon_counts[,!names(codon_counts) %in% c('Genus','Sample')])
codon_counts = aggregate(codon_counts[,!names(codon_counts) %in% c('Sample','Genus')], by=list(Genus=codon_counts$Genus), FUN=mean)
rownames(codon_counts) = codon_counts$Genus
codon_counts$Genus = NULL

freq_comp = expand.grid(Taxonomy=unique(ancom_table$taxa_id), Codon=colnames(codon_counts))
freq_comp$Genus = vapply(freq_comp$Taxonomy, function(x) as.character(strsplit(as.character(x),split = 'g__') [[1]][2]),FUN.VALUE = character(1))
freq_comp = freq_comp[freq_comp$Genus %in% rownames(codon_counts),]
freq_comp$GC_codon = (lengths(regmatches(freq_comp$Codon, gregexpr("G", freq_comp$Codon))) + lengths(regmatches(freq_comp$Codon, gregexpr("C", freq_comp$Codon)))) / 3
freq_comp$Frequency = vapply(1:nrow(freq_comp), function(x){codon_counts[rownames(codon_counts) == freq_comp$Genus[x], freq_comp$Codon[x]]}, FUN.VALUE = numeric(1))
freq_comp$gc_freq = freq_comp$GC_codon * freq_comp$Frequency

# Get aggregates by genus
library(ggpubr)
genus_gc_genes = freq_comp %>% group_by(Genus) %>% summarise(GC=sum(gc_freq))
genus_gc_genes$DA = vapply(genus_gc_genes$Genus, function(x) ancom_table$DA[ancom_table$Genus==x], FUN.VALUE = character(1))
genus_gc_genes$Genome_GC = vapply(genus_gc_genes$Genus, function(x) mean(refseq_table$GC.[refseq_table$Genus == x]), FUN.VALUE = numeric(1))
compare_means(data=genus_gc_genes, GC ~ DA, ref.group='Others')
# .y.   group1  group2      p   p.adj p.format p.signif method  
# GC    Others Cryo.  0.0219 0.022 0.022    *        Wilcoxon
genus_gc_genes %>% group_by(DA) %>% summarise(median(GC)) # Others=0.534, Cryo.=0.584
ggplot(genus_gc_genes, aes(y=GC*100,x=DA,fill=DA)) + geom_boxplot() + theme_linedraw() + ylab('Genes GC content [%]') + xlab('') + theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", ref.group='Others') + scale_fill_manual(values=c('#3C5488FF','#F39B7FFF'))
ggsave('2_Genus_analysis/2_2_Refseq/2_2_Genes_GC.pdf', width=2, height=3)

###############################################################
# Amino acids
library(Biostrings)
aa_counts_ints = as.data.frame(t(as.data.frame(codon_counts)))
aa_counts_ints$AminoAcid = vapply(rownames(aa_counts_ints), function(x){GENETIC_CODE[[x]]}, FUN.VALUE = character(1))
aa_counts_ints = aggregate(aa_counts_ints[,!names(aa_counts_ints) %in% c('AminoAcid')], by=list(AA=aa_counts_ints$AminoAcid), FUN=sum)
colnames(aa_counts_ints) = vapply(colnames(aa_counts_ints), function(x) gsub('\\[|\\]','',x), FUN.VALUE = character(1))
colnames(aa_counts_ints) = vapply(colnames(aa_counts_ints), function(x) gsub("'",'',x), FUN.VALUE = character(1))
AA_list = aa_counts_ints$AA
aa_counts_ints$AA = NULL
aa_counts_ints = aa_counts_ints * 100000
aa_counts_ints = as.data.frame(sapply(aa_counts_ints, as.integer))
aa_counts_ints[is.na(aa_counts_ints)] = 0

metadata = data.frame(Genus = colnames(aa_counts_ints))
metadata$Category = 'Others'
metadata$Category[metadata$Genus %in% ancom_table$Genus[ancom_table$DA == 'Cryo.']] = 'Cryo.'

aa_dds <- DESeqDataSetFromMatrix(countData=aa_counts_ints, colData=metadata, design=~Category)
aa_dds$Category = factor(aa_dds$Category, levels = c("Others","Cryo."))
aa_deseq <- DESeq(aa_dds)
aa_res <- results(aa_deseq)
aa_res$padj[is.na(aa_res$padj)] = 1
aa_res$AA = AA_list
aa_res$POS = FALSE
aa_res$POS[aa_res$log2FoldChange > 0] = TRUE

ggplot(as.data.frame(aa_res),aes(x=reorder(AA, log2FoldChange),y=log2FoldChange)) + scale_fill_manual(values=c('#F39B7FFF','#3C5488FF')) + 
  geom_bar(stat='identity',aes(fill=POS)) + coord_flip() + theme_linedraw() + xlab('') + ylab('Log2 Fold Change') + theme(legend.position = "none")
ggsave('2_Genus_analysis/2_2_Refseq/2_2_Amino_acid_comparison.pdf', width = 6, height = 5)


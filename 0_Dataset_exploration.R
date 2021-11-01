library(data.table)
library(ggplot2)
setwd('/Users/mabourqu/Documents/PhD/C1/')

#####################################################################
# Create abundance informed fasta files
library(seqinr)
# PP1
PP1_tab = fread('Data/PP1_table.tsv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep ='\t')
PP1_samples = as.vector(PP1_metadata$Sample)
# get colsums
colsums = PP1_tab[, lapply(.SD, sum), .SDcols = as.character(PP1_metadata$Sample)]
max_count = max(colsums)
# multiply by factor (based on max value in columns)
for (sample in as.character(PP1_metadata$Sample)){
  factor = max_count / colsums[[sample]]
  print(factor)
  PP1_tab = set(PP1_tab, i=NULL, j=sample, value=PP1_tab[[sample]]*factor)}
# verify the colsums are well normalised
norm_colsums = as.numeric(PP1_tab[, lapply(.SD, sum), .SDcols = as.character(PP1_metadata$Sample)])
mean(norm_colsums)
sd(norm_colsums)
# calculate the sum for the fasta files, rename the fasta entries
PP1_tab[, `:=`(ASV_sum = as.integer(rowSums(.SD))), .SDcols=as.character(PP1_metadata$Sample)]
PP1_fasta = read.fasta(file='Data/trees/PP1_ASV_seqs.fasta')
sum(names(PP1_fasta) == PP1_tab$ASV) == nrow(PP1_tab)
PP1_tab = PP1_tab[, abund_name := paste(ASV, ASV_sum, sep = "_")]
names(PP1_fasta) = PP1_tab$abund_name
write.fasta(PP1_fasta, names = names(PP1_fasta), file.out='Data/OTU_clustering/PP1_abundance_ASV_seqs.fasta')

###############################
# PP2
# Create abundance informed fasta files
PP2_tab = fread('Data/PP2_table.tsv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep ='\t')
PP2_samples = as.vector(PP2_metadata$Sample)
# get colsums
colsums = PP2_tab[, lapply(.SD, sum), .SDcols = as.character(PP2_metadata$Sample)]
max_count = max(colsums)
# multiply by factor (based on max value in columns)
for (sample in as.character(PP2_metadata$Sample)){
  factor = max_count / colsums[[sample]]
  print(factor)
  PP2_tab = set(PP2_tab, i=NULL, j=sample, value=PP2_tab[[sample]]*factor)}
# verify the colsums are well normalised
norm_colsums = as.numeric(PP2_tab[, lapply(.SD, sum), .SDcols = as.character(PP2_metadata$Sample)])
mean(norm_colsums)
sd(norm_colsums)
# calculate the sum for the fasta files, rename the fasta entries
PP2_tab[, `:=`(ASV_sum = as.integer(rowSums(.SD))), .SDcols=as.character(PP2_metadata$Sample)]
PP2_fasta = read.fasta(file='Data/trees/PP2_ASV_seqs.fasta')
sum(names(PP2_fasta) == PP2_tab$ASV) == nrow(PP2_tab)
PP2_tab = PP2_tab[, abund_name := paste(ASV, ASV_sum, sep = "_")]
names(PP2_fasta) = PP2_tab$abund_name
write.fasta(PP2_fasta, names = names(PP2_fasta), file.out='Data/OTU_clustering/PP2_abundance_ASV_seqs.fasta')

#####################################################################
# Create most abundant ASVs dataset
library(seqinr)
# PP1
PP1_tab = fread('Data/PP1_table.tsv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep ='\t')
PP1_samples = as.vector(PP1_metadata$Sample)
PP1_tab = PP1_tab[, (PP1_samples) := lapply(.SD, function(x){x/sum(x)}), .SDcols = PP1_samples]
PP1_tab = PP1_tab[, (PP1_samples) := lapply(.SD, function(x){ifelse(x > 0.005, 1, 0)}), .SDcols = PP1_samples]
PP1_tab_abund = as.data.frame(PP1_tab[rowSums(PP1_tab[,..PP1_samples]) > 0, ])
write.csv(PP1_tab_abund, file='Data/trees/PP1_005_ASVs_table.csv',quote = FALSE,row.names = FALSE)

PP1_fasta = read.fasta(file='Data/trees/PP1_ASV_seqs.fasta')
PP1_subset = PP1_fasta[names(PP1_fasta) %in% PP1_tab_abund$ASV]
write.fasta(PP1_subset, names = names(PP1_subset), file.out='Data/trees/PP1_005_ASV_seqs.fasta')

# PP2
PP2_tab = fread('Data/PP2_table.tsv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep ='\t')
PP2_samples = as.vector(PP2_metadata$Sample)
PP2_tab = PP2_tab[, (PP2_samples) := lapply(.SD, function(x){x/sum(x)}), .SDcols = PP2_samples]
PP2_tab = PP2_tab[, (PP2_samples) := lapply(.SD, function(x){ifelse(x > 0.005, 1, 0)}), .SDcols = PP2_samples]
PP2_tab_abund = as.data.frame(PP2_tab[rowSums(PP2_tab[,..PP2_samples]) > 0, ])
write.csv(PP2_tab_abund, file='Data/trees/PP2_005_ASVs_table.csv',quote = FALSE,row.names = FALSE)

PP2_fasta = read.fasta(file='Data/trees/PP2_ASV_seqs.fasta')
PP2_subset = PP2_fasta[names(PP2_fasta) %in% PP2_tab_abund$ASV]
write.fasta(PP2_subset, names = names(PP2_subset), file.out='Data/trees/PP2_005_ASV_seqs.fasta')

#####################################################################
# Load data, prepare counts
PP1_tab = fread('Data/PP1_table.tsv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')

PP2_tab = fread('Data/PP2_table.tsv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')

MTG_tab = fread('Data/MTG_table.tsv')
KEGG_tab = fread('Data/MTG_KEGG_counts.tsv')
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')

PP1_metadata$seq_number = vapply(PP1_metadata$Sample, function(x) sum(PP1_tab[,..x]), FUN.VALUE = numeric(1))
PP2_metadata$seq_number = vapply(PP2_metadata$Sample, function(x) sum(PP2_tab[,..x]), FUN.VALUE = numeric(1))
PP1_metadata$ASV_number = vapply(PP1_metadata$Sample, function(x) sum(PP1_tab[,..x] > 0), FUN.VALUE = numeric(1))
PP2_metadata$ASV_number = vapply(PP2_metadata$Sample, function(x) sum(PP2_tab[,..x] > 0), FUN.VALUE = numeric(1))
PP1_metadata$Dataset = 'PP1'
PP2_metadata$Dataset = 'PP2'
full_metadata = rbind(PP1_metadata, PP2_metadata)

# 1.1 Sample sequence number
ggplot(full_metadata) +
  geom_boxplot(aes(x = Ecosystem, y = seq_number, colour = Dataset)) +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  ylab('Sequence number') +
  scale_y_log10() + 
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('0_Dataset_exploration/SEQn.pdf', width = 6, height = 4)

# 1.2 Sample ASV number
ggplot(full_metadata) +
  geom_boxplot(aes(x = Ecosystem, y = ASV_number, colour = Dataset)) +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  ylab('Observed ASVs') +
  scale_y_log10() + 
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('0_Dataset_exploration/ASVn.pdf', width = 6, height = 4)

# 1.2 ASV number against sequencing depth
PP1_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP1',], formula = seq_number ~ a + b * log(ASV_number), start = list(a = 0, b = 0))
PP2_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP2',], formula = seq_number ~ a + b * log(ASV_number), start = list(a = 0, b = 0))
fit_data = data.frame(ASV_number=seq(from=2, to=10000, length.out=100))
fit_data$PP1_fit = predict(PP1_fit, list(ASV_number = fit_data$ASV_number), interval="prediction", type="l")
fit_data$PP2_fit = predict(PP2_fit, list(ASV_number = fit_data$ASV_number), interval="prediction", type="l")

ggplot(full_metadata) +
  geom_point(aes(x = ASV_number, y = seq_number, colour = Dataset), alpha=0.2) +
  geom_line(data=fit_data, aes(x=ASV_number, y=PP1_fit), color='#1F9DDB') +
  geom_line(data=fit_data, aes(x=ASV_number, y=PP2_fit), color='#DB0600') +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  xlab('Observed ASVs') + xlim(0,10000) +
  ylab('Sequencing depth') + ylim(2000, 185000) +
  scale_y_log10() + theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('0_Dataset_exploration/ASVn_SEQn.pdf', width = 7, height = 6)

#####################################################################
# 2 Number of ASVs, per category, taxonomic classification level, create Genus-level count table
library(dplyr)

PP1_tab_genus = PP1_tab[grepl('; g__',PP1_tab$Taxonomy),]
PP1_tab_genus = PP1_tab_genus[!grepl('; g__uncultured', PP1_tab_genus$Taxonomy),]
PP1_tab_genus = PP1_tab_genus[!grepl('; g__Unknown_Family', PP1_tab_genus$Taxonomy),]

PP2_tab_genus = PP2_tab[grepl('; g__',PP2_tab$Taxonomy),]
PP2_tab_genus = PP2_tab_genus[!grepl('; g__uncultured', PP2_tab_genus$Taxonomy),]
PP2_tab_genus = PP2_tab_genus[!grepl('; g__Unknown_Family', PP2_tab_genus$Taxonomy),]

MTG_tab_genus = MTG_tab[grepl('; g__',MTG_tab$Taxonomy),]
MTG_tab_genus = MTG_tab_genus[!grepl('; g__uncultured', MTG_tab_genus$Taxonomy),]
MTG_tab_genus = MTG_tab_genus[!grepl('incertae sedis', MTG_tab_genus$Taxonomy),]

#PP1_tab_genus$Genus = vapply(PP1_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
#PP2_tab_genus$Genus = vapply(PP2_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
#MTG_tab_genus$Genus = vapply(MTG_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
PP1_tab_genus$Genus = PP1_tab_genus$Taxonomy
PP2_tab_genus$Genus = PP2_tab_genus$Taxonomy
MTG_tab_genus$Genus = MTG_tab_genus$Taxonomy

PP1_tab_genus$ASV = NULL
PP2_tab_genus$ASV = NULL
MTG_tab_genus$mOTU = NULL

PP1_tab_genus$Taxonomy = NULL
PP2_tab_genus$Taxonomy = NULL
MTG_tab_genus$Taxonomy = NULL

PP1_tab_genus_out = PP1_tab_genus %>% group_by(Genus) %>% summarise_at(vars(PP1_metadata$Sample), sum)
write.csv(PP1_tab_genus_out, 'Data/PP1_genus.csv')
PP2_tab_genus_out = PP2_tab_genus %>% group_by(Genus) %>% summarise_at(vars(PP2_metadata$Sample), sum)
write.csv(PP2_tab_genus_out, 'Data/PP2_genus.csv')
MTG_tab_genus_out = MTG_tab_genus %>% group_by(Genus) %>% summarise_at(vars(MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]), sum)
# find smaller non zero value to multiply 1/value for pseudocounts
MTG_mat = MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]]
min_value = min(MTG_mat[MTG_mat > 0])
pseudocount_factor = 1/min_value
MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] = MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] * pseudocount_factor
MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] = apply(MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]], 2 , as.integer)
write.csv(MTG_tab_genus_out, 'Data/MTG_genus.csv')
# 3105 genera in PP1
# 3015 genera in PP2
# 1911 genera in MTG

#####################################################################
# 3. Number of samples per category
# Cryo vs others
PP1_Cryo = as.data.frame(table(PP1_metadata$Cryo))
PP2_Cryo = as.data.frame(table(PP2_metadata$Cryo))
MTG_Cryo = as.data.frame(table(MTG_metadata$Cryo[MTG_metadata$Sample %in% colnames(MTG_tab)]))
FUN_Cryo = as.data.frame(table(MTG_metadata$Cryo[MTG_metadata$Sample %in% colnames(KEGG_tab)]))
PP1_Cryo$Dataset = 'PP1'
PP2_Cryo$Dataset = 'PP2'
MTG_Cryo$Dataset = 'MTG'
FUN_Cryo$Dataset = 'FUN'
Cryo_data = rbind(PP1_Cryo, PP2_Cryo, MTG_Cryo, FUN_Cryo)
levels(Cryo_data$Var1) = c('Others', 'Cryosphere')

ggplot(Cryo_data) +
  geom_bar(aes(x = Dataset, y = Freq, fill = Dataset), stat = 'identity') +
  facet_wrap(~Var1, scales = 'free_y') +
  scale_fill_manual(values = c('#23A671','#5CAD0A', '#1A87BD', '#BA0600')) +
  ylab('Sample size') + xlab('') + scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"), legend.position = 'none')
ggsave('0_Dataset_exploration/FULL_sample_sizes.pdf', width = 5, height = 4)


#####################################################################
# 4. World maps
# Load the metadata again
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')

# load the functional data to subset only metaG with functional information
func_data = read.csv('Data/MTG_KEGG_counts.tsv', sep='\t')
mtg_samples = colnames(func_data)[colnames(func_data) != 'geneid']
MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% mtg_samples,]

# Merge metadata
MTG_metadata$Comments = NULL
PP1_metadata$Dataset = 'PP1'
PP2_metadata$Dataset = 'PP2'
# create mtg and functional metadata
MTG_metadata$Dataset = 'MTG'
MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_tab),]
full_metadata = rbind(PP1_metadata, PP2_metadata, MTG_metadata)
library(ggplot2)
library(maps)
library(ggrepel)
WorldData <- ggplot2::map_data('world') %>% fortify

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "#F2F0F7", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) + xlab('') + ylab('') +
  geom_point(data=mod_metadata, aes(x=Longitude, y=Latitude, colour = Dataset, shape=Cryosphere, size=`Sample n.`), alpha=0.9) + 
  scale_size_continuous(range = c(2,5), trans = 'log10') + 
  theme_bw() + theme(legend.title = element_text(size=9), 
                     legend.text=element_text(size=6), 
                     axis.title=element_text(size=8), 
                     legend.position="bottom",
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     legend.margin=margin(t = 0, unit='cm')) +
  scale_colour_manual(values = c("#00A087FF","#3C5488FF","#DC0000FF")) + 
  guides(colour = guide_legend(override.aes = list(alpha=1, size=3)), shape = guide_legend(override.aes = list(size=3))) 
ggsave('0_Dataset_exploration/Worldmap_all_datasets.pdf', width = 7, height = 4)

#####################################################################
# 4. Ecosystems sample sizes
PP1_eco = as.data.frame(table(PP1_metadata$Ecosystem[PP1_metadata$Cryosphere == 'Yes']))
PP2_eco = as.data.frame(table(PP2_metadata$Ecosystem[PP2_metadata$Cryosphere == 'Yes']))
MTG_eco = as.data.frame(table(MTG_metadata$Ecosystem[(MTG_metadata$Sample %in% colnames(MTG_tab)) & (MTG_metadata$Cryosphere == 'Yes')]))
FUN_eco = as.data.frame(table(MTG_metadata$Ecosystem[(MTG_metadata$Sample %in% colnames(KEGG_tab)) & (MTG_metadata$Cryosphere == 'Yes')]))
PP1_eco$Dataset = 'PP1'
PP2_eco$Dataset = 'PP2'
MTG_eco$Dataset = 'MTG'
FUN_eco$Dataset = 'FUN'
ECO_data = rbind(PP1_eco, PP2_eco, MTG_eco, FUN_eco)

ggplot(ECO_data) +
  geom_bar(aes(x = Dataset, y = Freq, fill = Dataset), stat = 'identity') +
  facet_wrap(~Var1, scales = 'free_y') +
  scale_fill_manual(values = c('#8AD8E6','#5CAD0A', '#1A87BD', '#BA0600')) +
  ylab('Sample size') + xlab('') +
  theme_bw() + scale_y_log10() +
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"), legend.position = 'none')
ggsave('0_Dataset_exploration/Cryo_sample_sizes.pdf', width = 5, height = 4)


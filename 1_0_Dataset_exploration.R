library(data.table)
library(ggplot2)
setwd('SET_THE_WORKING_DIRECTORY_HERE')
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
ggsave('1_ASV_analysis/1_0_Dataset_exploration/SEQn.pdf', width = 6, height = 4)

# 1.2 Sample ASV number
ggplot(full_metadata) +
  geom_boxplot(aes(x = Ecosystem, y = ASV_number, colour = Dataset)) +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  ylab('Observed ASVs') +
  scale_y_log10() + 
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('1_ASV_analysis/1_0_Dataset_exploration/ASVn.pdf', width = 6, height = 4)

# 1.2 ASV number against sequencing depth
PP1_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP1',], formula = seq_number ~ a + b * log(ASV_number), start = list(a = 0, b = 0))
PP2_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP2',], formula = seq_number ~ a + b * log(ASV_number), start = list(a = 0, b = 0))
fit_data = data.frame(ASV_number=seq(from=2, to=10000, length.out=100))
fit_data$PP1_fit = predict(PP1_fit, list(ASV_number = fit_data$ASV_number), interval="prediction", type="l")
fit_data$PP2_fit = predict(PP2_fit, list(ASV_number = fit_data$ASV_number), interval="prediction", type="l")

ggplot(full_metadata) +
  geom_point(aes(x = log(ASV_number), y = log(seq_number), colour = Dataset), alpha=0.2) +
  #geom_line(data=fit_data, aes(x=ASV_number, y=PP1_fit), color='#1F9DDB') +
  #geom_line(data=fit_data, aes(x=ASV_number, y=PP2_fit), color='#DB0600') +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  xlab('log(Observed ASVs)') + 
  ylab('log(Sequencing depth)') + geom_smooth(aes(x = log(ASV_number), y = log(seq_number), colour = Dataset)) +
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('1_ASV_analysis/1_0_Dataset_exploration/ASVn_SEQn.pdf', width = 7, height = 6)

#####################################################################
# 2 Number of ASVs, per category, taxonomic classification level, create Genus-level count table
library(dplyr)

PP1_tab_genus = PP1_tab[grepl('; g__',PP1_tab$Taxonomy),]
PP1_tab_genus = PP1_tab_genus[!grepl('; g__uncultured', PP1_tab_genus$Taxonomy),]
PP1_tab_genus = PP1_tab_genus[!grepl('; g__Unknown_Family', PP1_tab_genus$Taxonomy),]

PP2_tab_genus = PP2_tab[grepl('; g__',PP2_tab$Taxonomy),]
PP2_tab_genus = PP2_tab_genus[!grepl('; g__uncultured', PP2_tab_genus$Taxonomy),]
PP2_tab_genus = PP2_tab_genus[!grepl('; g__Unknown_Family', PP2_tab_genus$Taxonomy),]

#MTG_tab_genus = MTG_tab[grepl('; g__',MTG_tab$Taxonomy),]
#MTG_tab_genus = MTG_tab_genus[!grepl('; g__uncultured', MTG_tab_genus$Taxonomy),]
#MTG_tab_genus = MTG_tab_genus[!grepl('incertae sedis', MTG_tab_genus$Taxonomy),]

#PP1_tab_genus$Genus = vapply(PP1_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
#PP2_tab_genus$Genus = vapply(PP2_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
#MTG_tab_genus$Genus = vapply(MTG_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
PP1_tab_genus$Genus = PP1_tab_genus$Taxonomy
PP2_tab_genus$Genus = PP2_tab_genus$Taxonomy
#MTG_tab_genus$Genus = MTG_tab_genus$Taxonomy

PP1_tab_genus$ASV = NULL
PP2_tab_genus$ASV = NULL
#MTG_tab_genus$mOTU = NULL

PP1_tab_genus$Taxonomy = NULL
PP2_tab_genus$Taxonomy = NULL
#MTG_tab_genus$Taxonomy = NULL


PP1_tab_genus_out = PP1_tab_genus %>% group_by(Genus) %>% summarise_at(vars(PP1_metadata$Sample), sum)
write.csv(PP1_tab_genus_out, 'Data/PP1_genus.csv')
PP2_tab_genus_out = PP2_tab_genus %>% group_by(Genus) %>% summarise_at(vars(PP2_metadata$Sample), sum)
write.csv(PP2_tab_genus_out, 'Data/PP2_genus.csv')
#MTG_tab_genus_out = MTG_tab_genus %>% group_by(Genus) %>% summarise_at(vars(MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]), sum)
# find smaller non zero value to multiply 1/value for pseudocounts
#MTG_mat = MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]]
#min_value = min(MTG_mat[MTG_mat > 0])
#pseudocount_factor = 1/min_value
#MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] = MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] * pseudocount_factor
#MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]] = apply(MTG_tab_genus_out[,MTG_metadata$Sample[MTG_metadata$Sample %in% colnames(MTG_tab)]], 2 , as.integer)
#write.csv(MTG_tab_genus_out, 'Data/MTG_genus.csv')
# 3105 genera in PP1
# 3015 genera in PP2
# 1911 genera in MTG

#####################################################################
# 3. World maps
# Load the metadata again
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
MTG_metadata = read.csv('Metadata/MTG_metadata.txt', sep = '\t')

# Merge metadata
MTG_metadata$Comments = NULL
PP1_metadata$Dataset = 'PP1'
PP2_metadata$Dataset = 'PP2'
# create mtg and functional metadata
MTG_metadata$Dataset = 'MTG'
MTG_metadata$Habitat = ''
MTG_metadata$Date = ''
#MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_tab),]
full_metadata = rbind(PP1_metadata, PP2_metadata, MTG_metadata)

full_metadata$Latitude_round = round(full_metadata$Latitude)
full_metadata$Longitude_round = round(full_metadata$Longitude)
full_metadata$Sample_n = 1

library(dplyr)
count_metadata = full_metadata %>% group_by(Longitude_round, Latitude_round, Cryosphere, Dataset) %>% summarise(Sample_n = sum(Sample_n))
# Amplicon: others:3569   Cryo:768

library(ggplot2)
library(maps)
library(ggrepel)
WorldData <- ggplot2::map_data('world') %>% fortify

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id = region),
           fill = "grey", colour = "grey", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) + xlab('') + ylab('') +
  geom_point(data=count_metadata, aes(x=Longitude_round, y=Latitude_round, colour = Dataset, shape=Cryosphere, size=Sample_n), alpha=0.7) + 
  scale_size_continuous(range = c(2,6), trans = 'log10') + 
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
ggsave('1_ASV_analysis/1_0_Dataset_exploration/Worldmap_all_datasets.pdf', width = 7, height = 4)

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
ggsave('1_ASV_analysis/1_0_Dataset_exploration/Cryo_sample_sizes.pdf', width = 5, height = 4)


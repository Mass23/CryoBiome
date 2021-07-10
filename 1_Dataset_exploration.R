library(data.table)
library(ggplot2)
library(ggpubr)

PP1_tab = fread('Data/PP1_table.tsv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')

PP2_tab = fread('Data/PP2_table.tsv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')

MTG_tab = fread('Data/MTG_table.tsv')
KEGG_tab = fread('Data/MTG_KEGG_counts.tsv')
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')

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
#####################################################################
# 1.1 Sample sequence number
PP1_metadata$seq_number = vapply(PP1_metadata$Sample, function(x) sum(PP1_tab[,..x]), FUN.VALUE = numeric(1))
PP2_metadata$seq_number = vapply(PP2_metadata$Sample, function(x) sum(PP2_tab[,..x]), FUN.VALUE = numeric(1))
PP1_metadata$Dataset = 'PP1'
PP2_metadata$Dataset = 'PP2'
full_metadata = rbind(PP1_metadata, PP2_metadata)

ggplot(full_metadata) + 
  geom_boxplot(aes(x = Ecosystem, y = seq_number, colour = Dataset)) + 
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  ylab('Sequence number') +
  theme_bw()
ggsave('1_Dataset_exploration/SEQn.pdf', width = 6, height = 4)
  
# 1.2 Sample ASV number 
PP1_metadata$ASV_number = vapply(PP1_metadata$Sample, function(x) sum(PP1_tab[,..x] > 0), FUN.VALUE = numeric(1))
PP2_metadata$ASV_number = vapply(PP2_metadata$Sample, function(x) sum(PP2_tab[,..x] > 0), FUN.VALUE = numeric(1))
full_metadata = rbind(PP1_metadata, PP2_metadata)

ggplot(full_metadata) + 
  geom_boxplot(aes(x = Ecosystem, y = ASV_number, colour = Dataset)) + 
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  ylab('ASV number') +
  theme_bw()
ggsave('1_Dataset_exploration/ASVn.pdf', width = 6, height = 4)

# 1.2 ASV number against sequencing depth
PP1_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP1',], formula = ASV_number ~ a + b * log(seq_number), start = list(a = 0, b = 0))
PP2_fit <- nls(data=full_metadata[full_metadata$Dataset == 'PP2',], formula = ASV_number ~ a + b * log(seq_number), start = list(a = 0, b = 0))
fit_data = data.frame(seq_number=seq(from=0, to=2000000, length.out=1000))
fit_data$PP1_fit = predict(PP1_fit, list(seq_number = fit_data$seq_number), interval="prediction", type="l")
fit_data$PP2_fit = predict(PP2_fit, list(seq_number = fit_data$seq_number), interval="prediction", type="l")

ggplot(full_metadata) + 
  geom_point(aes(x = seq_number, y = ASV_number, colour = Dataset)) +
  scale_colour_manual(values=c('#1A87BD','#BA0600')) +
  xlab('Sequence number') +
  ylab('ASV number') +
  xlim(c(0,2000000)) + 
  ylim(c(0,15000)) +
  geom_line(data=fit_data, aes(x=seq_number, y=PP1_fit), color='#1F9DDB') +
  geom_line(data=fit_data, aes(x=seq_number, y=PP2_fit), color='#DB0600') +
  theme_bw()
ggsave('1_Dataset_exploration/ASVn_SEQn.pdf', width = 7, height = 6)

#####################################################################
# 2 Number of ASVs, per category, taxonomic classification level, create Genus-level count table
library(dplyr)

PP1_tab_genus = PP1_tab[grepl('; g__',PP1_tab$Taxonomy),]
PP1_tab_genus = PP1_tab_genus[!grepl('; g__uncultured', PP1_tab_genus$Taxonomy),]

PP2_tab_genus = PP2_tab[grepl('; g__',PP2_tab$Taxonomy),]
PP2_tab_genus = PP2_tab_genus[!grepl('; g__uncultured', PP2_tab_genus$Taxonomy),]

MTG_tab_genus = MTG_tab[grepl('; g__',MTG_tab$Taxonomy),]
MTG_tab_genus = MTG_tab_genus[!grepl('; g__uncultured', MTG_tab_genus$Taxonomy),]
MTG_tab_genus = MTG_tab_genus[!grepl('incertae sedis', MTG_tab_genus$Taxonomy),]

PP1_tab_genus$Genus = vapply(PP1_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
PP2_tab_genus$Genus = vapply(PP2_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))
MTG_tab_genus$Genus = vapply(MTG_tab_genus$Taxonomy, function(x) strsplit(x,'; g__')[[1]][2], FUN.VALUE = character(1))

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

#####################################################################
# 3. Number of samples per category
# EUCI vs non-EUCI
PP1_EUCI = as.data.frame(table(PP1_metadata$EUCI))
PP2_EUCI = as.data.frame(table(PP2_metadata$EUCI))
MTG_EUCI = as.data.frame(table(MTG_metadata$EUCI[MTG_metadata$Sample %in% colnames(MTG_tab)]))
FUN_EUCI = as.data.frame(table(MTG_metadata$EUCI[MTG_metadata$Sample %in% colnames(KEGG_tab)]))
PP1_EUCI$Dataset = 'PP1'
PP2_EUCI$Dataset = 'PP2'
MTG_EUCI$Dataset = 'MTG'
FUN_EUCI$Dataset = 'FUN'
EUCI_data = rbind(PP1_EUCI, PP2_EUCI, MTG_EUCI, FUN_EUCI)
levels(EUCI_data$Var1) = c('Others', 'EUCI')

ggplot(EUCI_data) + 
  geom_bar(aes(x = Dataset, y = Freq, fill = Dataset), stat = 'identity') + 
  facet_wrap(~Var1, scales = 'free_y') +
  scale_fill_manual(values = c('#23A671','#5CAD0A', '#1A87BD', '#BA0600')) +
  ylab('Sample size') + xlab('') + scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"), legend.position = 'none')
ggsave('1_Dataset_exploration/FULL_sample_sizes.pdf', width = 5, height = 4)


#####################################################################
# 4. World maps
# Load the metadata again
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')

# Merge metadata
MTG_metadata$Comments = NULL
PP1_metadata$Dataset = 'PP1'
PP2_metadata$Dataset = 'PP2'
# create mtg and functional metadata 
MTG_metadata$Dataset = 'MTG'
FUN_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(KEGG_tab),]
FUN_metadata$Dataset = 'FUN'
MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_tab),]
full_metadata = rbind(PP1_metadata, PP2_metadata, MTG_metadata, FUN_metadata)
full_metadata$colours = paste(full_metadata$Dataset, full_metadata$EUCI,sep='_')

library(ggplot2)
library(maps)
WorldData <- ggplot2::map_data('world') %>% fortify

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata), aes(x = Longitude, y = Latitude, colour = Dataset), size = 2) + 
  scale_colour_manual(values = c('#23A671','#5CAD0A','#1A87BD','#BA0600')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_all_datasets.pdf', width = 7, height = 5)

# Fig 1 map, shuffle rows for having both PP1 and PP2 displayed in front and not covering each other
full_metadata <- full_metadata[sample(nrow(full_metadata)), ]
ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[full_metadata$Dataset %in% c('PP1','PP2'),]), aes(x = Longitude, y = Latitude, alpha = EUCI, colour = colours, size = EUCI)) + 
  scale_size_manual(values=c(2,3)) +
  scale_colour_manual(values = c('#90ABD6','#0062A6','#E68294','#C20000')) + 
  scale_alpha_manual(values = c(0.5,1)) +
  theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_Amplicons.pdf', width = 7, height = 6)

#####################################################################
# 4. Ecosystem maps
PP1_eco = as.data.frame(table(PP1_metadata$Ecosystem[PP1_metadata$EUCI == 'Yes']))
PP2_eco = as.data.frame(table(PP2_metadata$Ecosystem[PP2_metadata$EUCI == 'Yes']))
MTG_eco = as.data.frame(table(MTG_metadata$Ecosystem[(MTG_metadata$Sample %in% colnames(MTG_tab)) & (MTG_metadata$EUCI == 'Yes')]))
FUN_eco = as.data.frame(table(MTG_metadata$Ecosystem[(MTG_metadata$Sample %in% colnames(KEGG_tab)) & (MTG_metadata$EUCI == 'Yes')]))
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
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"), legend.position = 'none')
ggsave('1_Dataset_exploration/EUCI_sample_sizes.pdf', width = 5, height = 4)


# maps
ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[(full_metadata$EUCI == 'Yes') & (full_metadata$Ecosystem == 'Snow'),]), aes(x = Longitude, y = Latitude, colour = Dataset), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c('#1A87BD')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_snow.pdf', width = 7, height = 5)

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[(full_metadata$EUCI == 'Yes') & (full_metadata$Ecosystem == 'Ice'),]), aes(x = Longitude, y = Latitude, colour = Dataset), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c('#8AD8E6','#5CAD0A','#1A87BD','#BA0600')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_ice.pdf', width = 7, height = 5)

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[(full_metadata$EUCI == 'Yes') & (full_metadata$Ecosystem == 'Terrestrial'),]), aes(x = Longitude, y = Latitude, colour = Dataset), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c('#8AD8E6','#5CAD0A','#1A87BD','#BA0600')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_terrestrial.pdf', width = 7, height = 5)

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[(full_metadata$EUCI == 'Yes') & (full_metadata$Ecosystem == 'Freshwater'),]), aes(x = Longitude, y = Latitude, colour = Dataset), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c('#8AD8E6','#5CAD0A','#1A87BD','#BA0600')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_freshwater.pdf', width = 7, height = 5)

ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(full_metadata[(full_metadata$EUCI == 'Yes') & (full_metadata$Ecosystem == 'Marine'),]), aes(x = Longitude, y = Latitude, colour = Dataset), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c('#8AD8E6','#5CAD0A','#1A87BD','#BA0600')) + theme(legend.position = 'none')
ggsave('1_Dataset_exploration/Worldmap_marine.pdf', width = 7, height = 5)

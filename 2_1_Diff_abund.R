library(ggplot2)
library(nlme)
library(compositions)
library(dplyr)

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

setwd('SET_THE_WORKING_DIRECTORY_HERE')
source('Scripts/ancom_v2.1.R')
############################################################################################################
# Differential abundance analysis
PP1_genus_tab = read.csv('Data/PP1_genus.csv')
PP2_genus_tab = read.csv('Data/PP2_genus.csv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP1_metadata$Dataset = 'PP1'
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
PP2_metadata$Dataset = 'PP2'

# Amplicon analysis
PP1_data = PP1_genus_tab[,PP1_metadata$Sample[PP1_metadata$Sample %in% colnames(PP1_genus_tab)]]
rownames(PP1_data) = PP1_genus_tab$Genus

PP2_data = PP2_genus_tab[,PP2_metadata$Sample[PP2_metadata$Sample %in% colnames(PP2_genus_tab)]]
rownames(PP2_data) = PP2_genus_tab$Genus

merged_table = merge(PP1_data,PP2_data,by='row.names',all=TRUE)
rownames(merged_table) = merged_table$Row.names
merged_table$Row.names = NULL
merged_table$X.x = NULL
merged_table$Genus.x = NULL
merged_table$X.y = NULL
merged_table$Genus.y = NULL
merged_table[is.na(merged_table)] = 0
merged_metadata = rbind(PP1_metadata, PP2_metadata)
length(unique(merged_metadata$Study))
length(unique(merged_metadata$Study[merged_metadata$Cryosphere == 'Yes']))
# 104 studies, 31 from cryospheric ecosystems
length(unique(merged_metadata$Sample))
length(unique(merged_metadata$Sample[merged_metadata$Cryosphere == 'Yes']))
# 4247 samples, 695 cryospheric (185 more cryo samples than ver. 1)

# COUNT NORMALISATION
merged_metadata$Cryosphere = as.character(merged_metadata$Cryosphere)
merged_metadata$Sample = as.character(merged_metadata$Sample)
merged_metadata$Dataset = as.character(merged_metadata$Dataset)
merged_metadata$Study = as.character(merged_metadata$Study)

# zero_cut = 0.99 means that only genera present in at least 21 samples are taken into account
feature_table = merged_table; sample_var = "Sample"; group_var = 'Cryosphere'
out_cut = 0.05; zero_cut = 0.995; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, merged_metadata, sample_var, group_var, out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
main_var = 'Cryosphere'; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | Dataset"
control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha, adj_formula, rand_formula, control=control)

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data W=1277
colnames(res$fig$data)[colnames(res$fig$data) == 'x'] = 'CLR_mean_diff'
colnames(res$fig$data)[colnames(res$fig$data) == 'y'] = 'W'
dat_ann = data.frame(x = min(res$fig$data$CLR_mean_diff), y = cut_off["detected_0.7"], label = "W[0.7]")

res$fig$data$DA = 'Others'
res$fig$data$DA[(res$fig$data$W>dat_ann$y)&(res$fig$data$CLR_mean_diff>0)] = 'Overrepresented'

ggplot(res$fig$data) + geom_point(aes(x=CLR_mean_diff,y=W,color=DA)) + scale_color_manual(values=c('dimgrey','#1A7A7F')) + theme_linedraw() +
  geom_hline(yintercept = dat_ann$y, linetype = "dashed") + ylab('W statistic') + xlab('CLR mean difference') + theme(legend.position = 'none') +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  geom_vline(xintercept = 0,size=0.2)
ggsave('2_Genus_analysis/2_1_Diff_abund/2_1_Amplicon_ancom_res.pdf',width = 4,height = 4)

over_genera = as.character(res$fig$data$taxa_id[(res$fig$data$CLR_mean_diff>0) & (res$fig$data$W > dat_ann$y)])
over_data = res$fig$data[(res$fig$data$CLR_mean_diff>0)&(res$fig$data$W>dat_ann$y),]

write.csv(res$fig$data, file='Data/Amplicon_ancom_res.csv')

# Most overrepresented genera barplot w=1251 for 0.7
over_data = read.csv('Data/Amplicon_ancom_res.csv')
over_data$Genus = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; g__')[[1]][2], FUN.VALUE = character(1))
over_data$Family = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][5], FUN.VALUE = character(1))
over_data$Order = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][4], FUN.VALUE = character(1))
over_data$Class = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][3], FUN.VALUE = character(1))
over_data$Phylum = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][2], FUN.VALUE = character(1))
over_data$Phylum = gsub('p__','',over_data$Phylum)

table(over_data$Phylum[over_data$DA == 'Overrepresented'])

# phylum: c('Acidobacteriota','Actinobacteriota','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Others,'Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota')
# colors: c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')
most_over = over_data[(over_data$CLR_mean_diff > 0.35) & (over_data$W >= 1269),]
ggplot(most_over) + 
  geom_bar(aes(x=CLR_mean_diff,y=reorder(as.factor(Genus), CLR_mean_diff),fill=Phylum), stat='identity') + ylab('') + xlab('CLR mean difference') +
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'), axis.text.y = element_text(face = "italic")) + 
  scale_fill_manual(values=c('#F39B7FFF','#DC0000FF','#00A087FF','grey','dimgrey','#3C5488FF','#8491B4FF'))
ggsave('2_Genus_analysis/2_1_Diff_abund/2_1_Ancom_most_over.pdf', width=6.5, height = 6)


# Create Taxonomic tree at the family level with the edge size as the number of overrepresented genera
library(metacoder)
tree_data = parse_tax_data(data.frame(Taxonomy= over_data$taxa_id[(over_data$CLR_mean_diff > 0) & (over_data$W >= 1269)]),
                           class_cols = "Taxonomy", # the column that contains taxonomic information
                           class_sep = "; ", # The character used to separate taxa in the classification
                           class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c(tax_rank = "info", # A key describing each regex capture group
                                         tax_name = "taxon_name"))
tree_data %>% filter_taxa(tree_data$data$class_data$taxon_id[tree_data$data$class_data$tax_rank == 'c'], supertaxa = TRUE) %>%
  filter_taxa(n_obs > 1, supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_color_range = c('grey','#B09C85FF','#00A087FF'),
            node_size=n_obs,
            node_color=n_obs,
            node_color_trans='log10',
            node_size_trans='log10',
            node_color_axis_label = "Cryo. genera #",
            margin_size = c(0.1,0.1,0.1,0.1),
            initial_layout = "re",
            overlap_avoidance = 10,
            output_file = "2_Genus_analysis/2_1_Diff_abund/2_1_Ancom_over_genera.pdf")

library(data.table)
library(ggplot2)
library(metacoder)
library(tibble)
library(tidyr)
library(foreach)
library(doMC)
library(plyr)
library(ggpubr)

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
setwd('/Users/admin/Documents/Academia/PhD/Chapter I/')

###########################################################################################################################
# 1. Data loading, we load only the EUCI samples data
# PP1
PP1_tab  = fread('Data/PP1_table.tsv')
PP1_metadata = fread('Metadata/PP1_metadata.tsv')
PP1_EUCI_cols = c(PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes'], 'Taxonomy', 'ASV')
PP1_EUCI_samples = PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes']
PP1_EUCI_tab = PP1_tab[,..PP1_EUCI_cols]
PP1_EUCI_tab = PP1_EUCI_tab[rowSums(PP1_EUCI_tab[,..PP1_EUCI_samples]) > 0,]
PP1_data = parse_tax_data(PP1_EUCI_tab,
                          class_cols = "Taxonomy", # the column that contains taxonomic information
                          class_sep = "; ", # The character used to separate taxa in the classification
                          class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                          class_key = c(tax_rank = "info", # A key describing each regex capture group
                                        tax_name = "taxon_name"))
PP1_data$data$tax_data <- calc_obs_props(PP1_data, "tax_data")
PP1_data$data$tax_abund <- calc_taxon_abund(PP1_data, "tax_data")

# PP2
PP2_tab  = fread('Data/PP2_table.tsv')
PP2_metadata = fread('Metadata/PP2_metadata.tsv')

PP2_metadata = PP2_metadata[!PP2_metadata$Study == 'PRJNA518106',]
PP2_EUCI_cols = c(PP2_metadata$Sample[(PP2_metadata$EUCI == 'Yes')], 'Taxonomy', 'ASV')
PP2_EUCI_samples = PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes']
PP2_EUCI_tab = PP2_tab[,..PP2_EUCI_cols]
PP2_EUCI_tab = PP2_EUCI_tab[rowSums(PP2_EUCI_tab[,..PP2_EUCI_samples]) > 0,]
PP2_data = parse_tax_data(PP2_EUCI_tab,
                          class_cols = "Taxonomy", # the column that contains taxonomic information
                          class_sep = "; ", # The character used to separate taxa in the classification
                          class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                          class_key = c(tax_rank = "info", # A key describing each regex capture group
                                        tax_name = "taxon_name"))
PP2_data$data$tax_data <- calc_obs_props(PP2_data, "tax_data")
PP2_data$data$tax_abund <- calc_taxon_abund(PP2_data, "tax_data")

# Load deseq results
PP1_res = read.csv('3_Specific_analysis/PP1_deseq_res.csv')
PP2_res = read.csv('3_Specific_analysis/PP2_deseq_res.csv')
genus_list = union(PP1_res$X,PP2_res$X)

###########################################################################################################################
# 2. Prevalence analysis
PP1_data$data$tax_occ <- calc_n_samples(PP1_data, "tax_abund", more_than = 1/2000)
PP2_data$data$tax_occ <- calc_n_samples(PP2_data, "tax_abund", more_than = 1/2000)
PP1_data$data$tax_occ$n_samples = PP1_data$data$tax_occ$n_samples / 238
PP2_data$data$tax_occ$n_samples = PP2_data$data$tax_occ$n_samples / 300

PP1_res$taxon_id = vapply(PP1_res$X, function(x) {ifelse(x %in% PP1_res$X, PP1_data$data$class_data$taxon_id[(PP1_data$data$class_data$tax_name == x) & (PP1_data$data$class_data$tax_rank == 'g')], NA)}, FUN.VALUE = character(1))
PP1_res$Dataset = 'PP1'
PP1_res$Prevalence = vapply(PP1_res$taxon_id, function(x) {ifelse(x %in% PP1_res$taxon_id, PP1_data$data$tax_occ$n_samples[PP1_data$data$tax_occ$taxon_id == x], NA)}, FUN.VALUE = numeric(1))

PP2_res$taxon_id = vapply(PP2_res$X, function(x) {ifelse(x %in% PP2_res$X, PP2_data$data$class_data$taxon_id[(PP2_data$data$class_data$tax_name == x) & (PP2_data$data$class_data$tax_rank == 'g')], NA)}, FUN.VALUE = character(1))
PP2_res$Dataset = 'PP2'
PP2_res$Prevalence = vapply(PP2_res$taxon_id, function(x) {ifelse(x %in% PP2_res$taxon_id, PP2_data$data$tax_occ$n_samples[PP2_data$data$tax_occ$taxon_id == x], NA)}, FUN.VALUE = numeric(1))

preval_deseq = rbind(PP1_res, PP2_res)

preval_deseq$Group = 'Others'
for (genus in unique(preval_deseq$X)){
  # Relaxed evidenc genera are: significant (p < 0.05) and differentially abundant (log2Fc > 1) in at least one of the two datasets
  if (sum((preval_deseq$padj[(preval_deseq$X == genus) & (preval_deseq$Dataset == 'PP1')] < 0.05)  & (preval_deseq$log2FoldChange[(preval_deseq$X == genus) & (preval_deseq$Dataset == 'PP1')] > 1)) > 0 |
      sum((preval_deseq$padj[(preval_deseq$X == genus) & (preval_deseq$Dataset == 'PP2')] < 0.05)  & (preval_deseq$log2FoldChange[(preval_deseq$X == genus) & (preval_deseq$Dataset == 'PP2')] > 1)) > 0 ){
    # Strict genera are significant and differentially abundant in both datasets, with a min. prevalence of 5%
    if ((sum(preval_deseq$padj[preval_deseq$X == genus] < 0.05) == 2) &
        (sum(preval_deseq$log2FoldChange[preval_deseq$X == genus] > 1) == 2) &
        (sum(preval_deseq$Prevalence[preval_deseq$X == genus] > 0.05) == 2)) {preval_deseq$Group[preval_deseq$X == genus] = 'Core'}
    else {preval_deseq$Group[preval_deseq$X == genus] = 'Ancillary'}}}

PP1_gen = unique(PP1_data$data$class_data$tax_name[PP1_data$data$class_data$tax_rank == 'g'])
PP2_gen = unique(PP2_data$data$class_data$tax_name[PP2_data$data$class_data$tax_rank == 'g'])
length(union(PP1_gen, PP2_gen))
length(unique(preval_deseq$X[preval_deseq$Group == 'Ancillary']))
length(unique(preval_deseq$X[preval_deseq$Group == 'Core']))
# 2227 genera, 571 relaxed, 71 strict

preval_deseq$Group  <- factor(preval_deseq$Group , levels = c("Core","Ancillary", "Others"))
ggplot(preval_deseq) +
  geom_boxplot(aes(x = Group, colour = Dataset, y = Prevalence)) +
  scale_color_manual(values = c('#1A87BD', '#BA0600')) + scale_y_log10() +
  xlab('') + guides(colour=FALSE) + theme_bw() + ylab('Prevalence')
ggsave('4_Taxonomic_analysis/Genus_prevalence.pdf', width = 3, height = 3)

######## STRICT #########
mean(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')])
median(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')])
# Mean prevalence, PP1, strict = 0.198485
# Median prevalence, PP1, strict = 0.1722689
mean(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')])
median(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')])
# Mean prevalence, PP2, strict = 0.1974648
# Median prevalence, PP2, strict = 0.1566667
wilcox.test(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')], preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')])
wilcox.test(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')], preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')])
# W = 30294, p-value < 2.2e-16
# W = 30102, p-value < 2.2e-16

######## RELAXED #########
mean(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
median(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
# Mean prevalence, PP1, relaxed = 0.08226153
# Median prevalence, PP1, relaxed = 0.04621849
mean(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')]))
median(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
# Mean prevalence, PP2, relaxed = 0.07442029
# Median prevalence, PP2, relaxed = 0.04621849
wilcox.test(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')], preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')])
wilcox.test(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')], preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')])
# W = 422655, p-value < 2.2e-16
# W = 484134, p-value < 2.2e-16

######## OTHERS #########
mean(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')]))
median(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')]))
# Mean prevalence, PP1, others = 0.01640796
# Median prevalence, PP1, others = 0.004201681
mean(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')]))
median(na.omit(preval_deseq$Prevalence[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')]))
# Mean prevalence, PP2, others = 0.02192139
# Median prevalence, PP2, others = 0.003333333

relaxed_genera = unique(preval_deseq$X[preval_deseq$Group == 'Ancillary'])
write.csv(relaxed_genera, '4_Taxonomic_analysis/relaxed_genera.csv', row.names = F)
strict_genera = unique(preval_deseq$X[preval_deseq$Group == 'Core'])
write.csv(strict_genera, '4_Taxonomic_analysis/strict_genera.csv', row.names = F)


###########################################################################################################################
# 2.1 Abundance analysis: abundance of the genera in EUCI
library(dplyr)
PP1_EUCI_samples = PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes']
PP2_EUCI_samples = PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes']
preval_deseq$mean_abund = vapply(1:nrow(preval_deseq), function(x){ifelse(preval_deseq$Dataset[x] == 'PP1', 
                                                                          sum(PP1_data$data$tax_abund[PP1_data$data$tax_abund$taxon_id == preval_deseq$taxon_id[x],] %>% summarise(across(PP1_EUCI_samples, mean))) / 238, 
                                                                          sum(PP2_data$data$tax_abund[PP2_data$data$tax_abund$taxon_id == preval_deseq$taxon_id[x],] %>% summarise(across(PP2_EUCI_samples, mean))) / 300)}, FUN.VALUE = numeric(1))

######## STRICT #########
mean(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')])
median(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')])
# Mean abundance, PP1, strict = 0.003885491
# Median abundance, PP1, strict = 0.001645029
mean(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')])
median(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')])
# Mean abundance, PP2, strict = 0.003190746
# Median abundance, PP2, strict = 0.001349517
wilcox.test(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Core')], preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')])
wilcox.test(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Core')], preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')])
# W = 31594, p-value < 2.2e-16
# W = 30837, p-value < 2.2e-16

######## RELAXED #########
mean(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
median(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
# Mean abundance, PP1, relaxed = 0.0008132448
# Median abundance, PP1, relaxed = 0.0001552926
mean(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')]))
median(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')]))
# Mean abundance, PP2, relaxed = 0.0007386248
# Median abundance, PP2, relaxed = 0.0001552926
wilcox.test(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Ancillary')], preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')])
wilcox.test(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Ancillary')], preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')])
# W = 428880, p-value < 2.2e-16
# W = 488646, p-value < 2.2e-16

######## OTHERS #########
mean(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')]))
median(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == 'Others')]))
# Mean abundance, PP1, others = 0.0001099757
# Median abundance, PP1, others = 7.442073e-06
mean(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')]))
median(na.omit(preval_deseq$mean_abund[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == 'Others')]))
# Mean abundance, PP2, others = 0.0001367853
# Median abundance, PP2, others = 8.431297e-06


preval_deseq$Group  <- factor(preval_deseq$Group , levels = c("Core","Ancillary", "Others"))
ggplot(preval_deseq) +
  geom_boxplot(aes(y = mean_abund, colour = Dataset, x = Group)) +
  scale_color_manual(values = c('#1A87BD', '#BA0600')) + scale_y_log10() +
  xlab('') + ylab('Relative abundance') + guides(colour=FALSE) + 
  theme_bw()
ggsave('4_Taxonomic_analysis/Genus_abundance.pdf', width = 3, height = 3)


# 2.2 Abundance analysis: Ecosystems comparison in genera contribution to communities
PP1_strict_ids = unique(PP1_data$data$class_data$taxon_id[(PP1_data$data$class_data$tax_name %in% strict_genera) & (PP1_data$data$class_data$tax_rank == 'g')])
PP1_relaxed_ids = unique(PP1_data$data$class_data$taxon_id[(PP1_data$data$class_data$tax_name %in% relaxed_genera) & (PP1_data$data$class_data$tax_rank == 'g')])
PP2_strict_ids = unique(PP2_data$data$class_data$taxon_id[(PP2_data$data$class_data$tax_name %in% strict_genera) & (PP2_data$data$class_data$tax_rank == 'g')])
PP2_relaxed_ids = unique(PP2_data$data$class_data$taxon_id[(PP2_data$data$class_data$tax_name %in% relaxed_genera) & (PP2_data$data$class_data$tax_rank == 'g')])


# Mean rel. abundance PP1
PP1_prop_assigned_genus = colSums(PP1_data$data$tax_abund[PP1_data$data$tax_abund$taxon_id %in% PP1_data$data$class_data$taxon_id[PP1_data$data$class_data$tax_rank == 'g'], PP1_metadata$Sample[PP1_metadata$EUCI=='Yes']])
PP1_prop_strict_genus = colSums(PP1_data$data$tax_abund[PP1_data$data$tax_abund$taxon_id %in% PP1_strict_ids, PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes']]) 
PP1_prop_relaxed_genus = colSums(PP1_data$data$tax_abund[PP1_data$data$tax_abund$taxon_id %in% PP1_relaxed_ids, PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes']]) 
# Mean rel. abundance PP1
PP2_prop_assigned_genus = colSums(PP2_data$data$tax_abund[PP2_data$data$tax_abund$taxon_id %in% PP2_data$data$class_data$taxon_id[PP2_data$data$class_data$tax_rank == 'g'], PP2_metadata$Sample[PP2_metadata$EUCI=='Yes']])
PP2_prop_strict_genus = colSums(PP2_data$data$tax_abund[PP2_data$data$tax_abund$taxon_id %in% PP2_strict_ids, PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes']]) 
PP2_prop_relaxed_genus = colSums(PP2_data$data$tax_abund[PP2_data$data$tax_abund$taxon_id %in% PP2_relaxed_ids, PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes']]) 

mean(PP1_prop_strict_genus / PP1_prop_assigned_genus)
# Mean rel. abund, strict = 29.31191 %
mean(PP2_prop_strict_genus / PP2_prop_assigned_genus)
# Mean rel. abund, strict = 25.35672 %

mean(PP1_prop_relaxed_genus / PP1_prop_assigned_genus)
# Mean rel. abund, strict = 44.82479 %
mean(PP2_prop_relaxed_genus / PP2_prop_assigned_genus)
# Mean rel. abund, strict = 41.53972 %

abund_deseq = data.frame(Dataset = c(), Ecosystem = c(), Group = c(), mean_abund = c(), se_low = c(), se_high = c())

ecosystems = c('Terrestrial', 'Snow', 'Ice', 'Marine', 'Freshwater')
groups = c('Others', 'Core', 'Ancillary')

for (ecosystem in ecosystems){
  for (group in groups){
    PP1_tax_ids = preval_deseq$taxon_id[(preval_deseq$Dataset == 'PP1') & (preval_deseq$Group == group)]
    PP2_tax_ids = preval_deseq$taxon_id[(preval_deseq$Dataset == 'PP2') & (preval_deseq$Group == group)]
    
    pp1_data = colSums(PP1_data$data$tax_abund[PP1_data$data$tax_abund$taxon_id %in% PP1_tax_ids, PP1_metadata$Sample[(PP1_metadata$Ecosystem == ecosystem) & (PP1_metadata$EUCI == 'Yes')]]) / PP1_prop_assigned_genus
    pp2_data = colSums(PP2_data$data$tax_abund[PP2_data$data$tax_abund$taxon_id %in% PP2_tax_ids, PP2_metadata$Sample[(PP2_metadata$Ecosystem == ecosystem) & (PP2_metadata$EUCI == 'Yes')]]) / PP2_prop_assigned_genus
    
    pp1_sd = sd(pp1_data)
    pp2_sd = sd(pp2_data)
    
    abund_deseq = rbind(abund_deseq, data.frame(Dataset = 'PP1', Ecosystem = ecosystem, Group = group, mean_abund = mean(pp1_data), sd_low = mean(pp1_data) - pp1_sd, sd_high = mean(pp1_data) + pp1_sd))
    abund_deseq = rbind(abund_deseq, data.frame(Dataset = 'PP2', Ecosystem = ecosystem, Group = group, mean_abund = mean(pp2_data), sd_low = mean(pp2_data) - pp2_sd, sd_high = mean(pp2_data) + pp2_sd))
  }
}

to_delete = which(abund_deseq$Dataset == 'PP2' & abund_deseq$Ecosystem == 'Snow')
abund_deseq = abund_deseq[-to_delete,]

abund_deseq$Group  <- factor(abund_deseq$Group , levels = c("Others", "Ancillary", "Core"))
ggplot(abund_deseq) + 
  geom_bar(aes(x=Dataset, y=mean_abund, fill = Dataset, alpha = Group), size = 0.9, stat = 'identity', position = "fill") +
  coord_flip() + #scale_y_reverse() + scale_y_continuous(breaks=c(1,0.75,0.5,0.25,0), labels=c(0,0.25,0.5,0.75,1), limits=c(0,1)) +
  geom_hline(yintercept = 0, lty=2) +
  facet_grid(Ecosystem ~ ., scales = "free",  space = 'free') +
  theme_bw() +
  xlab('') + ylab('Community contribution') +
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_fill_manual(values = c('#1A87BD', '#BA0600')) +
  scale_alpha_manual(values = c(0, 0.5, 1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))
ggsave('4_Taxonomic_analysis/Barplot_relative_abundance.pdf', width = 4, height = 6)


###########################################################################################################################
# 3. Heat trees of the strict and relaxed genera
relaxed_genera = read.table('4_Taxonomic_analysis/relaxed_genera.csv', sep = ',', header = TRUE)
strict_genera = read.table('4_Taxonomic_analysis/strict_genera.csv', sep = ',', header = TRUE)

PP1_strict_ids = unique(PP1_data$data$class_data$taxon_id[(PP1_data$data$class_data$tax_name %in% strict_genera$x) & (PP1_data$data$class_data$tax_rank == 'g')])
PP1_relaxed_ids = unique(PP1_data$data$class_data$taxon_id[(PP1_data$data$class_data$tax_name %in% relaxed_genera$x) & (PP1_data$data$class_data$tax_rank == 'g')])
PP2_strict_ids = unique(PP2_data$data$class_data$taxon_id[(PP2_data$data$class_data$tax_name %in% strict_genera$x) & (PP2_data$data$class_data$tax_rank == 'g')])
PP2_relaxed_ids = unique(PP2_data$data$class_data$taxon_id[(PP2_data$data$class_data$tax_name %in% relaxed_genera$x) & (PP2_data$data$class_data$tax_rank == 'g')])

# verify that it worked (should be 71 for both)
length(unique(PP1_data$data$class_data$tax_name[PP1_data$data$class_data$taxon_id %in% PP1_strict_ids]))
length(unique(PP2_data$data$class_data$tax_name[PP2_data$data$class_data$taxon_id %in% PP2_strict_ids]))

# Taxonomic tree of the strict IDs
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_subtaxa,
  node_color = n_supertaxa,
  overlap_avoidance = 1.5,
  aspect_ratio = 1.1,
  node_size_axis_label = "Subtaxa N",
  node_color_axis_label = "Supertaxa N",
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Strict_tree.pdf')

# Prevalence for both datasets
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = n_samples, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/PP1_strict_tree.pdf')

PP1_data %>% filter_taxa(na.omit(PP1_relaxed_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = n_samples, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/PP1_relaxed_tree.pdf')

PP2_data %>% filter_taxa(na.omit(PP2_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = n_samples, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/PP2_strict_tree.pdf')

PP2_data %>% filter_taxa(na.omit(PP2_relaxed_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = n_samples, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/PP2_relaxed_tree.pdf')

###########################################################################################################################
# 3. Heat trees of the different ecosystems
PP1_data$data$tax_occ_eco <- calc_n_samples(PP1_data, "tax_abund", more_than = 1/2000, cols = PP1_metadata$Sample[PP1_metadata$EUCI == 'Yes'], groups = PP1_metadata$Ecosystem[PP1_metadata$EUCI == 'Yes'])
PP2_data$data$tax_occ_eco <- calc_n_samples(PP2_data, "tax_abund", more_than = 1/2000, cols = PP2_metadata$Sample[PP2_metadata$EUCI == 'Yes'], groups = PP2_metadata$Ecosystem[PP2_metadata$EUCI == 'Yes'])

######## SNOW #########
PP1_data %>% filter_taxa(PP1_strict_ids, supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Snow, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Snow/PP1_snow_tree.pdf')

######## ICE #########
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Ice, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Ice/PP1_ice_tree.pdf')
PP2_data %>% filter_taxa(na.omit(PP2_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Ice, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Ice/PP2_ice_treet.pdf')

######## TERRESTRIAL #########
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Terrestrial, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Terrestrial/PP1_terr_tree.pdf')
PP2_data %>% filter_taxa(na.omit(PP2_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Terrestrial, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Terrestrial/PP2_terr_tree.pdf')

######## MARINE #########
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Marine, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Marine/PP1_mari_tree.pdf')
PP2_data %>% filter_taxa(na.omit(PP2_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Marine, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Marine/PP2_mari_tree.pdf')

######## FRESHWATER #########
PP1_data %>% filter_taxa(na.omit(PP1_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Freshwater, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#1A87BD"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Freshwater/PP1_frwa_tree.pdf')
PP2_data %>% filter_taxa(na.omit(PP2_strict_ids), supertaxa = TRUE) %>% heat_tree(
  node_label = taxon_names,
  node_size = n_obs,
  node_color = Freshwater, 
  node_size_axis_label = "ASV count",
  node_color_axis_label = "Prevalence",
  node_color_range = c("grey", "#BA0600"),
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  output_file = '4_Taxonomic_analysis/Freshwater/PP2_frwa_tree.pdf')


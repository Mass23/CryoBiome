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
# MTG: #6BC90C
# colors dark:
# PP1: #1A87BD
# PP2: #BA0600
# MTG: #5CAD0A
###########################################################################################################################
# 1. Data loading
# PP1
MTG_tab  = fread('Data/MTG_table.tsv')
MTG_metadata = fread('Metadata/MTG_metadata.tsv')
MTG_data = parse_tax_data(MTG_tab,
                          class_cols = "Taxonomy", # the column that contains taxonomic information
                          class_sep = "; ", # The character used to separate taxa in the classification
                          class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                          class_key = c(tax_rank = "info", # A key describing each regex capture group
                                        tax_name = "taxon_name"))
MTG_data$data$tax_data <- calc_obs_props(MTG_data, "tax_data")
MTG_data$data$tax_abund <- calc_taxon_abund(MTG_data, "tax_data")

strict_genera = read.csv('4_Taxonomic_analysis/strict_genera.csv')
relaxed_genera = read.csv('4_Taxonomic_analysis/relaxed_genera.csv')

MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')
MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_tab),]
MTG_EUCI_samples = MTG_metadata$Sample[MTG_metadata$EUCI == 'Yes']
MTG_COMP_samples = MTG_metadata$Sample[MTG_metadata$EUCI == 'No']

###########################################################################################################################
# 2. Prevalence
MTG_data$data$tax_occ <- calc_n_samples(MTG_data, "tax_abund", more_than = 1/2000, cols = MTG_EUCI_samples)
MTG_data$data$tax_occ$n_samples = MTG_data$data$tax_occ$n_samples / length(MTG_EUCI_samples)

MTG_prevalence_df = MTG_data$data$tax_occ[MTG_data$data$tax_occ$taxon_id %in% unique(MTG_data$data$class_data$taxon_id[MTG_data$data$class_data$tax_rank == 'g']),]
MTG_prevalence_df$taxon_name = vapply(MTG_prevalence_df$taxon_id, function(x) unique(MTG_data$data$class_data$tax_name[MTG_data$data$class_data$taxon_id == x]), FUN.VALUE = character(1))

MTG_prevalence_df$Group = 'Others'
MTG_prevalence_df$Group[MTG_prevalence_df$taxon_name %in% relaxed_genera$NCBI] = 'Ancillary'
MTG_prevalence_df$Group[MTG_prevalence_df$taxon_name %in% strict_genera$NCBI] = 'Core'

ggplot(MTG_prevalence_df) +
  geom_boxplot(aes(x = Group, y = n_samples), colour = '#5CAD0A') +
  xlab('') + ylab('Prevalence') + scale_y_log10() +
  theme_bw()
ggsave('6_MTG_taxonomy/MTG_prevalence.pdf', width = 3, height = 4)

###########################################################################################################################
# 3. Abundance
library(dplyr)
MTG_prevalence_df$mean_abund = vapply(MTG_prevalence_df$taxon_id, function(x)  sum(MTG_data$data$tax_abund[MTG_data$data$tax_abund$taxon_id == x,] %>% summarise(across(MTG_EUCI_samples, mean))) / length(MTG_EUCI_samples), FUN.VALUE = numeric(1))

ggplot(MTG_prevalence_df) +
  geom_boxplot(aes(y = mean_abund, x = Group), colour = '#5CAD0A') +
  xlab('') + ylab('Mean abundance in EUCI samples') + scale_y_log10() +
  theme_bw()
ggsave('6_MTG_taxonomy/MTG_abundance.pdf', width = 3, height = 4)


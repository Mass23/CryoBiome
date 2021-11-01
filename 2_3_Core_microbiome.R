library(ggplot2)
library(plyr)
library(UpSetR)
library(doMC)
library(foreach)
library(dplyr)
library(ComplexUpset)

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

setwd('/Users/mabourqu/Documents/PhD/C1/')
############################################################################################################
# Data laoding
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP1_metadata$Dataset = 'PP1'
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
PP2_metadata$Dataset = 'PP2'
full_metadata = rbind(PP1_metadata,PP2_metadata)
levels(full_metadata$Ecosystem)[levels(full_metadata$Ecosystem)=='Snow'] <- 'Snow/Ice'
levels(full_metadata$Ecosystem)[levels(full_metadata$Ecosystem)=='Ice'] <- 'Snow/Ice'

PP1_genus_tab = read.csv('Data/PP1_genus.csv')
PP2_genus_tab = read.csv('Data/PP2_genus.csv')
PP1_data = PP1_genus_tab[,PP1_metadata$Sample[(PP1_metadata$Sample %in% colnames(PP1_genus_tab)) & (PP1_metadata$Cryosphere == 'Yes')]]
rownames(PP1_data) = PP1_genus_tab$Genus
PP1_data$X = NULL
PP1_data$Genus = NULL
PP2_data = PP2_genus_tab[,PP2_metadata$Sample[(PP2_metadata$Sample %in% colnames(PP2_genus_tab)) & (PP2_metadata$Cryosphere == 'Yes')]]
rownames(PP2_data) = PP2_genus_tab$Genus
PP2_data$X = NULL
PP2_data$Genus = NULL
merged_table = merge(PP1_data,PP2_data,by='row.names',all=TRUE)
rownames(merged_table) = merged_table$Row.names
merged_table$Row.names = NULL
merged_table[is.na(merged_table)] = 0
merged_table = merged_table[rowSums(merged_table) > 0,]
merged_table = sweep(merged_table,2,colSums(merged_table),'/')

ancom_res = read.csv('Data/Amplicon_ancom_res.csv')

########################################################################################
# NMDS analysis
library(vegan)
library(tidyr)
library(pairwiseAdonis)
library(ggforce)

# testing k values
test_k <- function(k1, k2, data){
  k = c()
  stress = c()
  for (i in k1:k2){
    k = c(k, i)
    nmds = metaMDS(data, distance = "bray", k=i, trymax = 20)
    stress = c(stress, nmds$stress)}
  return(list(k=k,stress=stress))}
test_k = test_k(1,10,log(t(merged_table)+1))

ggplot(as.data.frame(test_k)) + geom_line(aes(x=k,y=stress)) + theme_linedraw()
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_k_comparison.pdf')

nmds = metaMDS(log(t(merged_table)+1), distance = "bray", k=2, trymax = 10000)
nmds_data = as.data.frame(scores(nmds))
nmds_data$Sample = colnames(merged_table)
nmds_data$Ecosystem = vapply(nmds_data$Sample, function(x){as.character(full_metadata$Ecosystem[full_metadata$Sample == x])}, FUN.VALUE = character(1))
nmds_data$Dataset = vapply(nmds_data$Sample, function(x){as.character(full_metadata$Dataset[full_metadata$Sample == x])}, FUN.VALUE = character(1))
nmds$stress # 0.208002

ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, shape = Ecosystem)) + geom_point(aes(colour = Dataset)) + scale_color_manual(values=c('#3C5488FF', '#DC0000FF')) +
  stat_ellipse(linetype = 'dashed') + theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey')) + geom_text(x=-2.3,y=3.8,label='k=2, stress=0.208', size=2)
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_bray.pdf', width=6, height=4)

pdf('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_stressplot.pdf')
stressplot(nmds)
dev.off()

merged_table_meta = full_metadata[full_metadata$Sample %in% colnames(merged_table),]
rownames(merged_table_meta) = merged_table_meta$Sample
# Statistical testing
div <- adonis2(t(merged_table) ~ Ecosystem + Dataset, data = merged_table_meta, permutations = 999, method="bray")
div # p < 0.001 for datasets (r2 = 0.027) and ecosystems (r2 = 0.201)
#               Df  SumOfSqs      R2      F       Pr(>F)    
#    Ecosystem   3   45.265  0.20135      43.742  0.001 ***
#    Dataset     1    6.032  0.02683      17.486  0.001 ***
#    Residual  503  173.507  0.77182                  
#    Total     507  224.803  1.00000   

pair_a = pairwise.adonis2(t(merged_table) ~ Ecosystem + Dataset , data = merged_table_meta, permutations = 999, method="bray")
pair_a # p < 0.001 for all comparisons, for datasets and ecosystems

r2_df = data.frame(Group1=c(),Group2=c(),r2=c())

# $`Snow/Ice_vs_Marine`
#               Df SumsOfSqs MeanSqs  F.Model      R2  Pr(>F)    
#   Ecosystem   1    14.331  14.3308  41.890  0.15622  0.001 ***
#   Dataset     1     4.876   4.8761  14.253  0.05316  0.001 ***
#   Residuals 212    72.526   0.3421          0.79062                 
#   Total     214    91.733                   1.00000           
r2_df = rbind(r2_df, data.frame(Group1='Snow/Ice',Group2='Marine',r2=0.15622))
# 
# $`Snow/Ice_vs_Terrestrial`
#               Df SumsOfSqs MeanSqs  F.Model      R2  Pr(>F)    
#   Ecosystem   1    15.121  15.1209  46.732  0.14506  0.001 ***
#   Dataset     1     5.639   5.6393  17.428  0.05410  0.001 ***
#   Residuals 258    83.481   0.3236          0.80084           
#   Total     260   104.241                   1.00000
r2_df = rbind(r2_df, data.frame(Group1='Snow/Ice',Group2='Terrestrial',r2=0.14506))
# 
# $`Snow/Ice_vs_Freshwater`
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    10.109 10.1092 26.4271 0.09031  0.001 ***
# Dataset     1     2.753  2.7533  7.1977 0.02460  0.001 ***
# Residuals 259    99.076  0.3825         0.88509           
# Total     261   111.938                 1.00000     
r2_df = rbind(r2_df, data.frame(Group1='Snow/Ice',Group2='Freshwater',r2=0.09031))
# 
# $Marine_vs_Terrestrial
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    20.101 20.1014  69.297 0.20551  0.001 ***
# Dataset     1     7.220  7.2203  24.891 0.07382  0.001 ***
# Residuals 243    70.489  0.2901         0.72067           
# Total     245    97.811                 1.00000  
r2_df = rbind(r2_df, data.frame(Group1='Marine',Group2='Terrestrial',r2=0.20551))
# 
# $Marine_vs_Freshwater
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    14.341 14.3414  41.164 0.13690  0.001 ***
# Dataset     1     5.408  5.4082  15.523 0.05162  0.001 ***
# Residuals 244    85.010  0.3484         0.81148           
# Total     246   104.759                 1.00000    
r2_df = rbind(r2_df, data.frame(Group1='Marine',Group2='Freshwater',r2=0.13690))
# 
# $Terrestrial_vs_Freshwater
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    16.488 16.4883  49.936 0.13900  0.001 ***
# Dataset     1     6.382  6.3815  19.327 0.05380  0.001 ***
# Residuals 290    95.755  0.3302         0.80721           
# Total     292   118.625                 1.00000   
r2_df = rbind(r2_df, data.frame(Group1='Terrestrial',Group2='Freshwater',r2=0.13900))

ggplot(r2_df) + geom_tile(aes(x=Group1, y=Group2, fill=r2)) + xlab('') + ylab('') + scale_fill_viridis_c() + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_r2_pairs.pdf', width=3.2, height=2.5)




# Bray-Curtis distances
snow_ice_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Snow/Ice']]
freshwater_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']]
marine_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']]
terrestrial_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']]

library(reshape2)
beta_df = melt(as.matrix(vegdist(t(subset), method = 'bray')))
PP1_snow_ice = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Snow/Ice') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Snow/Ice') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_snow_ice$Ecosystem = 'Snow/Ice'
PP1_terrestrial = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_terrestrial$Ecosystem = 'Terrestrial'
PP1_marine = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_marine$Ecosystem = 'Marine'
PP1_freshwater = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_freshwater$Ecosystem = 'Freshwater'

PP2_snow_ice = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Snow/Ice') & (merged_table_meta$Dataset == 'PP2')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Snow/Ice') & (merged_table_meta$Dataset == 'PP2')]),]
PP2_snow_ice$Ecosystem = 'Snow/Ice'
PP2_terrestrial = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP2')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP2')]),]
PP2_terrestrial$Ecosystem = 'Terrestrial'
PP2_marine = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP2')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP2')]),]
PP2_marine$Ecosystem = 'Marine'
PP2_freshwater = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP2')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP2')]),]
PP2_freshwater$Ecosystem = 'Freshwater'

beta_box = do.call("rbind",list(PP1_freshwater,PP1_marine,PP1_terrestrial,PP1_snow_ice,PP2_freshwater,PP2_marine,PP2_terrestrial,PP2_snow_ice))

ggplot(beta_box) + geom_boxplot(aes(fill=Ecosystem,y=value,x=Ecosystem)) + theme_linedraw() + ylab('Bray-Curtis distance') + xlab('')
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_bc_distances_eco.pdf')




########################################################################################
# Core analysis
# Calculating the binomial for the proba of occurrence
prob_occ <- function(gen_table, metadata, genus, ab_thr){
  subset = as.data.frame(t(as.matrix(gen_table[rownames(gen_table) == genus,])))
  
  if (sum(subset) > 0){
  subset[subset < ab_thr] = 0
  subset[subset > 0] = 1
  subset$Ecosystem = vapply(1:nrow(subset), function(x){as.character(full_metadata$Ecosystem[full_metadata$Sample == rownames(subset)[x]])}, FUN.VALUE = character(1))
  subset$Dataset = vapply(1:nrow(subset), function(x){as.character(full_metadata$Dataset[full_metadata$Sample == rownames(subset)[x]])}, FUN.VALUE = character(1))
  colnames(subset)[colnames(subset) == genus] = 'genus'
  fit = glm(genus ~ Ecosystem + Dataset, data=subset, family = binomial())
  pred_df = expand.grid(Ecosystem=c('Terrestrial','Snow/Ice','Marine','Freshwater'), Dataset=c('PP1','PP2'))
  pred_df$pred = predict(fit, newdata=pred_df, type='response')
  return(data.frame(Genus = genus, Ab_thr = ab_thr, P_cryosphere = mean(pred_df$pred), P_snow_ice = mean(pred_df$pred[pred_df$Ecosystem == 'Snow/Ice']), P_terrestrial = mean(pred_df$pred[pred_df$Ecosystem == 'Terrestrial']), P_marine = mean(pred_df$pred[pred_df$Ecosystem == 'Marine']), P_freshwater = mean(pred_df$pred[pred_df$Ecosystem == 'Freshwater'])))}
  
  else{return(data.frame(Genus = genus, Ab_thr = ab_thr, P_cryosphere = NA, P_snow_ice = NA, P_terrestrial = NA, P_marine = NA, P_freshwater = NA))}}

seq_log <- function(v1, v2, n) {exp(seq(from = log(v1), to = log(v2), length.out = n))}

registerDoMC(10)
proba_df = foreach(abundance=seq_log(1/5000,1,30), .combine='rbind') %:% foreach(genus=rownames(merged_table), .combine='rbind') %dopar% {prob_occ(merged_table, full_metadata, genus, abundance)}
proba_df = na.omit(proba_df)

# min sample counts
min(colSums(PP1_genus_tab[,colnames(PP1_genus_tab) %in% PP1_metadata$Sample[PP1_metadata$Cryosphere == 'Yes']])) # 3549
min(colSums(PP2_genus_tab[,colnames(PP2_genus_tab) %in% PP2_metadata$Sample[PP2_metadata$Cryosphere == 'Yes']])) # 1847

# Core definition
core_size = expand.grid(Abundance=seq_log(1/1847,1,30),Prevalence=seq_log(1/1847,1,30))
core_size$N = vapply(1:nrow(core_size), function(x){length(unique(proba_df$Genus[(proba_df$Ab_thr > core_size$Abundance[x]) & (proba_df$P_cryosphere > core_size$Prevalence[x])]))}, FUN.VALUE = numeric(1))
ggplot(core_size, aes(x=Abundance,y=Prevalence)) + geom_tile(aes(fill=N)) + 
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)), aes(x=x,y=y), color='darkgrey', linetype = "dashed") + geom_point(aes(x=0.0008685272, y=0.2), color='red',size=5) +
  scale_x_log10() + scale_y_log10() + scale_fill_gradient(name = 'Core size', trans = "log", breaks = c(1,10,100,1000),na.value = 'transparent') + theme_linedraw()
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_size.pdf', width=6, height=5)

proba_df$DA = 'Others'
proba_df$DA[proba_df$Genus %in% deseq_res$X[deseq_res$DA == 'Overrepresented']] = 'Cryo.'

proba_df$Core = 'Ancillary'
proba_df$Core[proba_df$Genus %in% proba_df$Genus[(proba_df$P_cryosphere > 0.2) & (proba_df$Ab_thr > 0.001)]] = 'Core'
length(unique(proba_df$Genus[proba_df$Core=='Core'])) # 52
# removing from the core taxa not detected in all of the cryospheric ecosystem
proba_df$Core[(proba_df$P_snow_ice < 0.0001) | (proba_df$P_terrestrial < 0.0001) | (proba_df$P_marine < 0.0001) | (proba_df$P_freshwater < 0.0001)] = 'Ancillary'
length(unique(proba_df$Genus[proba_df$Core=='Core'])) # 40

ggplot(proba_df) + geom_line(aes(x=Ab_thr,y=P_cryosphere,group=Genus, color=Core), alpha=0.5) + theme(panel.grid = element_line(colour = 'darkgrey')) +
  scale_x_log10() + theme_linedraw() + ylab('Prob[Presence] in the cryosphere') + xlab('Abundance threshold') + geom_hline(yintercept = 0.2, linetype='dashed') + geom_vline(xintercept = 0.001, linetype='dashed')
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_prob_presence.pdf', width=6, height=5)
# core defined at ab_thr > 0.001 present in at least 20% of the samples in the cryosphere, we removed 12 genera that were not detected in all ecosystems but only in 3/4

# Prevalence of the core & overrepresented genera in the cryospheric ecosystems
preval_hm = expand.grid(Taxonomy=unique(proba_df$Genus[proba_df$Core == 'Core']), Ecosystem=c('Terrestrial','Snow/Ice','Marine','Freshwater'))
GenusPreval <- function(genus, ecosystem, proba_df, thr){
  if (ecosystem == 'Terrestrial'){return(max(proba_df$P_terrestrial[(proba_df$Genus == genus) & (proba_df$Ab_thr >= thr)]))}
  else if (ecosystem == 'Marine'){return(max(proba_df$P_marine[(proba_df$Genus == genus) & (proba_df$Ab_thr >= thr)]))}
  else if (ecosystem == 'Freshwater'){return(max(proba_df$P_freshwater[(proba_df$Genus == genus) & (proba_df$Ab_thr >= thr)]))}
  else if (ecosystem == 'Snow/Ice'){return(max(proba_df$P_snow_ice[(proba_df$Genus == genus) & (proba_df$Ab_thr >= thr)]))}}
preval_hm$Prevalence = vapply(1:nrow(preval_hm), function(x){GenusPreval(preval_hm$Taxonomy[x], preval_hm$Ecosystem[x], proba_df, 1e-07)}, FUN.VALUE = numeric(1))
preval_hm$Phylum = vapply(preval_hm$Taxonomy, function(x) gsub('d__Bacteria; p__', '', strsplit(as.character(x), split='; c__')[[1]][1]), FUN.VALUE = character(1))
preval_hm$Genus = vapply(preval_hm$Taxonomy, function(x) strsplit(as.character(x),split = '; g__')[[1]][2], FUN.VALUE = character(1))

preval_hm$Phylum[!(preval_hm$Phylum %in% c('Proteobacteria','Bacteroidota','Actinobacteriota','Acidobacteriota','Verrucomicrobiota','Chloroflexi','Planctomycetota'))] = 'Others'

ggplot(preval_hm, aes(x=Ecosystem, y=Genus, size=Prevalence, color=Phylum)) + theme_linedraw() + labs(size='Prob(Presence)') +
  scale_color_manual(values = c('#B09C85FF', '#F39B7FFF', '#DC0000FF', '#91D1C2FF', 'grey', '#4DBBD5FF', '#3C5488FF', '#8491B4FF')) + 
  theme(panel.grid = element_line(colour = 'darkgrey'), strip.text = element_blank(), axis.text.y = element_text(face = "italic")) +
  geom_point() + facet_grid(Phylum~., scales = 'free_y', space = 'free_y') + xlab('') + ylab('')
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_genera_prevalence_eco.pdf', width = 7, height = 8)

########################################################################################
# Core and overrepresented genera overlap
core_genera = unique(proba_df$Genus[(proba_df$Core == 'Core')])
ancom_res$CORE = vapply(1:nrow(ancom_res), function(x) ifelse(ancom_res$X[x] %in% core_genera, 'Core', 'Ancillary'), FUN.VALUE = character(1))
table(ancom_res$CORE, ancom_res$DA)
#             Others  Overrepresented
# Ancillary   1472    558
# Core        11      29

fisher.test(table(ancom_res$DA, ancom_res$CORE))
# p-value = 6.383e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval: 3.341659 15.517322
# sample estimates: odds ratio = 6.946169 

mat_df = data.frame('DA'=c('Others','Cryo.','Others','Cryo.'), 'CORE'=c('Core','Core','Ancillary','Ancillary'), 'N'=c(11,29,1472,558))
mat_df$CORE <- factor(mat_df$CORE, levels = c('Core','Ancillary'))
ggplot(mat_df) + 
  geom_tile(aes(y=DA,x=CORE, alpha=N),fill='#F0C21A',show.legend = FALSE) + geom_text(aes(y=DA,x=CORE, label=N)) + xlab('') + ylab('') + scale_alpha_continuous(trans = "log") +
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_Over_intersect.pdf', width = 2.5, height = 1.5)

########################################################################################
# Core genera relative abundance
# Relative abundance of the core and DA genera in the different ecosystems
eco_abund = expand.grid(Dataset=c('PP1','PP2'), Ecosystem=c('Terrestrial','Snow/Ice','Freshwater','Marine'), DA = c('Others','Cryo.'), CORE = c('Core','Ancillary'))
ancom_res$DA = as.character(ancom_res$DA)
ancom_res$DA[ancom_res$DA == 'Overrepresented'] = 'Cryo.'
eco_abund$Abundance = vapply(1:nrow(eco_abund), function(x) mean(colSums(merged_table[rownames(merged_table) %in% ancom_res$X[(ancom_res$CORE == eco_abund$CORE[x]) & (ancom_res$DA == eco_abund$DA[x])],colnames(merged_table) %in% full_metadata$Sample[(full_metadata$Dataset == eco_abund$Dataset[x]) & (full_metadata$Ecosystem == eco_abund$Ecosystem[x])]])), numeric(1))

ggplot(eco_abund) + 
  geom_tile(aes(alpha=Abundance, y=DA, x=CORE, fill=Dataset)) + geom_text(aes( y=DA, x=CORE, label=paste0(round(Abundance,3)*100,'%'))) +
  theme_linedraw() + scale_fill_manual(values=c('#3C5488FF', '#DC0000FF')) + xlab('') + facet_grid(Ecosystem ~ Dataset) + ylab('') +
  theme(panel.grid = element_line(colour = 'darkgrey')) + guides(alpha = FALSE, fill = FALSE) 
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_Over_abundance.pdf', width=4, height=5)

# Taxonomic tree of the core microbiome
library(metacoder)
core_data = parse_tax_data(data.frame(Genus = ancom_res$taxa_id[ancom_res$CORE == 'Core']),
                           class_cols = "Genus", # the column that contains taxonomic information
                           class_sep = "; ", # The character used to separate taxa in the classification
                           class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c(tax_rank = "info", # A key describing each regex capture group
                                         tax_name = "taxon_name"))

heat_tree(core_data,
          node_label = taxon_names,
          node_size = n_subtaxa,
          node_color = n_supertaxa, 
          overlap_avoidance = 1.5,
          aspect_ratio = 1.1,
          node_size_axis_label = "Subtaxa N",
          node_color_axis_label = "Supertaxa N",
          node_color_range = c('#F39B7FFF','#8491B4FF'),
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          output_file = '2_Genus_analysis/2_3_Core_microbiome/2_3_Core_tax_tree.pdf')


# Comparison of the ecosystems overlap in most prevalent genera
core_snow_ice = as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.1)]))
core_terrestrial = as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.1)]))
core_freshwater = as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.1)]))
core_marine = as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.1)]))

eco_top_df = data.frame(Genus = union(core_snow_ice,union(core_terrestrial,union(core_marine,core_freshwater))))
eco_top_df$`Snow/Ice` = ifelse(eco_top_df$Genus %in% core_snow_ice, 1, 0)
eco_top_df$Marine = ifelse(eco_top_df$Genus %in% core_marine, 1, 0)
eco_top_df$Freshwater = ifelse(eco_top_df$Genus %in% core_freshwater, 1, 0)
eco_top_df$Terrestrial = ifelse(eco_top_df$Genus %in% core_terrestrial, 1, 0)

pdf('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_intersect.pdf', width=5, height = 5)
upset(eco_top_df, intersect = c('Snow/Ice','Freshwater','Terrestrial','Marine'), base_annotations=list('Intersection size'=intersection_size(counts=FALSE,fill='#F0C21A')))
dev.off()

ggplot(proba_df) + geom_density(aes(x=P_marine)) + geom_density(aes(x=P_terrestrial)) + geom_density(aes(x=P_freshwater)) + geom_density(aes(x=P_snow_ice)) + scale_x_log10() + theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))  
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Prevalence_thr_for_presence.pdf')



########################################################################################
# Alpha Diversity analysis
obs_gen_sample <- function(table, metadata, ab_thr, sample){
  subset = table[,colnames(table)==sample]
  subset[subset < ab_thr] = 0
  subset[subset > 0] = 1
  return(sum(subset))}

shannon_i <- function(table, metadata,sample){
  subset = table[,colnames(table)==sample]
  subset = vapply(na.omit(subset[subset > 0]), function(x){x * log(x)},FUN.VALUE = numeric(1))
  return(-sum(subset))}

merged_table_meta = full_metadata[full_metadata$Sample %in% colnames(merged_table),]
merged_table_meta$obs_genera = vapply(merged_table_meta$Sample, function(x){obs_gen_sample(merged_table, merged_table_meta, 1/1847,x)}, FUN.VALUE = numeric(1))
merged_table_meta$shannon_i = vapply(merged_table_meta$Sample, function(x){shannon_i(merged_table, merged_table_meta,x)}, FUN.VALUE = numeric(1))
sample_sizes <- merged_table_meta %>% group_by(Ecosystem, Dataset) %>% tally()
sample_sizes$n = paste0('N=',sample_sizes$n)

library(GFD)
merged_table_meta$Dataset = as.factor(merged_table_meta$Dataset)
model = GFD::GFD(shannon_i ~ Ecosystem * Dataset, data = merged_table_meta)
summary(model)
# Wald-Type Statistic (WTS):
# WTest                 statistic    df      p-value p-value WTPS
# WEcosystem            112.4423405  3  0.000000e+00       0.0000
# WDataset              0.4832448    1  4.869562e-01       0.4866
# WEcosystem:Dataset    53.7331209   3  1.279177e-11       0.0000

# WANOVA-Type Statistic (ATS):
# W                      Test statistic      df1        df2       p-value
# WEcosystem             36.9602555          2.811997   159.7556  0.000000e+00
# WDataset               0.4832448           1.000000   159.7556  4.879661e-01
# WEcosystem:Dataset     17.1786135          2.811997   159.7556  2.574276e-09

plot(model, factor = "Ecosystem", legendpos = "topleft", col = 3:5, pch = 17)

ggplot(merged_table_meta) + geom_boxplot(aes(y=obs_genera,x=Ecosystem)) + theme_linedraw() + xlab('') + ylab('Observed genera') +
  theme(panel.grid = element_line(colour = 'darkgrey')) + 
  geom_text(data = sample_sizes, aes(x=Ecosystem, fill=Dataset, label = n, y=250), vjust = 1, position = position_dodge(width = 1))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Div_observed.pdf')

ggplot(merged_table_meta) + geom_boxplot(aes(y=shannon_i,x=Ecosystem)) + theme_linedraw() + xlab('') + ylab("Shannon's Index") +
  theme(panel.grid = element_line(colour = 'darkgrey')) + 
  geom_text(data = sample_sizes, aes(x=Ecosystem, fill=Dataset, label = n, y=5), vjust = 1, position = position_dodge(width = 1))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Div_shannon.pdf')





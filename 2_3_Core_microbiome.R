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

setwd('/Users/mabourqu/Desktop/cryobiome_revisions/')
############################################################################################################
# Data laoding
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep = '\t')
PP1_metadata$Dataset = 'PP1'
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep = '\t')
PP2_metadata$Dataset = 'PP2'

PP1_genus_tab = read.csv('Data/PP1_genus.csv')
PP2_genus_tab = read.csv('Data/PP2_genus.csv')
PP1_metadata$seq_n = vapply(PP1_metadata$Sample, function(x) sum(as.integer(PP1_genus_tab[,colnames(PP1_genus_tab) == x])), FUN.VALUE = numeric(1))
PP2_metadata$seq_n = vapply(PP2_metadata$Sample, function(x) sum(as.integer(PP2_genus_tab[,colnames(PP2_genus_tab) == x])), FUN.VALUE = numeric(1))
# keep only samples with 2000 reads assigned to genus-taxonomy level ASVs
PP1_metadata = PP1_metadata[PP1_metadata$seq_n >= 2000,]
PP2_metadata = PP2_metadata[PP2_metadata$seq_n >= 2000,]

full_metadata = rbind(PP1_metadata,PP2_metadata)

PP1_data = PP1_genus_tab[,colnames(PP1_genus_tab) %in% PP1_metadata$Sample[(PP1_metadata$Sample %in% colnames(PP1_genus_tab)) & (PP1_metadata$Cryosphere == 'Yes')]]
rownames(PP1_data) = PP1_genus_tab$Genus
PP1_data$X = NULL
PP1_data$Genus = NULL
PP2_data = PP2_genus_tab[,colnames(PP2_genus_tab) %in% PP2_metadata$Sample[(PP2_metadata$Sample %in% colnames(PP2_genus_tab)) & (PP2_metadata$Cryosphere == 'Yes')]]
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

nmds = metaMDS(log(t(merged_table)+1), distance = "bray", k=2, trymax = 1000)
nmds_data = as.data.frame(scores(nmds))
nmds_data$Sample = colnames(merged_table)
nmds_data$Ecosystem = vapply(nmds_data$Sample, function(x){as.character(full_metadata$Ecosystem[full_metadata$Sample == x])}, FUN.VALUE = character(1))
nmds_data$Dataset = vapply(nmds_data$Sample, function(x){as.character(full_metadata$Dataset[full_metadata$Sample == x])}, FUN.VALUE = character(1))
nmds$stress # 0.206265

ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, shape = Ecosystem)) + geom_point(aes(colour = Dataset)) + scale_color_manual(values=c('#3C5488FF', '#DC0000FF')) +
  stat_ellipse(linetype = 'dashed') + theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey')) + geom_text(x=-2.1,y=3,label='k=2, stress=0.206', size=2)
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_bray.pdf', width=6, height=4)

pdf('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_stressplot.pdf')
stressplot(nmds)
dev.off()

merged_table_meta = full_metadata[full_metadata$Sample %in% colnames(merged_table),]
rownames(merged_table_meta) = merged_table_meta$Sample
# Statistical testing
div <- adonis2(t(merged_table) ~ Ecosystem + Dataset, data = merged_table_meta, permutations = 999, method="bray")
div # p < 0.001 for datasets (r2 = 0.01846) and ecosystems (r2 = 0.18269)
#             Df SumOfSqs       R2      F  Pr(>F)    
# Ecosystem   3    55.774  0.18319 52.702  0.001 ***
# Dataset     1     5.633  0.01850 15.969  0.001 ***
# Residual  689   243.055  0.79831                  
# Total     693   304.462  1.00000  

pair_a = pairwise.adonis2(t(merged_table) ~ Ecosystem + Dataset , data = merged_table_meta, permutations = 999, method="bray")
pair_a # p < 0.001 for all comparisons, for datasets and ecosystems

r2_df = data.frame(Group1=c(),Group2=c(),r2=c())
r2_df = rbind(r2_df, data.frame(Group1='Ice/Snow', Group2='Terrestrial', mean_sqs=18.8587))
# $`Ice/Snow_vs_Terrestrial`
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    18.859 18.8587   54.36 0.11779  0.001 ***
# Dataset     1     4.906  4.9056   14.14 0.03064  0.001 ***
# Residuals 393   136.340  0.3469         0.85157           
# Total     395   160.105                 1.00000

r2_df = rbind(r2_df, data.frame(Group1='Ice/Snow', Group2='Marine', mean_sqs=23.2800))
# $`Ice/Snow_vs_Marine`
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    23.280 23.2800  66.918 0.15758  0.001 ***
# Dataset     1     4.087  4.0873  11.749 0.02767  0.001 ***
# Residuals 346   120.370  0.3479         0.81476           
# Total     348   147.737                 1.00000   

r2_df = rbind(r2_df, data.frame(Group1='Ice/Snow',Group2='Freshwater',mean_sqs=10.9912))
# $`Ice/Snow_vs_Freshwater`
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    10.991 10.9912 28.6254 0.06768  0.001 ***
# Dataset     1     3.195  3.1947  8.3202 0.01967  0.001 ***
# Residuals 386   148.211  0.3840         0.91265           
# Total     388   162.397                 1.00000   

r2_df = rbind(r2_df, data.frame(Group1='Marine',Group2='Terrestrial',mean_sqs=25.4016))
# $Terrestrial_vs_Marine
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    25.402 25.4016  85.019 0.20705  0.001 ***
# Dataset     1     7.052  7.0523  23.604 0.05748  0.001 ***
# Residuals 302    90.230  0.2988         0.73547           
# Total     304   122.684                 1.00000

r2_df = rbind(r2_df, data.frame(Group1='Marine',Group2='Freshwater',mean_sqs=16.2002))
# $Marine_vs_Freshwater
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    18.464 18.4637  53.295 0.14665  0.001 ***
# Dataset     1     5.240  5.2400  15.125 0.04162  0.001 ***
# Residuals 295   102.202  0.3464         0.81173           
# Total     297   125.906                 1.00000 

r2_df = rbind(r2_df, data.frame(Group1='Terrestrial',Group2='Freshwater',mean_sqs=18.4637))
# $Terrestrial_vs_Freshwater
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Ecosystem   1    16.200 16.2002  47.332 0.11536  0.001 ***
# Dataset     1     7.175  7.1750  20.963 0.05109  0.001 ***
# Residuals 342   117.056  0.3423         0.83355           
# Total     344   140.431                 1.00000




ggplot(r2_df) + geom_tile(aes(x=Group1, y=Group2, fill=mean_sqs)) + xlab('') + ylab('') + scale_fill_viridis_c() + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_nmds_mean_squares_pairs.pdf', width=3.2, height=2.5)


# Bray-Curtis distances
snow_ice_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Ice/Snow']]
freshwater_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']]
marine_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']]
terrestrial_tab = merged_table[,colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']]
full_tab = as.matrix(cbind(snow_ice_tab,freshwater_tab,marine_tab,terrestrial_tab))

library(reshape2)
beta_df = melt(as.matrix(vegdist(t(full_tab), method = 'bray')))
PP1_snow_ice = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Ice/Snow') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Ice/Snow') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_snow_ice$Ecosystem = 'Ice/Snow'
PP1_terrestrial = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Terrestrial') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_terrestrial$Ecosystem = 'Terrestrial'
PP1_marine = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Marine') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_marine$Ecosystem = 'Marine'
PP1_freshwater = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP1')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Freshwater') & (merged_table_meta$Dataset == 'PP1')]),]
PP1_freshwater$Ecosystem = 'Freshwater'

PP2_snow_ice = beta_df[(beta_df$Var1 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Ice/Snow') & (merged_table_meta$Dataset == 'PP2')]) & (beta_df$Var2 %in% merged_table_meta$Sample[(merged_table_meta$Ecosystem == 'Ice/Snow') & (merged_table_meta$Dataset == 'PP2')]),]
PP2_snow_ice$Ecosystem = 'Ice/Snow'
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
  pred_df = expand.grid(Ecosystem=c('Terrestrial','Ice/Snow','Marine','Freshwater'), Dataset=c('PP1','PP2'))
  pred_df$pred = predict(fit, newdata=pred_df, type='response')
  return(data.frame(Genus = genus, Ab_thr = ab_thr, P_cryosphere = mean(pred_df$pred), P_snow_ice = mean(pred_df$pred[pred_df$Ecosystem == 'Ice/Snow']), P_terrestrial = mean(pred_df$pred[pred_df$Ecosystem == 'Terrestrial']), P_marine = mean(pred_df$pred[pred_df$Ecosystem == 'Marine']), P_freshwater = mean(pred_df$pred[pred_df$Ecosystem == 'Freshwater'])))}
  
  else{return(data.frame(Genus = genus, Ab_thr = ab_thr, P_cryosphere = NA, P_snow_ice = NA, P_terrestrial = NA, P_marine = NA, P_freshwater = NA))}}

seq_log <- function(v1, v2, n) {exp(seq(from = log(v1), to = log(v2), length.out = n))}

registerDoMC(10)
proba_df = foreach(abundance=seq_log(1/5000,1,30), .combine='rbind') %:% foreach(genus=rownames(merged_table), .combine='rbind') %dopar% {prob_occ(merged_table, full_metadata, genus, abundance)}
proba_df = na.omit(proba_df)

# min sample counts
min(colSums(PP1_genus_tab[,colnames(PP1_genus_tab) %in% PP1_metadata$Sample[PP1_metadata$Cryosphere == 'Yes']])) # 3549
min(colSums(PP2_genus_tab[,colnames(PP2_genus_tab) %in% PP2_metadata$Sample[PP2_metadata$Cryosphere == 'Yes']])) # 3141

# Core definition
core_size = expand.grid(Abundance=seq_log(1/3141,1,30),Prevalence=seq_log(1/3141,1,30))
core_size$N = vapply(1:nrow(core_size), function(x){length(unique(proba_df$Genus[(proba_df$Ab_thr > core_size$Abundance[x]) & (proba_df$P_cryosphere > core_size$Prevalence[x])]))}, FUN.VALUE = numeric(1))
ggplot(core_size, aes(x=Abundance,y=Prevalence)) + geom_tile(aes(fill=N)) + 
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)), aes(x=x,y=y), color='darkgrey', linetype = "dashed") + geom_point(aes(x=0.0008685272, y=0.2), color='red',size=5) +
  scale_x_log10() + scale_y_log10() + scale_fill_gradient(name = 'Core size', trans = "log", breaks = c(1,10,100,1000),na.value = 'transparent') + theme_linedraw()
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_size.pdf', width=6, height=5)

proba_df$DA = 'Others'
proba_df$DA[proba_df$Genus %in% ancom_res$X[ancom_res$DA == 'Overrepresented']] = 'Cryo.'

proba_df$Core = 'Ancillary'
proba_df$Core[proba_df$Genus %in% proba_df$Genus[(proba_df$P_cryosphere > 0.2) & (proba_df$Ab_thr > 0.001)]] = 'Core'
length(unique(proba_df$Genus[proba_df$Core=='Core'])) # 56
# removing from the core taxa not detected in all of the cryospheric ecosystem
proba_df$Core[(proba_df$P_snow_ice < 0.0001) | (proba_df$P_terrestrial < 0.0001) | (proba_df$P_marine < 0.0001) | (proba_df$P_freshwater < 0.0001)] = 'Ancillary'
length(unique(proba_df$Genus[proba_df$Core=='Core'])) # 37

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
# Ancillary   1445    562
# Core        10      27

fisher.test(table(ancom_res$DA, ancom_res$CORE))
# p-value = 2.544e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval: 3.227802 16.159243
# sample estimates: odds ratio = 6.933864 

mat_df = data.frame('DA'=c('Others','Cryo.','Others','Cryo.'), 'CORE'=c('Core','Core','Ancillary','Ancillary'), 'N'=c(10,27,1445,562))
mat_df$CORE <- factor(mat_df$CORE, levels = c('Core','Ancillary'))
ggplot(mat_df) + 
  geom_tile(aes(y=DA,x=CORE, alpha=N),fill='#F0C21A',show.legend = FALSE) + geom_text(aes(y=DA,x=CORE, label=N)) + xlab('') + ylab('') + scale_alpha_continuous(trans = "log") +
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_Over_intersect.pdf', width = 2.5, height = 1.5)

########################################################################################
# Core genera relative abundance
# Relative abundance of the core and DA genera in the different ecosystems
eco_abund = expand.grid(Dataset=c('PP1','PP2'), Ecosystem=c('Terrestrial','Ice/Snow','Freshwater','Marine'), DA = c('Others','Cryo.'), CORE = c('Core','Ancillary'))
ancom_res$DA = as.character(ancom_res$DA)
ancom_res$DA[ancom_res$DA == 'Overrepresented'] = 'Cryo.'
eco_abund$Abundance = vapply(1:nrow(eco_abund), function(x) mean(colSums(merged_table[rownames(merged_table) %in% ancom_res$X[(ancom_res$CORE == eco_abund$CORE[x]) & (ancom_res$DA == eco_abund$DA[x])],colnames(merged_table) %in% full_metadata$Sample[(full_metadata$Dataset == eco_abund$Dataset[x]) & (full_metadata$Ecosystem == eco_abund$Ecosystem[x])]])), numeric(1))

ggplot(eco_abund) + 
  geom_tile(aes(alpha=Abundance, y=DA, x=CORE, fill=Dataset)) + geom_text(aes( y=DA, x=CORE, label=paste0(round(Abundance,3)*100,'%'))) +
  theme_linedraw() + scale_fill_manual(values=c('#3C5488FF', '#DC0000FF')) + xlab('') + facet_grid(Ecosystem ~ Dataset) + ylab('') +
  theme(panel.grid = element_line(colour = 'darkgrey')) + guides(alpha = FALSE, fill = FALSE) 
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_Over_abundance.pdf', width=4, height=5)

# The core of the ice/snow in other cryospheric ecosystems
PP1_genus_tab = read.csv('Data/PP1_genus.csv')
PP2_genus_tab = read.csv('Data/PP2_genus.csv')
PP1_other_data = PP1_genus_tab[,colnames(PP1_genus_tab) %in% PP1_metadata$Sample[(PP1_metadata$Sample %in% colnames(PP1_genus_tab)) & (PP1_metadata$Cryosphere == 'No')]]
rownames(PP1_other_data) = PP1_genus_tab$Genus
PP1_other_data$X = NULL
PP1_other_data$Genus = NULL
PP2_other_data = PP2_genus_tab[,colnames(PP2_genus_tab) %in% PP2_metadata$Sample[(PP2_metadata$Sample %in% colnames(PP2_genus_tab)) & (PP2_metadata$Cryosphere == 'No')]]
rownames(PP2_other_data) = PP2_genus_tab$Genus
PP2_other_data$X = NULL
PP2_other_data$Genus = NULL
merged_other_table = merge(PP1_other_data,PP2_other_data,by='row.names',all=TRUE)
rownames(merged_other_table) = merged_other_table$Row.names
merged_other_table$Row.names = NULL
merged_other_table[is.na(merged_other_table)] = 0
merged_other_table = sweep(merged_other_table,2,colSums(merged_other_table),'/')

eco_abund_others = expand.grid(Dataset=c('PP1','PP2'), Ecosystem=c('Terrestrial','Freshwater','Marine','Ice/Snow'), DA = c('Others','Cryo.'), CORE = c('Core','Ancillary'))
ancom_res$DA = as.character(ancom_res$DA)
ancom_res$DA[ancom_res$DA == 'Overrepresented'] = 'Cryo.'
ancom_res$CORE = 'Ancillary'
ancom_res$CORE[ancom_res$X %in% proba_df$Genus[proba_df$Core=='Core']] = 'Core'
eco_abund_others$Abundance = vapply(1:nrow(eco_abund_others), function(x) mean(colSums(merged_other_table[rownames(merged_other_table) %in% ancom_res$X[(ancom_res$CORE == eco_abund_others$CORE[x]) & (ancom_res$DA == eco_abund_others$DA[x])],colnames(merged_other_table) %in% full_metadata$Sample[(full_metadata$Dataset == eco_abund_others$Dataset[x]) & (full_metadata$Ecosystem == eco_abund_others$Ecosystem[x])]])), numeric(1))

mean(eco_abund_others$Abundance[(eco_abund_others$DA == 'Others') & (eco_abund_others$CORE == 'Ancillary')], na.rm = T) # 0.7303518
mean(eco_abund_others$Abundance[(eco_abund_others$DA == 'Others') & (eco_abund_others$CORE == 'Core')], na.rm = T) # 0.03740859
mean(eco_abund_others$Abundance[(eco_abund_others$DA == 'Cryo.') & (eco_abund_others$CORE == 'Ancillary')], na.rm = T) # 0.1536504
mean(eco_abund_others$Abundance[(eco_abund_others$DA == 'Cryo.') & (eco_abund_others$CORE == 'Core')], na.rm = T) # 0.06825328


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
core_snow_ice = as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)]))
core_terrestrial = as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)]))
core_freshwater = as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)]))
core_marine = as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)]))

eco_top_df = data.frame(Genus = union(core_snow_ice,union(core_terrestrial,union(core_marine,core_freshwater))))
eco_top_df$`Ice/Snow` = ifelse(eco_top_df$Genus %in% core_snow_ice, 1, 0)
eco_top_df$Marine = ifelse(eco_top_df$Genus %in% core_marine, 1, 0)
eco_top_df$Freshwater = ifelse(eco_top_df$Genus %in% core_freshwater, 1, 0)
eco_top_df$Terrestrial = ifelse(eco_top_df$Genus %in% core_terrestrial, 1, 0)

pdf('2_Genus_analysis/2_3_Core_microbiome/2_3_Core_intersect.pdf', width=5, height = 5)
upset(eco_top_df, intersect = c('Ice/Snow','Freshwater','Terrestrial','Marine'), base_annotations=list('Intersection size'=intersection_size(counts=FALSE,fill='#F0C21A')))
dev.off()

ggplot(proba_df) + geom_density(aes(x=P_marine)) + geom_density(aes(x=P_terrestrial)) + geom_density(aes(x=P_freshwater)) + geom_density(aes(x=P_snow_ice)) + scale_x_log10() + theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'))  
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Prevalence_thr_for_presence.pdf')


# abundance of snow/ice core (20% prevalence thr.) in other cryospheric ecosystems
length(as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)]))) # 78
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Ice/Snow']], na.rm = T)) # 0.6862067
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']], na.rm = T)) # 0.02720503
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']], na.rm = T)) # 0.2594294
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_snow_ice > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']], na.rm = T)) # 0.2449747

length(as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)]))) # 134
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Ice/Snow']], na.rm = T)) # 0.3723059
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']], na.rm = T)) # 0.05472161
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']], na.rm = T)) # 0.3162096
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_terrestrial > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']], na.rm = T)) # 0.730957

length(as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)]))) # 72
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Ice/Snow']], na.rm = T)) # 0.3106653
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']], na.rm = T)) # 0.09880818
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']], na.rm = T)) # 0.5314461
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_freshwater > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']], na.rm = T)) # 0.2024502

length(as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)]))) # 128
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Ice/Snow']], na.rm = T)) # 0.07335711
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Marine']], na.rm = T)) # 0.8809816
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Freshwater']], na.rm = T)) # 0.120688
mean(colSums(merged_table[rownames(merged_table) %in% as.character(unique(proba_df$Genus[(proba_df$P_marine > 0.2) & (proba_df$Ab_thr >= 0.001)])), colnames(merged_table) %in% merged_table_meta$Sample[merged_table_meta$Ecosystem == 'Terrestrial']], na.rm = T)) # 0.02567832


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
model = GFD::GFD(shannon_i ~ Ecosystem + Dataset:Ecosystem, data = merged_table_meta)
summary(model)
# Call: 
# shannon_i ~ Ecosystem + Dataset:Ecosystem
# 
# Descriptive:
#   Ecosystem Dataset   n    Means Variances Lower 95 % CI Upper 95 % CI
# 1  Freshwater     PP1  29 3.723131 0.2619025      3.528467      3.917796
# 5  Freshwater     PP2 140 2.842509 0.4813978      2.726569      2.958449
# 2    Ice/Snow     PP1  92 2.760108 0.6895707      2.588136      2.932080
# 6    Ice/Snow     PP2 128 2.925377 0.6782041      2.781337      3.069416
# 3      Marine     PP1  88 3.028792 0.1780076      2.939398      3.118186
# 7      Marine     PP2  41 3.709338 0.4086679      3.507559      3.911116
# 4 Terrestrial     PP1  92 3.693610 0.7978177      3.508633      3.878588
# 8 Terrestrial     PP2  84 3.642556 0.2837374      3.526959      3.758152
weighted.mean(c(3.723131, 2.835726), c(29,141)) # 2.987107
weighted.mean(c(2.760108, 2.925377), c(92,128)) # 2.856265
weighted.mean(c(3.028792, 3.709338), c(88,41))  # 3.24509
weighted.mean(c(3.693610, 3.642556), c(92,84))  # 3.669243
# 
# Wald-Type Statistic (WTS):
#                   Test statistic df p-value p-value WTPS
# Ecosystem               112.0236  3       0            0
# Ecosystem:Dataset       103.1681  4       0            0
# 
# ANOVA-Type Statistic (ATS):
#   Test statistic      df1      df2 p-value
# Ecosystem               37.78357 2.998591 331.5813       0
# Ecosystem:Dataset       25.72954 3.997182 331.5813       0

pdf('2_Genus_analysis/2_3_Core_microbiome/2_3_Div_shannon.pdf')
plot(model, factor = "Ecosystem", legendpos = "topleft", col = 3:5, pch = 17)
dev.off()

ggplot(merged_table_meta) + geom_boxplot(aes(y=obs_genera,x=Ecosystem, fill = Dataset)) + theme_linedraw() + xlab('') + ylab('Observed genera') +
  theme(panel.grid = element_line(colour = 'darkgrey'))
ggsave('2_Genus_analysis/2_3_Core_microbiome/2_3_Div_observed.pdf')







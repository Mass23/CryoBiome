library(ggplot2)
library(dplyr)
library(stringr)
library(gt)
library(formattable)
library(webshot)
library(htmltools)
webshot::install_phantomjs(force = T)

# Data loading
setwd('/Users/mabourqu/Desktop/cryobiome_revisions')
stats_file = read.csv('Data/Functional_clusters/Unassigned_clusters_stats.tsv', sep='\t')
clusters_file = read.csv('Data/Functional_clusters/clusters_min_30_seq_2_samp.tsv', sep='\t')
diamond_uniprot = read.csv('Data/Functional_clusters/diamond_clusters.tsv',sep='\t',header = FALSE)
metadata = read.csv('Metadata/MTG_metadata.tsv')

Cryosphere_samples = as.character(metadata$Sample[metadata$Cryosphere == 'Yes'])

##################################################################################################
# Append other data to the stats file
diamond_uniprot$V1 = gsub("_", "", diamond_uniprot$V1)

clusters_trna = read.csv('Data/Functional_clusters/cluster_tRNA_collection.txt',sep='\t', header = F)
clusters_trna$V1 = gsub('_tRNA.txt','',clusters_trna$V1)
clusters_trna$V1 = gsub('_','',clusters_trna$V1)

clusters_file$KEGG = str_split_fixed(clusters_file$Sequence, "_", 3)[,1]
clusters_file$Accession = str_split_fixed(clusters_file$Sequence, "_", 3)[,2]
clusters_file$SeqID =  str_split_fixed(clusters_file$Sequence, "_", 3)[,3]
stats_file$ClusterID = gsub('_','',stats_file$ClusterID)
stats_file$tRNA = 'No'
stats_file$Uniprot_match = 'No_match'

for (cluster in unique(stats_file$ClusterID)){
  cluster = gsub("_", "", cluster)
  print(cluster)
  # tRNA
  if (cluster %in% clusters_trna$V1){stats_file$tRNA[stats_file$ClusterID == cluster] <- as.character(clusters_trna$V2[clusters_trna$V1 == cluster])}
  
  # Proportion of Cryosphere sequences
  sample_list = clusters_file$Accession[clusters_file$ClusterID == cluster]
  n_cryo_samples = sum(sample_list %in% Cryosphere_samples)
  n_other_samples = length(sample_list) - n_cryo_samples

  
  if ((n_other_samples == 0) & (n_cryo_samples > 0)){stats_file$Cryosphere[stats_file$ClusterID == cluster] = 'Cryosphere'}
  else if ((n_other_samples > 0) & (n_cryo_samples == 0)){stats_file$Cryosphere[stats_file$ClusterID == cluster] = 'non-Cryosphere'}
  else if ((n_other_samples > 0) & (n_cryo_samples > 0)) {stats_file$Cryosphere[stats_file$ClusterID == cluster] <- 'Shared'}
  
  # Prop of assigned KEGG in cluster
  kegg_list = unique(clusters_file$KEGG[clusters_file$ClusterID == cluster])
  unassigned = FALSE
  if ('Unassigned' %in% kegg_list){unassigned = TRUE}
  
  if ((length(kegg_list) == 1) & (unassigned == FALSE)){stats_file$Unassigned[stats_file$ClusterID == cluster] <- 'KEGG'}
  else if((length(kegg_list) > 1)){stats_file$Unassigned[stats_file$ClusterID == cluster] <- 'Ambiguous'}
  else {stats_file$Unassigned[stats_file$ClusterID == cluster] <- 'Unassigned'}
  
  # Uniprot matches
  if (cluster %in% diamond_uniprot$V1){
    stats_file$Uniprot_match[stats_file$ClusterID == cluster] <- 'Uniprot_match'
    stats_file$Uniprot_identity[stats_file$ClusterID == cluster] <- diamond_uniprot$V3[diamond_uniprot$V1 == cluster]
  }}


############################## STATS SUMMARY ############################## 
stat_summary = data.frame(Annotation=c(), Category=c(), N_clusters=c(), Uniprot=c(), tRNA=c())
stat_summary = rbind(stat_summary, data.frame(Annotation='KEGG', Category='Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'KEGG'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'KEGG') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='KEGG', Category='Shared', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'KEGG'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'KEGG') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='KEGG', Category='Non-Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'KEGG'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'KEGG') & (stats_file$Uniprot_match != 'No_match'),])))

stat_summary = rbind(stat_summary, data.frame(Annotation='Ambiguous', Category='Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'Ambiguous'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'Ambiguous') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='Ambiguous', Category='Shared', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'Ambiguous'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'Ambiguous') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='Ambiguous', Category='Non-Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'Ambiguous'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'Ambiguous') & (stats_file$Uniprot_match != 'No_match'),])))

stat_summary = rbind(stat_summary, data.frame(Annotation='Unassigned', Category='Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'Unassigned'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Cryosphere') & (stats_file$Unassigned == 'Unassigned') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='Unassigned', Category='Shared', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'Unassigned'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'Shared') & (stats_file$Unassigned == 'Unassigned') & (stats_file$Uniprot_match != 'No_match'),])))
stat_summary = rbind(stat_summary, data.frame(Annotation='Unassigned', Category='Non-Cryosphere', N_clusters = nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'Unassigned'),]), Uniprot=nrow(stats_file[(stats_file$Cryosphere == 'non-Cryosphere') & (stats_file$Unassigned == 'Unassigned') & (stats_file$Uniprot_match != 'No_match'),])))

stat_summary$Uniprot = round((stat_summary$Uniprot / stat_summary$N_clusters) * 100,2)
colnames(stat_summary)[colnames(stat_summary) == 'N_clusters'] = 'Number of clusters'
colnames(stat_summary)[colnames(stat_summary) == 'Uniprot'] = 'Uniprot match [%]'
t = formattable(stat_summary, list('Category' = formatter("span", style = x ~ style(color = ifelse(x == 'Cryosphere', '#23A671', ifelse(x == 'Shared', '#A673D1', 'grey')),font.weight = "bold"))))#, 
                                #area(row = c(1,4), col = 'Uniprot match [%]') ~ color_tile("white","#23A671"),
                                #area(row = c(2,5), col = 'Uniprot match [%]') ~ color_tile("white","#A673D1"),
                                #area(row = c(3,6), col = 'Uniprot match [%]') ~ color_tile("white","lightgrey"),
                                #area(row = c(1,4), col = 'tRNA') ~ color_tile("white","#23A671"),
                                #area(row = c(2,5), col = 'tRNA') ~ color_tile("white","#A673D1"),
                                #area(row = c(3,6), col = 'tRNA') ~ color_tile("white","lightgrey")))
w <- as.htmlwidget(t, width = '600', height = '200')
path <- html_print(w)
webshot(url = path, file = '3_Functional_analysis/3_3_Unknown_space/Unknown_func_space_summary.jpg', selector = ".html-widget", expand = 100, zoom = 20)

############################## Unknown ############################## 
library(ggpubr)
stats_file$Cryosphere[stats_file$Cryosphere == 'Cryosphere'] = 'Cryo.'
stats_file$Cryosphere[stats_file$Cryosphere == 'non-Cryosphere'] = 'Others'
ggplot(stats_file,aes(x=as.factor(Cryosphere), fill=as.factor(Cryosphere),y=as.numeric(Uniprot_identity))) + geom_boxplot() + guides(fill=FALSE) + 
  facet_grid(~Unassigned) + theme_linedraw() + ylab('Uniprot match identity [%]') + xlab('')  + scale_fill_manual(values = c('#23A671','grey','#A673D1')) + 
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Cryo.","Others"), c("Cryo.","Shared"), c("Shared","Others")))  + theme(panel.grid = element_line(color='grey'))
ggsave('3_Functional_analysis/3_3_Unknown_space/Unknown_uniprot.pdf', width = 5, height = 4)

ggplot(stats_file,aes(x=as.factor(Cryosphere), fill=as.factor(Cryosphere),y=meanGC)) + geom_boxplot() + guides(fill=FALSE) + 
  facet_grid(~Unassigned) + theme_linedraw() + ylab('GC content [%]') + xlab('')  + scale_fill_manual(values = c('#23A671','grey','#A673D1')) + 
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Cryo.","Others"), c("Cryo.","Shared"), c("Shared","Others")))  + theme(panel.grid = element_line(color='grey'))
ggsave('3_Functional_analysis/3_3_Unknown_space/Unknown_meanGC.pdf', width = 5, height = 4)

ggplot(stats_file, aes(x=as.factor(Cryosphere), fill=as.factor(Cryosphere),y=meanPID)) + geom_boxplot() + guides(fill=FALSE) + 
  facet_grid(~Unassigned) + theme_linedraw() + ylab('Mean identity [%]') + xlab('')  + scale_fill_manual(values = c('#23A671','grey','#A673D1')) + 
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Cryo.","Others"), c("Cryo.","Shared"), c("Shared","Others"))) + theme(panel.grid = element_line(color='grey'))
ggsave('3_Functional_analysis/3_3_Unknown_space/Unknown_meanId.pdf', width = 5, height = 4)

# Logistics for the number of shared clusters
stats_file$Cryosphere <- relevel(as.factor(stats_file$Cryosphere),ref = "Others")
stats_file$Unassigned <- relevel(as.factor(stats_file$Unassigned),ref = "KEGG")
uniprot_m = glm(as.factor(Uniprot_match) ~ Cryosphere + Unassigned, data = stats_file, family = binomial())
summary(uniprot_m)
# Deviance Residuals: 
# Min       1Q   Median       3Q      Max  
# -1.3820  -1.1056  -0.7021   1.0472   1.7508  
# Coefficients:
#                      Estimate    Std. Error  z value  Pr(>|z|)    
# (Intercept)          0.46905     0.03652     12.843   < 2e-16 ***
# CryosphereCryo.      -0.65501    0.12386     -5.288   1.24e-07 ***
# CryosphereShared     -0.64040    0.04093     -15.646  < 2e-16 ***
# UnassignedAmbiguous  -0.15486    0.04515     -3.430   0.000603 ***
# UnassignedUnassigned -1.10333    0.04843     -22.784  < 2e-16 ***
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 16738  on 12124  degrees of freedom
# Residual deviance: 15791  on 12120  degrees of freedom
# AIC: 15801
# Number of Fisher Scoring iterations: 4
ors = as.data.frame(exp(cbind(coef(uniprot_m),confint(uniprot_m))))
colnames(ors) = c('Estimate', 'low', 'high')
ors$Group = c('Absent in the cryosphere and KEGG annotation', 'Detected in the cryosphere', 'Specific to the cryosphere', 'Ambiguous annotation', 'Unassigned annotation')
ors$facets = c('Intercept', 'Cryosphere', 'Cryosphere', 'Annotation','Annotation')
ggplot(ors) + geom_point(aes(y=Group, x=Estimate), size = 2) + geom_errorbarh(aes(y=Group, xmax = high, xmin = low)) + facet_grid(facets~., space = 'free_y', scales = 'free_y') +
  geom_vline(aes(xintercept = 1), linetype = 'dashed') + theme_linedraw() + ylab('') + xlab('Odds ratio (of Uniprot match)') + theme(panel.grid = element_line(color='grey'))
ggsave('3_Functional_analysis/3_3_Unknown_space/Uniprot_match_OR.pdf', width = 5, height = 4)




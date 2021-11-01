library(phytools)
library(phangorn)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(foreach)
library(doMC)
library(rstatix)
library(coin)
library(tidyverse)
library(dplyr)
setwd('/Users/mabourqu/Documents/PhD/C1/')

library(devtools)
#install_github("Russel88/MicEco")
library(MicEco)
library(picante)
####################################################################################################
PP1_tab_abund = read.csv('Data/trees/PP1_005_ASVs_table.csv')
PP1_ASV_coefs = read.csv('Data/PP1_Logistics_coefs.csv')
PP1_metadata = read.csv('Metadata/PP1_metadata.tsv', sep ='\t')
PP1_tree = read.newick(file='Data/trees/PP1_005_trim.fasta.treefile')
PP1_tree = midpoint(PP1_tree)
PP1_tab = PP1_tab_abund[,colnames(PP1_tab_abund) %in% PP1_metadata$Sample]
rownames(PP1_tab) = PP1_tab_abund$ASV

PP2_tab_abund = read.csv('Data/trees/PP2_005_ASVs_table.csv')
PP2_ASV_coefs = read.csv('Data/PP2_Logistics_coefs.csv')
PP2_metadata = read.csv('Metadata/PP2_metadata.tsv', sep ='\t')
PP2_tree = read.newick(file='Data/trees/PP2_005_trim.fasta.treefile')
PP2_tree = midpoint(PP2_tree)
PP2_tab = PP2_tab_abund[,colnames(PP2_tab_abund) %in% PP2_metadata$Sample]
rownames(PP2_tab) = PP2_tab_abund$ASV

registerDoMC(10)
PP1_cryo_samples = PP1_metadata$Sample[PP1_metadata$Cryosphere == 'Yes']
PP1_other_samples = PP1_metadata$Sample[PP1_metadata$Cryosphere == 'No']
PP2_cryo_samples = PP2_metadata$Sample[PP2_metadata$Cryosphere == 'Yes']
PP2_other_samples = PP2_metadata$Sample[PP2_metadata$Cryosphere == 'No']

PP1_metadata$Ecosystem[PP1_metadata$Ecosystem %in% c('Snow', 'Ice')] = 'Snow/Ice'
PP2_metadata$Ecosystem[PP2_metadata$Ecosystem %in% c('Snow', 'Ice')] = 'Snow/Ice'

####################################################################################################
# 1. Sorensen's PD
calc_sor <- function(table, tree, cryo_s, other_s, n_samp){
  cryo_sub = sample(cryo_s, n_samp)
  other_sub = sample(other_s, n_samp)
  samples = c(cryo_sub, other_sub)
  tab = table[,colnames(table) %in% samples]
  tab = tab[rowSums(tab) > 0,]
  tree = prune.sample(t(tab), tree)
  sor_pd = phylosor(t(tab), tree)
  melted_sor_pd = melt(as.matrix(sor_pd), varnames = c("S1", "S2"))
  melted_sor_pd$Group = vapply(1:nrow(melted_sor_pd), function(x){
    if ((melted_sor_pd$S1[x] %in% cryo_sub) & (melted_sor_pd$S2[x] %in% cryo_sub)){return('Cryo-Cryo')}
    else if ((melted_sor_pd$S1[x] %in% cryo_sub) | (melted_sor_pd$S2[x] %in% cryo_sub)){return('Cryo-Others')}
    else {return('Others-Others')}}, FUN.VALUE = character(1))
  return(melted_sor_pd)}


# PP1
PP1_sor_pd = foreach(i=1:50, .combine = 'rbind') %dopar% {calc_sor(PP1_tab, PP1_tree, PP1_cryo_samples, PP1_other_samples, 50)}
PP1_sor_pd = PP1_sor_pd[!duplicated(PP1_sor_pd),]
min_n_samples = min(c(sum(PP1_sor_pd$Group == 'Cryo-Cryo'), sum(PP1_sor_pd$Group == 'Cryo-Others'), sum(PP1_sor_pd$Group == 'Others-Others')))
PP1_sor_cc = PP1_sor_pd[PP1_sor_pd$Group == 'Cryo-Cryo',]
PP1_sor_co = PP1_sor_pd[PP1_sor_pd$Group == 'Cryo-Others',]
PP1_sor_oo = PP1_sor_pd[PP1_sor_pd$Group == 'Others-Others',]
PP1_sor_pd = rbind(PP1_sor_cc[sample(1:nrow(PP1_sor_cc), min_n_samples),], 
                   PP1_sor_co[sample(1:nrow(PP1_sor_co), min_n_samples),], 
                   PP1_sor_oo[sample(1:nrow(PP1_sor_oo), min_n_samples),])
write.csv(PP1_sor_pd, 'Data/PP1_sor_PD.csv', quote = F, row.names = F)

PP1_sor_pd = read.csv('Data/PP1_sor_PD.csv')
PP1_sor_pd %>% group_by(Group) %>% summarise(median = median(value)) 
# Cryo-Cryo     0.265
# Cryo-Others   0.234
# Others-Others 0.225

ggplot(PP1_sor_pd, aes(y=value,x=Group,fill=Group)) + geom_boxplot() + 
  theme_linedraw() + scale_fill_manual(values = c('#3C5488FF','#8491B4FF','#91D1C2FF')) + ylim(c(0,1.35)) +
  xlab('') + theme(axis.text.x = element_blank(), legend.position = "none") + ylab(expression(beta~'-Sor. PD')) +
  stat_compare_means(aes(label = ..p.signif..), comparisons=list(c(1,2), c(1,3), c(2,3)))
ggsave('1_ASV_analysis/1_4_Phylogenetics/1_PP1_sor_PD.pdf', width=2.5,height = 3)
kruskal.test(PP1_sor_pd$value, PP1_sor_pd$Group)
# Kruskal-Wallis chi-squared = 2500.8, df = 2, p-value < 2.2e-16
compare_means(data=PP1_sor_pd, formula = value ~ Group)
# .y.     group1      group2        p        p.adj  p.format  p.signif method  
# 1 value Cryo-Cryo   Cryo-Others   0        0      <2e-16    ****     Wilcoxon
# 2 value Cryo-Cryo   Others-Others 0        0      <2e-16    ****     Wilcoxon
# 3 value Cryo-Others Others-Others 5.47e-16 5.5e-16 5.5e-16  ****     Wilcoxon
PP1_sor_pd %>% wilcox_effsize(value ~ Group)
#   .y.   group1      group2         effsize    n1    n2    magnitude
# 1 value Cryo-Cryo   Cryo-Others    0.173      60456 60456 small
# 2 value Cryo-Cryo   Others-Others  0.120      60456 60456 small   
# 3 value Cryo-Others Others-Others  0.0233     60456 60456 small  


# PP2
PP2_sor_pd = foreach(i=1:50, .combine = 'rbind') %dopar% {calc_sor(PP2_tab, PP2_tree, PP2_cryo_samples, PP2_other_samples, 50)}
PP2_sor_pd = PP2_sor_pd[!duplicated(PP2_sor_pd),]
min_n_samples = min(c(sum(PP2_sor_pd$Group == 'Cryo-Cryo'), sum(PP2_sor_pd$Group == 'Cryo-Others'), sum(PP2_sor_pd$Group == 'Others-Others')))
PP2_sor_cc = PP2_sor_pd[PP2_sor_pd$Group == 'Cryo-Cryo',]
PP2_sor_co = PP2_sor_pd[PP2_sor_pd$Group == 'Cryo-Others',]
PP2_sor_oo = PP2_sor_pd[PP2_sor_pd$Group == 'Others-Others',]
PP2_sor_pd = rbind(PP2_sor_cc[sample(1:nrow(PP2_sor_cc), min_n_samples),], 
                   PP2_sor_co[sample(1:nrow(PP2_sor_co), min_n_samples),], 
                   PP2_sor_oo[sample(1:nrow(PP2_sor_oo), min_n_samples),])
write.csv(PP2_sor_pd, 'Data/PP2_sor_PD.csv', quote = F, row.names = F)

PP2_sor_pd = read.csv('Data/PP2_sor_PD.csv')
PP2_sor_pd %>% group_by(Group) %>% summarise(median = median(value)) 
# Cryo-Cryo     0.265
# Cryo-Others   0.212
# Others-Others 0.193

ggplot(PP2_sor_pd, aes(y=value,x=Group,fill=Group)) + geom_boxplot() + ylim(c(0,1.35)) +
  theme_linedraw() + scale_fill_manual(values = c('#DC0000FF','#F39B7FFF','#91D1C2FF')) +
  xlab('') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab(expression(beta~'-Sor. PD')) +
  stat_compare_means(aes(label = ..p.signif..), comparisons=list(c(1,2), c(1,3), c(2,3)))
ggsave('1_ASV_analysis/1_4_Phylogenetics/1_PP2_sor_PD.pdf', width=2.5,height = 3.7)
kruskal.test(PP2_sor_pd$value, PP2_sor_pd$Group)
# Kruskal-Wallis chi-squared = 12199, df = 2, p-value < 2.2e-16
compare_means(data=PP2_sor_pd, formula = value ~ Group)
# .y.     group1      group2        p         p.adj    p.format  p.signif  method  
# 1 value Cryo-Cryo   Cryo-Others   0         0        <2e-16    ****      Wilcoxon
# 2 value Cryo-Cryo   Others-Others 0         0        <2e-16    ****      Wilcoxon
# 3 value Cryo-Others Others-Others 3.64e-186 3.6e-186 <2e-16    ****      Wilcoxon
PP2_sor_pd %>% wilcox_effsize(value ~ Group)
#   .y.   group1      group2         effsize    n1    n2 magnitude
# 1 value Cryo-Cryo   Cryo-Others    0.217  73164 73164 small 
# 2 value Cryo-Cryo   Others-Others  0.269  73164 73164 small 
# 3 value Cryo-Others Others-Others  0.0761 73164 73164 small   




####################################################################################################
# 2. MNTD
calc_mntd <- function(table, tree, cryo_s, other_s, n_samp){
  cryo_sub = sample(cryo_s, n_samp)
  other_sub = sample(other_s, n_samp)
  samples = c(cryo_sub, other_sub)
  tab = table[,colnames(table) %in% samples]
  tab = tab[rowSums(tab) > 0,]
  tree = prune.sample(t(tab), tree)
  mntd = comdistnt(t(tab), cophenetic(tree), abundance.weighted = FALSE)
  melted_mntd = melt(as.matrix(mntd), varnames = c("S1", "S2"))
  melted_mntd$Group = vapply(1:nrow(melted_mntd), function(x){
    if ((melted_mntd$S1[x] %in% cryo_sub) & (melted_mntd$S2[x] %in% cryo_sub)){return('Cryo-Cryo')}
    else if ((melted_mntd$S1[x] %in% cryo_sub) | (melted_mntd$S2[x] %in% cryo_sub)){return('Cryo-Others')}
    else {return('Others-Others')}}, FUN.VALUE = character(1))
  return(melted_mntd)}

# PP1
PP1_mntd = foreach(i=1:50, .combine = 'rbind') %dopar% {calc_mntd(PP1_tab, PP1_tree, PP1_cryo_samples, PP1_other_samples, 50)}
PP1_mntd = PP1_mntd[!duplicated(PP1_mntd),]
min_n_samples = min(c(sum(PP1_mntd$Group == 'Cryo-Cryo'), sum(PP1_mntd$Group == 'Cryo-Others'), sum(PP1_mntd$Group == 'Others-Others')))
PP1_mntd_cc = PP1_mntd[PP1_mntd$Group == 'Cryo-Cryo',]
PP1_mntd_co = PP1_mntd[PP1_mntd$Group == 'Cryo-Others',]
PP1_mntd_oo = PP1_mntd[PP1_mntd$Group == 'Others-Others',]
PP1_mntd = rbind(PP1_mntd_cc[sample(1:nrow(PP1_mntd_cc), min_n_samples),], 
                 PP1_mntd_co[sample(1:nrow(PP1_mntd_co), min_n_samples),], 
                 PP1_mntd_oo[sample(1:nrow(PP1_mntd_oo), min_n_samples),])
write.csv(PP1_mntd, 'Data/PP1_mntd.csv', quote = F, row.names = F)

PP1_mntd = read.csv('Data/PP1_mntd.csv')
PP1_mntd %>% group_by(Group) %>% summarise(median = median(value)) 
# Cryo-Cryo     0.424
# Cryo-Others   0.439
# Others-Others 0.427

# PP2
PP2_mntd = foreach(i=1:50, .combine = 'rbind') %dopar% {calc_mntd(PP2_tab, PP2_tree, PP2_cryo_samples, PP2_other_samples, 50)}
PP2_mntd = PP2_mntd[!duplicated(PP2_mntd),]
min_n_samples = min(c(sum(PP2_mntd$Group == 'Cryo-Cryo'), sum(PP2_mntd$Group == 'Cryo-Others'), sum(PP2_mntd$Group == 'Others-Others')))
PP2_mntd_cc = PP2_mntd[PP2_mntd$Group == 'Cryo-Cryo',]
PP2_mntd_co = PP2_mntd[PP2_mntd$Group == 'Cryo-Others',]
PP2_mntd_oo = PP2_mntd[PP2_mntd$Group == 'Others-Others',]
PP2_mntd = rbind(PP2_mntd_cc[sample(1:nrow(PP2_mntd_cc), min_n_samples),], 
                 PP2_mntd_co[sample(1:nrow(PP2_mntd_co), min_n_samples),], 
                 PP2_mntd_oo[sample(1:nrow(PP2_mntd_oo), min_n_samples),])
write.csv(PP2_mntd, 'Data/PP2_mntd.csv', quote = F, row.names = F)

PP2_mntd = read.csv('Data/PP2_mntd.csv')
PP2_mntd %>% group_by(Group) %>% summarise(median = median(value)) 
# Cryo-Cryo     0.340
# Cryo-Others   0.371
# Others-Others 0.367

# Plotting
ggplot(PP1_mntd, aes(y=value,x=Group,fill=Group)) + geom_boxplot() + ylim(c(0,1.3)) +
  theme_linedraw() + scale_fill_manual(values = c('#3C5488FF','#8491B4FF','#91D1C2FF')) +
  xlab('') + theme(axis.text.x = element_blank(), legend.position = "none") + ylab(expression(beta~'-MNTD')) +
  stat_compare_means(aes(label = ..p.signif..), comparisons=list(c(1,2), c(1,3), c(2,3)))
ggsave('1_ASV_analysis/1_4_Phylogenetics/2_PP1_mntd.pdf', width=2.5,height = 3)
kruskal.test(PP1_mntd$value, PP1_mntd$Group)
# Kruskal-Wallis chi-squared = 1824.9, df = 2, p-value < 2.2e-16
compare_means(data=PP1_mntd, formula = value ~ Group)
# .y.     group1      group2        p          p.adj     p.format  p.signif method  
# 1 value Cryo-Cryo   Cryo-Others   0          0         <2e-16    ****     Wilcoxon
# 2 value Cryo-Cryo   Others-Others 7.85e-12   7.8e-2    7.8e-12   ****     Wilcoxon
# 3 value Cryo-Others Others-Others 71.71e-216 3.40e-216 <2e-16    ****     Wilcoxon
PP1_mntd %>% wilcox_effsize(value ~ Group)
#   .y.   group1      group2         effsize    n1    n2 magnitude
# 1 value Cryo-Cryo   Cryo-Others    0.123  56350 56350 small   
# 2 value Cryo-Cryo   Others-Others  0.0204 56350 56350 small     
# 3 value Cryo-Others Others-Others  0.0935 56350 56350 small

ggplot(PP2_mntd, aes(y=value,x=Group,fill=Group)) + geom_boxplot() + ylim(c(0,1.3)) +
  theme_linedraw() + scale_fill_manual(values = c('#DC0000FF','#F39B7FFF','#91D1C2FF')) +
  xlab('') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab(expression(beta~'-MNTD')) +
  stat_compare_means(aes(label = ..p.signif..), comparisons=list(c(1,2), c(1,3), c(2,3)))
ggsave('1_ASV_analysis/1_4_Phylogenetics/2_PP2_mntd.pdf', width=2.5,height = 3.7)
kruskal.test(PP2_mntd$value, PP2_mntd$Group)
# Kruskal-Wallis chi-squared = 5143.9, df = 2, p-value < 2.2e-16
compare_means(data=PP2_mntd, formula = value ~ Group)
#   .y.   group1      group2        p         p.adj      p.format p.signif method  
# 1 value Cryo-Cryo   Cryo-Others   0         0          <2e-16   ****     Wilcoxon
# 2 value Cryo-Cryo   Others-Others 0         0          <2e-16   ****     Wilcoxon
# 3 value Cryo-Others Others-Others 1.27e-37  1.3e-37    <2e-16   ****     Wilcoxon
PP2_mntd %>% wilcox_effsize(value ~ Group)
#   .y.   group1      group2         effsize    n1    n2 magnitude
# 1 value Cryo-Cryo   Cryo-Others    0.185  68682 68682 small      
# 2 value Cryo-Cryo   Others-Others  0.144  68682 68682 small     
# 3 value Cryo-Others Others-Others  0.0346 68682 68682 small 

####################################################################################################
# 3. alpha-Metrics
PP1_alpha = as.data.frame(pd(t(PP1_tab), PP1_tree, include.root = F))
PP2_alpha = as.data.frame(pd(t(PP2_tab), PP2_tree, include.root = F))
PP1_alpha$MPD = mpd(t(PP1_tab), cophenetic(PP1_tree), abundance.weighted = F)
PP2_alpha$MPD = mpd(t(PP2_tab), cophenetic(PP2_tree), abundance.weighted = F)
PP1_alpha$MNTD = mntd(t(PP1_tab), cophenetic(PP1_tree), abundance.weighted = F)
PP2_alpha$MNTD = mntd(t(PP2_tab), cophenetic(PP2_tree), abundance.weighted = F)

rownames(PP1_alpha) = colnames(PP1_tab)
rownames(PP2_alpha) = colnames(PP2_tab)
PP1_alpha$Cryosphere = vapply(rownames(PP1_alpha), function(x) PP1_metadata$Cryosphere[PP1_metadata$Sample == x], FUN.VALUE = character(1))
PP1_alpha$Ecosystem = vapply(rownames(PP1_alpha), function(x) PP1_metadata$Ecosystem[PP1_metadata$Sample == x], FUN.VALUE = character(1))
PP1_alpha$Dataset = 'PP1'
PP2_alpha$Cryosphere = vapply(rownames(PP2_alpha), function(x) PP2_metadata$Cryosphere[PP2_metadata$Sample == x], FUN.VALUE = character(1))
PP2_alpha$Ecosystem = vapply(rownames(PP2_alpha), function(x) PP2_metadata$Ecosystem[PP2_metadata$Sample == x], FUN.VALUE = character(1))
PP2_alpha$Dataset = 'PP2'
alpha_df = rbind(PP1_alpha, PP2_alpha)

ggplot(alpha_df, aes(y=MPD, x=log(SR), color = Cryosphere)) + facet_grid(Dataset~.) + geom_point() + geom_smooth(method = 'lm')

MPD_m = lm(MPD ~ Cryosphere + log(SR) + Dataset, data = na.omit(alpha_df))
summary(MPD_m)
# Residuals:
# Min       1Q   Median       3Q      Max 
# -0.45052 -0.07964 -0.00518  0.07446  0.65272 
# 
# Coefficients:
#                  Estimate    Std. Error  t value  Pr(>|t|)    
#   (Intercept)    0.490874    0.011985    40.956   <2e-16 ***
#   CryosphereYes  0.049507    0.005474    9.044    <2e-16 ***
#   log(SR)        0.048146    0.003460    13.917   <2e-16 ***
#   DatasetPP2    -0.075126    0.003721    -20.191  <2e-16 ***
# Residual standard error: 0.1145 on 3991 degrees of freedom
# Multiple R-squared:  0.1387,	Adjusted R-squared:  0.1381 
# F-statistic: 214.3 on 3 and 3991 DF,  p-value: < 2.2e-16

MNTD_m = lm(MNTD ~ Cryosphere + log(SR) + Dataset, data = na.omit(alpha_df))
summary(MNTD_m)
# Residuals:
# Min       1Q   Median       3Q      Max 
# -0.21661 -0.03887 -0.00707  0.03208  0.87751 
# 
# Coefficients:
#                  Estimate   Std. Error t value  Pr(>|t|)    
#   (Intercept)    0.341571   0.006835   49.973   < 2e-16 ***
#   CryosphereYes  0.006427   0.003122   2.059    0.0396 *  
#   log(SR)        -0.060753  0.001973   -30.793  < 2e-16 ***
#   DatasetPP2     0.012777   0.002122   6.021    1.88e-09 ***
# Residual standard error: 0.0653 on 3991 degrees of freedom
# Multiple R-squared:  0.1982,	Adjusted R-squared:  0.1976 
# F-statistic: 328.8 on 3 and 3991 DF,  p-value: < 2.2e-16

PD_m = lm(PD ~ Cryosphere + log(SR) + Dataset, data = na.omit(alpha_df))
summary(PD_m)
# Residuals:
# Min      1Q  Median      3Q     Max 
# -3.3278 -0.8812 -0.1408  0.7099  6.2195 
# 
# Coefficients:
#                 Estimate    Std. Error  t value  Pr(>|t|)    
#   (Intercept)   -5.33978    0.13353     -39.991  < 2e-16 ***
#   CryosphereYes  0.23945    0.06099     3.926    8.77e-05 ***
#   log(SR)        2.85003    0.03854     73.946   < 2e-16 ***
#   DatasetPP2     0.37850    0.04145     9.131    < 2e-16 ***
# Residual standard error: 1.276 on 3991 degrees of freedom
# Multiple R-squared:  0.5842,	Adjusted R-squared:  0.5839 
# F-statistic:  1869 on 3 and 3991 DF,  p-value: < 2.2e-16

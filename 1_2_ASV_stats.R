library(data.table)
library(ggplot2)
library(tibble)
library(tidyr)

setwd('SET_THE_WORKING_DIRECTORY_HERE')
###########################################################################################################################
# 1. Data loading, we load only the cryo samples data
# PP1
PP1_ASV_tab  = read.csv('Data/trees/PP1_005_ASVs_table.csv')
PP1_metadata = fread('Metadata/PP1_metadata.tsv')
PP1_cryo_samples = PP1_metadata$Sample[PP1_metadata$Cryosphere == 'Yes']
PP1_other_samples = PP1_metadata$Sample[PP1_metadata$Cryosphere == 'No']

# PP2
PP2_ASV_tab  = read.csv('Data/trees/PP2_005_ASVs_table.csv')
PP2_metadata = fread('Metadata/PP2_metadata.tsv')
PP2_cryo_samples = PP2_metadata$Sample[PP2_metadata$Cryosphere == 'Yes']
PP2_other_samples = PP2_metadata$Sample[PP2_metadata$Cryosphere == 'No']

###########################################################################################################################
# 2. Define Specific, Shared and Others ASVs.
PP1_ASV_tab$cryo_preval = rowSums(PP1_ASV_tab[,PP1_cryo_samples] > 0)
PP1_ASV_tab$Others_preval = rowSums(PP1_ASV_tab[,PP1_other_samples] > 0)
PP2_ASV_tab$cryo_preval = rowSums(PP2_ASV_tab[,PP2_cryo_samples] > 0)
PP2_ASV_tab$Others_preval = rowSums(PP2_ASV_tab[,PP2_other_samples] > 0)

PP1_ASV_tab$group = 'Others'
PP2_ASV_tab$group = 'Others'

PP1_ASV_tab$group[PP1_ASV_tab$cryo_preval > 0] = 'Shared'
PP2_ASV_tab$group[PP2_ASV_tab$cryo_preval > 0] = 'Shared'

PP1_ASV_tab$group[(PP1_ASV_tab$cryo_preval > 0) & (PP1_ASV_tab$Others_preval == 0)] = 'Cryosphere'
PP2_ASV_tab$group[(PP2_ASV_tab$cryo_preval > 0) & (PP2_ASV_tab$Others_preval == 0)] = 'Cryosphere'

write.csv(PP1_ASV_tab[,c('ASV','group')], file='1_ASV_analysis/1_2_ASV_stats/1_2_PP1_ASV_groups.csv',row.names=F)
write.csv(PP2_ASV_tab[,c('ASV','group')], file='1_ASV_analysis/1_2_ASV_stats/1_2_PP2_ASV_groups.csv',row.names=F)

# PP1
nrow(PP1_ASV_tab)                        # 19353
sum(PP1_ASV_tab$group == 'Others')       # 15647
sum(PP1_ASV_tab$group == 'Shared')       # 346
sum(PP1_ASV_tab$group == 'Cryosphere')   # 3360 

# PP2
nrow(PP2_ASV_tab)                        # 18296
sum(PP2_ASV_tab$group == 'Others')       # 13605 
sum(PP2_ASV_tab$group == 'Shared')       # 224 
sum(PP2_ASV_tab$group == 'Cryosphere')   # 4467 


ASV_number_df = data.frame('Dataset' = c('PP1','PP1','PP1','PP2','PP2','PP2'), 
                           'Group' = c('Others', 'Cryosphere','Shared', 'Others', 'Cryosphere', 'Shared'),
                           'Count' = c(15647,3360,346,13605,4467,224))

ggplot(ASV_number_df[ASV_number_df$Group != 'Others',]) + geom_bar(aes(x=Group,y=Count,fill=Dataset),stat='identity') + 
  scale_fill_manual(values = c("#3C5488FF", "#DC0000FF")) + theme_linedraw() + facet_grid(~Dataset) + xlab('') + ylab('Abund. (rel. ab. > 0.5%) ASV number') + guides(fill=FALSE)
ggsave('1_ASV_analysis/1_2_ASV_stats/1_2_ASV_numbers_cryo.pdf',width=5,height=5)

ggplot(ASV_number_df) + geom_bar(aes(x=Group,y=Count,fill=Dataset),stat='identity') + 
  scale_fill_manual(values = c("#3C5488FF", "#DC0000FF")) + theme_linedraw() + facet_grid(~Dataset) + xlab('') + ylab('Abund. (rel. ab. > 0.5%) ASV number') + guides(fill=FALSE)
ggsave('1_ASV_analysis/1_2_ASV_stats/1_2_ASV_numbers_all.pdf',width=5,height=5)

###########################################################################################################################
# 3. Contribution to the communities: we load full tables, to find ASV shared with other ecosystems
PP1_tab = fread('Data/PP1_table.tsv')
PP2_tab = fread('Data/PP2_table.tsv')

PP1_to_keep = c('ASV','Taxonomy',PP1_cryo_samples, PP1_other_samples)
PP1_rel_tab = as.data.frame(PP1_tab[,..PP1_to_keep])
PP2_to_keep = c('ASV','Taxonomy',PP2_cryo_samples, PP2_other_samples)
PP2_rel_tab = as.data.frame(PP2_tab[,..PP2_to_keep])

PP1_samples = c(PP1_cryo_samples, PP1_other_samples)
PP2_samples = c(PP2_cryo_samples, PP2_other_samples)
PP1_rel_tab[,PP1_samples] = apply(PP1_rel_tab[,PP1_samples],2,function(x){x/sum(x)})
PP2_rel_tab[,PP2_samples] = apply(PP2_rel_tab[,PP2_samples],2,function(x){x/sum(x)})

PP1_rel_tab$cryo_preval = rowSums(PP1_rel_tab[,PP1_cryo_samples] > 0)
PP1_rel_tab$Others_preval = rowSums(PP1_rel_tab[,PP1_other_samples] > 0)
PP2_rel_tab$cryo_preval = rowSums(PP2_rel_tab[,PP2_cryo_samples] > 0)
PP2_rel_tab$Others_preval = rowSums(PP2_rel_tab[,PP2_other_samples] > 0)

PP1_rel_tab$group = 'Others'
PP2_rel_tab$group = 'Others'

PP1_rel_tab$group[(PP1_rel_tab$cryo_preval > 0) & (PP1_rel_tab$Others_preval > 0)] = 'Shared'
PP2_rel_tab$group[(PP2_rel_tab$cryo_preval > 0) & (PP2_rel_tab$Others_preval > 0)] = 'Shared'

PP1_rel_tab$group[(PP1_rel_tab$cryo_preval > 0) & (PP1_rel_tab$Others_preval == 0)] = 'Cryosphere'
PP2_rel_tab$group[(PP2_rel_tab$cryo_preval > 0) & (PP2_rel_tab$Others_preval == 0)] = 'Cryosphere'

PP1_specific_abundance = colSums(PP1_rel_tab[PP1_rel_tab$group=='Cryosphere',PP1_cryo_samples]) # mean: 0.6240109
PP1_shared_abundance = colSums(PP1_rel_tab[PP1_rel_tab$group=='Shared',PP1_cryo_samples]) # mean: 0.3759891
nrow(PP1_rel_tab[PP1_rel_tab$group=='Cryosphere',PP1_cryo_samples])   # 99336
nrow(PP1_rel_tab[PP1_rel_tab$group=='Shared',PP1_cryo_samples]) # 12876

PP2_specific_abundance = colSums(PP2_rel_tab[PP2_rel_tab$group=='Cryosphere',PP2_cryo_samples]) # mean: 0.7795496
PP2_shared_abundance = colSums(PP2_rel_tab[PP2_rel_tab$group=='Shared',PP2_cryo_samples]) # mean: 0.2204504
nrow(PP2_rel_tab[PP2_rel_tab$group=='Cryosphere',PP2_cryo_samples])   # 76699
nrow(PP2_rel_tab[PP2_rel_tab$group=='Shared',PP2_cryo_samples]) # 4867

ASV_contribution_df = data.frame('Dataset' = c(rep('PP1',301),rep('PP1',301),rep('PP2',394),rep('PP2',394)), 
                                 'Group' = c(rep('Cryosphere',301),rep('Shared',301),rep('Cryosphere',394),rep('Shared',394)), 
                                 'Abundance' = c(PP1_specific_abundance,PP1_shared_abundance,PP2_specific_abundance,PP2_shared_abundance))

ggplot(ASV_contribution_df) + geom_boxplot(aes(x=Group,y=Abundance,fill=Dataset)) + 
  scale_fill_manual(values = c("#3C5488FF", "#DC0000FF")) + theme_linedraw() + facet_grid(~Dataset) + xlab('') + ylab('ASV number') + guides(fill=FALSE)
ggsave('1_ASV_analysis/1_2_ASV_stats/1_2_ASV_groups_relative_abundance.pdf',width=3,height=4)


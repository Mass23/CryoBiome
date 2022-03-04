library(ggplot2)
library(EnhancedVolcano)
library(DESeq2)

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

setwd('SET_THE_WORKING_DIRECTORY_HERE')
############################################################################################################
# 1. Loading data

MTG_KEGG_tab = read.csv('Data/MTG_KEGG_counts.tsv', sep = '\t')
rownames(MTG_KEGG_tab) = MTG_KEGG_tab$Geneid
MTG_KEGG_tab$Geneid = NULL
MTG_metadata = read.csv('Metadata/MTG_metadata.tsv', sep = '\t')
MTG_metadata = MTG_metadata[MTG_metadata$Sample %in% colnames(MTG_KEGG_tab),]

MTG_KEGG_tab = MTG_KEGG_tab[,colnames(MTG_KEGG_tab) %in% MTG_metadata$Sample]
write.csv(MTG_metadata, file = 'Metadata/MTG_metadata.tsv', quote = F, row.names = F)

############################################################################################################
# 2. DESeq2 analysis
MTG_metadata$Cryosphere = as.factor(MTG_metadata$Cryosphere)
MTG_metadata$Cryosphere <- relevel(MTG_metadata$Cryosphere, "Yes" )
KEGG_dds <- DESeqDataSetFromMatrix(countData=MTG_KEGG_tab+1, colData=MTG_metadata, design=~Cryosphere)
KEGG_deseq <- DESeq(KEGG_dds)
KEGG_res <- results(KEGG_deseq)
KEGG_res$padj[is.na(KEGG_res$padj)] = 1
KEGG_significant = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]

write.csv(as.data.frame(KEGG_res), file = 'Data/KEGG_deseq_results.csv')

############################################################################################################
# 1. DESeq2 volcano plot

EnhancedVolcano(KEGG_res,
                lab = rownames(KEGG_res),
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey30", "grey30", "#23A671"),
                x = 'log2FoldChange',
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                y = 'padj')
ggsave('3_Functional_analysis/3_1_Enrichment/KEGG_enriched_Cryosphere.pdf', width = 7, height = 7)

############################################################################################################
# 2. KEGG decoder
MTG_Cryosphere_samples = MTG_metadata$Sample[MTG_metadata$Cryosphere == 'Yes']

kegg_decoder_df = data.frame(Sample=c(), KEGG=c())
for (kegg_id in KEGG_significant) {
  for (sample in MTG_Cryosphere_samples) {
    if(MTG_KEGG_tab[rownames(MTG_KEGG_tab) == kegg_id, sample] > 0) {
      kegg_decoder_df = rbind(kegg_decoder_df, data.frame(Sample=sample, KEGG=kegg_id))
    }
  }
}
write.table(kegg_decoder_df, file='Data/KEGG_decoder_df_in.tsv', sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)

# in the Data folder using the CryoBiome-kegg-decoder env: KEGG-decoder -i KEGG_decoder_df_in.tsv -o KEGG_decoder_output

KEGG_decoder = read.table('Data/KEGG_decoder_out.tsv', sep='\t', header=TRUE)
colsums = colSums(KEGG_decoder[,!colnames(KEGG_decoder) %in% c('Function')])
KEGG_decoder = KEGG_decoder[,c('Function', names(colsums[colsums > 0]))]
names(KEGG_decoder)[names(KEGG_decoder)=="Function"] <- "Sample"
KEGG_decoder = reshape2::melt(KEGG_decoder, variable.name = 'Function', value.name = 'Completion')

# Renaming pathways for plotting
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="NAD.P.H.quinone.oxidoreductase"] <- "NAD(P)H quinone oxidoreductase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="chitinase"] <- "Chitinase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Cytochrome.bd.complex"] <- "Cytochrome-bd-complex"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Cytochrome.o.ubiquinol.oxidase"] <- "Cytochrome-o-ubiquinol-oxidase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="NiFe.hydrogenase.Hyd.1"] <- "NiFe hydrogenase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="hydrogen.quinone.oxidoreductase"] <- "Hydrogen-quinone-oxidoreductase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="sulfite.dehydrogenase"] <- "Sulfite-dehydrogenase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="nitrogen.fixation"] <- "Nitrogen fixation"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Photosystem.I"] <- "Photosystem I"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="alcohol.oxidase"] <- "Alcohol oxidase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Naphthalene.degradation.to.salicylate"] <- "Naphtalene degradation to salicylate"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Mixed.acid..Formate.to.CO2...H2"] <- "Mixed acid: Formate to CO2 & H2"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Mixed.acid..Acetate"] <- "Mixed acid: Acetate"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Entner.Doudoroff.Pathway"] <- "Entner-Doudoroff Pathway"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Competence.related.core.components"] <- "Competence related core components"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Biofilm.PGA.Synthesis.protein"] <- "Biofilm PGA-synthesis protein"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Type.III.Secretion"] <- "Type III Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Type.VI.Secretion"] <- "Type VI Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Type.IV.Secretion"] <- "Type IV Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Type.II.Secretion"] <- "Type II Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function)=="Fe.Mn.transporter.MntH"] <- "Fe/Mn transporter"

KEGG_decoder$mean_comp = vapply(1:nrow(KEGG_decoder), function(x) mean(KEGG_decoder$Completion[KEGG_decoder$Function == KEGG_decoder$Function[x]]), FUN.VALUE = numeric(1))
KEGG_decoder$Function <- factor(KEGG_decoder$Function, levels=unique((KEGG_decoder$Function)[order(KEGG_decoder$mean_comp)]))
KEGG_decoder$Ecosystem = vapply(1:nrow(KEGG_decoder), function(x) MTG_metadata$Ecosystem[MTG_metadata$Sample == KEGG_decoder$Sample[x]], FUN.VALUE = character(1))
ggplot(KEGG_decoder) +
  geom_tile(aes(x = as.factor(Sample), y = as.factor(Function), fill = Completion)) +
  xlab('') + ylab('') + theme_linedraw() + guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10)) +
  scale_fill_gradientn(colours = c('white','#00A087FF'), values = c(0,1)) + facet_grid(~Ecosystem, scales = 'free_x', space = 'free_x') +
  theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), panel.spacing = unit(0, "lines"), strip.background = element_blank())
ggsave('3_Functional_analysis/3_1_Enrichment/KEGG_Cryosphere_specific.pdf', width = 8, height = 5)

############################################################################################################
# 3. World map with samples used in this analysis
library(dplyr)
library(mapproj)
WorldData <- map_data('world') %>% fortify
ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id=region),
           fill = "white", colour = "grey10", size=0.2) +
  coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) +
  xlab('') + ylab('') +
  theme_bw() + 
  geom_point(data = na.omit(MTG_metadata), aes(x = Longitude, y = Latitude, colour = Cryosphere), size = 4) + 
  scale_colour_manual(values = c('#23A671','grey')) + theme(legend.position = 'none')
ggsave('3_Functional_analysis/3_1_Enrichment/Worldmap_functional.pdf', width = 7.5, height = 5)

############################################################################################################
# 4. DESeq2 to taxonomy
KEGG_sign_samp = kegg_decoder_df[kegg_decoder_df$KEGG %in% KEGG_significant,]
write.table(KEGG_sign_samp, file='Data/KEGG_sign_sample.tsv', sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)



import pandas as pd
import glob

kegg_files = glob.glob('Data/Raw/Functional/*KEGG_counts.tsv')
out_df = pd.DataFrame(columns=['Geneid'])

for file in kegg_files:
    with open(file, 'r') as kegg_file:
        data = pd.read_csv(kegg_file, sep = '\t')
        sample_name = file.split('/')[-1].replace('_KEGG_counts.tsv','')
        data = data[['Geneid', 'Assembly/mg.reads.sorted.bam']]
        data.rename(columns = {'Assembly/mg.reads.sorted.bam':sample_name}, inplace = True)
        out_df = pd.merge(out_df, data, on=['Geneid'], how='outer')

out_df = out_df.fillna(0)
out_df.to_csv('Data/MTG_KEGG_counts.tsv', sep = '\t', index = False)

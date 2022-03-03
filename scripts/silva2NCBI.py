import pandas as pd
import numpy as np

deseq_table = pd.read_csv('../Data/Deseq_results_amplicon.csv',sep=',')
print(deseq_table)
deseq_table[['Tax','Genus']] = deseq_table['Unnamed: 0'].str.split('g__',expand=True,)

ncbi_table = pd.DataFrame(columns=['silva_genus','ncbi_genus','DA','CORE'])

for i, r in deseq_table.iterrows():
    ncbi_genus = r['Genus'].replace('[','').replace(']','').replace('_',' ') 
    sign = r['DA']
    group= r['CORE']
    ncbi_table = ncbi_table.append({'silva_genus':r['Genus'],'ncbi_genus':ncbi_genus,'DA':r['DA'],'CORE':r['CORE']}, ignore_index=True)
    print(r['Genus'], ncbi_genus, r['DA'], r['CORE'])

ncbi_table.to_csv('../Data/deseq_ncbi_table.csv',sep=',')
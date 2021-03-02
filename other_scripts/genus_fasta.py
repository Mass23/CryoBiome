import pandas as pd
from Bio import SeqIO

PP1_tax = pd.read_csv('Data/Raw/PP1/PP1_raw_taxonomy.tsv', sep = '\t')
PP2_tax = pd.read_csv('Data/Raw/PP2/PP2_raw_taxonomy.tsv', sep = '\t')

# Keep onlyASVs assigned to genus
PP1_tax = PP1_tax[PP1_tax['Taxon'].str.contains('; g__')]
PP2_tax = PP2_tax[PP2_tax['Taxon'].str.contains('; g__')]
# Remove species level classification
PP1_tax['Taxon'] = PP1_tax['Taxon'].str.replace(r';\ s__\w+','',regex=True)
PP2_tax['Taxon'] = PP2_tax['Taxon'].str.replace(r';\ s__\w+','',regex=True)
# Extract different level data
PP1_tax[['domain', 'phylum', 'class', 'order',' family', 'genus']] = PP1_tax['Taxon'].str.split('; ', n=5, expand=True)
PP2_tax[['domain', 'phylum', 'class', 'order',' family', 'genus']] = PP2_tax['Taxon'].str.split('; ', n=5, expand=True)

PP1_genera = list(set(PP1_tax['genus'].to_list()))
PP1_genera = [i for i in PP1_genera if 'g__uncultured' not in i]
PP1_genera = [i for i in PP1_genera if 'g__Mitochondria' not in i]
PP1_genera = [i.replace('/','_').replace('.','') for i in PP1_genera]
PP2_genera = list(set(PP2_tax['genus'].to_list()))
PP2_genera = [i for i in PP2_genera if 'g__uncultured' not in i]
PP2_genera = [i for i in PP2_genera if 'g__Mitochondria' not in i]
PP2_genera = [i.replace('/','_').replace('.','') for i in PP2_genera]

print('----------------------\n          PP1\n----------------------')
PP1_seqs = list(SeqIO.parse('Data/Raw/PP1/PP1_raw_seqs.fasta', 'fasta'))
for genus in PP1_genera:
    ASV_genus = PP1_tax['Feature ID'][PP1_tax['genus'] == genus].to_list()
    if len(ASV_genus) > 9:
        print(genus + '\t' + str(len(ASV_genus)))
        SeqIO.write([seq for seq in PP1_seqs if seq.id in ASV_genus] ,'Data/Filtered/Genera/PP1_' + genus.lstrip('g__') + '.fasta' , 'fasta')

print('----------------------\n          PP2\n----------------------')
PP2_seqs = list(SeqIO.parse('Data/Raw/PP2/PP2_raw_seqs.fasta', 'fasta'))
for genus in PP2_genera:
    ASV_genus = PP2_tax['Feature ID'][PP2_tax['genus'] == genus].to_list()
    if len(ASV_genus) > 9:
        print(genus + '\t' + str(len(ASV_genus)))
        SeqIO.write([seq for seq in PP2_seqs if seq.id in ASV_genus] ,'Data/Filtered/Genera/PP2_' + genus.lstrip('g__') + '.fasta' , 'fasta')

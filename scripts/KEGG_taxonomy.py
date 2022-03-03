import pandas as pd
from Bio import Entrez
Entrez.email = "massimo.bourquin@epfl.ch"  # Always tell NCBI who you are
Entrez.api_key = "1b37cec2fdc64d05c5cc5337b17b09b4c609"

significant_kegg = pd.read_csv('../Data/KEGG_sign_sample.tsv', sep='\t', header = None)
significant_kegg.rename(columns = {0:'Sample'}, inplace = True)
significant_kegg.rename(columns = {1:'KEGG'}, inplace = True)

# if you want to parse taxonomies with python's NCBI entrez
def GetTaxonomy(ids_list):
    out_list = []
    for tax_id in ids_list:
        if tax_id != '0':
            handle = Entrez.efetch(db="Taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            if len(records) > 0:
                for i in records:
                    if i['Rank'] in ['genus','species','subspecies','varietas','forma']:
                        taxonomy = i["Lineage"]
                        out_list.append(taxonomy)
            else:
                continue
    return(out_list)

significant_kegg['Taxon name'] = ''
significant_kegg['Taxid'] = ''
significant_kegg['Contigs'] = ''

samples_list = list(set(significant_kegg['Sample'].values))
kegg_list = list(set(significant_kegg['KEGG'].values))
for sample in samples_list:
    sample_file = pd.read_csv('Data/Raw_func/' + sample + '_KEGG_counts.tsv', sep='\t')

    tax_file = pd.read_csv('Data/Raw_func/' + sample + '.labels.txt', sep='\t', header = None)

    sample_file_sign_sub = sample_file[sample_file['Geneid'].isin(kegg_list)]
    kegg_sub_list = list(set(list(sample_file_sign_sub['Geneid'].values)))
    for kegg in kegg_sub_list:
        contigs_list = sample_file_sign_sub.loc[sample_file_sign_sub['Geneid'] == kegg,'Chr'].to_list()
        contigs_list = [i for l in contigs_list for i in l.split(';')]
        significant_kegg.loc[(significant_kegg['KEGG'] == kegg) & (significant_kegg['Sample'] == sample), 'Contigs'] = ','.join(contigs_list)

        contigs_tax_name_l = tax_file[tax_file[1].isin(contigs_list)]
        contigs_tax_name = contigs_tax_name_l[2].to_list()
        significant_kegg.loc[(significant_kegg['KEGG'] == kegg) & (significant_kegg['Sample'] == sample), 'Taxon name'] = ','.join(contigs_tax_name)

        contigs_tax_ids = [i.split('(taxid ')[-1].split(')')[0] for i in contigs_tax_name if i != 'unclassified (taxid 0)']
        significant_kegg.loc[(significant_kegg['KEGG'] == kegg) & (significant_kegg['Sample'] == sample), 'Taxid'] = ','.join(contigs_tax_ids)



significant_kegg.to_csv('../Data/KEGG_sign_tax.tsv', sep='\t')

from Bio import SeqIO
import argparse
import pandas as pd
import random
import numpy as np
import itertools
import glob

def GetKmerSpectrum(seq, k):
    counts_dict = dict()
    for char in range(0,len(seq)-k+1):
        if seq[char,char+k] in counts_dict.keys():
            counts_dict[seq[char,char+k]] += 1
        else:
            counts_dict[seq[char,char+k]] = 1
    return(counts_dict)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--k_size',      type=int, help="k size of the kmers.")

    args = parser.parse_args()
    k = args.k_size

    kmers = GetAllKmers(k)
    counts_dict = pd.DataFrame(columns=['Group'] + kmers)
    count = 0

    fasta_files = glob.glob('*/*.fasta')
    for file in fasta_files:
        with open(file, "rt") as handle:
            for read in list(SeqIO.parse(handle, "fasta")):
                print(count)
                read_group = read.description.split('_MiSeq')[0]
                read_dict = GetKmerSpectrum(read.seq, k)
                read_dict['Gene'] = read.description
                read_dict['KEGG'] = str(file).split('/')[0]
                counts_dict = counts_dict.append(read_dict, ignore_index=True)
        count += 1

    counts_dict.to_csv('EUCI_kmer_counts_k' + str(k) + ''  + '.csv', sep=',')

if __name__ == "__main__":
    main()

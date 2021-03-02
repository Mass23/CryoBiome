import subprocess
import glob
from multiprocessing import Pool


def mafft(fasta_file):
    mafft_args = ['ginsi' + fasta_file + '>' + fasta_file.replace('.fasta','_aln.fasta')]
    subprocess.call(' '.join(mafft_args), shell=True, stdout=None, stderr=None)

def trimal(fasta_file):
    trimal_args = ['trimal', '-in', fasta_file.replace('.fasta','_aln.fasta'), '-out', fasta_file.replace('.fasta','_trim.fasta'), '-automated1']
    subprocess.call(' '.join(trimal_args), shell=True, stdout=None, stderr=None)

def iqtree(fasta_file):
    iqtree_args = ['iqtree', '-s', fasta_file.replace('.fasta','_trim.fasta'), '-m', 'GTR']
    subprocess.call(' '.join(iqtree_args), shell=True, stdout=None, stderr=None)

def process_file(input_file):
    print('File: ' + input_file)
    try:
        mafft(input_file)
        trimal(input_file)
        iqtree(input_file)
        print(input_file + ' processed!')
        return()
    except:
        print('Error with: ' + input_file)
        return()

def main():
    genus_fasta_files = glob.glob('Data/Genera/*.fasta')
    for genus_file in genus_fasta_files:
        process_file(genus_file)

if __name__ == '__main__':
    main()

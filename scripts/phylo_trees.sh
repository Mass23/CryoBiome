mafft --thread 8 PP1_005_ASV_seqs.fasta > PP1_005_aln_seqs.fasta
mafft --thread 8 PP2_005_ASV_seqs.fasta > PP2_005_aln_seqs.fasta

trimal -gt 0.95 -in PP1_005_aln_seqs.fasta -out PP1_005_trim_seqs.fasta
trimal -gt 0.95 -in PP2_005_aln_seqs.fasta -out PP2_005_trim_seqs.fasta

iqtree -s PP1_005_trim_seqs.fasta -fast -m GTR -T AUTO
iqtree -s PP2_005_trim_seqs.fasta -fast -m GTR -T AUTO
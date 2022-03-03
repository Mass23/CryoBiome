# CryoBiome
- All environments, scripts and jupyter notebooks required to replicate the Ecosystems Under Cryospheric Influence analyses are here within

## Overview
- Below are the individual folders, scripts and `Snakemake` files required to accomplish all the analyses in the manuscript
### Folders
- `0_preprocessing`: All scripts required for preprocessing both amplicon and metagenomic (metag) data
- `envs`: conda environments required for all analyses using `Snakemake`
- `scripts`: custom and exisiting scripts for individual steps and certain preprocessing steps including utility and accessory scripts for downstream analyses

### Steps
#### Amplicon analyses
- Scripts and jupyter notebooks involved in the amplicon analyses include the following
```
  - 0_preprocessing/PP1_classifer.sh
  - 0_preprocessing/PP2_classifer.sh
  - 0_preprocessing.R
  - 0_Dataset_exploration.R
  - 1_1_Logistice_classifiers.ipynb
  - 1_2_ASV_stats.R
  - 1_3_ASV_coefficients.R
  - 1_4_Phylogenetics.R
  - 2_1_Diff_abund.R
  - 2_2_refseq_cryobiome.R
  - 2_3_Core_microbiome.R
```

#### Metagenomic (metaG) analyses
- Scripts and `Snakemake` files required for metagenomic preprocessing and downstream analyses include the following:
```
  - 0_preprocessing/IMP_config.yaml: example config file for running IMP analyses
  - 0_preprocessing/metag_config.yaml: config file for preprocessing samples using IMP
  - 0_preprocessing/slurm.yaml: SLURM configuration file for running the analyses on a HPC-environment
  - 0_preprocessing/snakemake_parsable_launcher.sh: SBATCH launcher for running initial IMP preprocessing
  - 3_1_Enrichment.R
  - 3_2_Taxonomy.R
  - 3.3_Unknown_space.R
  - cluster_list: list of all clusters of KEGG+Unassigned genes based on "mmseqs2" analyses
  - config.yaml: config file for running the initial gene fasta extraction and "mmseqs2" clustering
  - coverage_Snakefile: to estimate coverage for each contig/gene in the assemblies
  - coverage_config.yaml: config file required to run "coverage_Snakefile"
  - extraction_Snakefile: to extract fasta sequences based on cluster IDs, get cluster stats, MAFFT alignment, consensus sequence selection and tRNA analyses.
  - extraction_config.yaml: config file required to run "extraction_Snakefile"
  - get_clusterID.sh: bash script to collect all contig IDs belonging to individual cluster IDs
  - kegg_plus_unassigned_list: list of all unassigned genes, and KEGG IDs
  - list: test list for testing Snakefiles
  - sample_list: list of samples used for metaG analyses
  - unassigned_clusters_EUCI.tsv: Tab-delimited file indicating which contig IDs belong to which respective clusters.
```

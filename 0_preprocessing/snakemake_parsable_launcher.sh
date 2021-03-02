#!/bin/bash -l

##############################
# SLURM
# NOTE: used for this script only, NOT for the snakemake call below

#SBATCH -J nomis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=5-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

# conda env name
IMP_ENV="snakemake"
# number of cores for snakemake
IMP_CORES=20
# number of jobs
IMP_JOBS=20
# IMP config file
IMP_CONFIG="config.yaml" # USER INPUT REQUIRED
# slurm config file
IMP_SLURM="slurm.yaml"
# slurm cluster call
IMP_CLUSTER="{cluster.call} {cluster.runtime} {cluster.threads}{threads} {cluster.partition} {cluster.nodes} {cluster.quality} --job-name={cluster.job-name}"
# IMP_CLUSTER="{cluster.call} {cluster.runtime} {cluster.mem_per_cpu} {cluster.threads}{threads} {cluster.partition} {cluster.nodes} {cluster.quality} --job-name={cluster.job-name}"

##############################
# IMP

# activate the env
conda activate ${IMP_ENV}
export PYTHONNOUSERSITE=TRUE

# run the pipeline
snakemake -s snakefile -rp --jobs ${IMP_JOBS} --local-cores 1 --configfile ${IMP_CONFIG} \
--use-conda --conda-prefix ${CONDA_PREFIX}/IMP --cluster-config ${IMP_SLURM} --cluster "${IMP_CLUSTER}" 

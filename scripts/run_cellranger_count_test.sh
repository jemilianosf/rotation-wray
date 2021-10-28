#!/bin/bash
#SBATCH -J cellranger_count # A single job name for the array
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH --mem 64GB

# Move to working dir
cd $SLURM_SUBMIT_DIR

# Parse sample info 

# Run cellranger_count
cellranger count --id=run_count_1kpbmcs \
--fastqs=/home/jes157/data_storage/cellranger_count_test/pbmc_1k_v3_fastqs \
--sample=pbmc_1k_v3 \
--transcriptome=/home/jes157/data_storage/cellranger_refs/refdata-gex-GRCh38-2020-A

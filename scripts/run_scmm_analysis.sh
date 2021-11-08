#!/bin/bash
#SBATCH -J scmm # A single job name for the array
#SBATCH -n 4 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH --mem 64GB
#SBATCH --array=1-4%2 

# Move to working dir
cd $SLURM_SUBMIT_DIR

module load Anaconda3
conda activate seurat

# Parse file info
file_name=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' /home/jes157/data_storage/data_seurat/rds_list.txt);

# Run analysis
Rscript scmm_analysis.R $file_name

conda deactivate
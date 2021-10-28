#!/bin/bash
#SBATCH -J cellranger_count # A single job name for the array
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH --mem 64GB
#SBATCH --array=1-4%2 

# Move to working dir
cd $SLURM_SUBMIT_DIR

# Parse sample info
sample_name=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' /home/jes157/data_storage/execution_files/sample_id_list.txt);
sample_id=run_count_${sample_name};

# Run cellranger_count
cellranger count --id=$sample_id \
--fastqs=/home/jes157/data_storage/HYTVFDMXX \
--sample=$sample_name \
--transcriptome=/home/jes157/data_storage/cellranger_refs/refdata-gex-mm10-2020-A

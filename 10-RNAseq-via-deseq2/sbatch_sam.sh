#!/bin/bash

#SBATCH --job-name=sbatch_sam
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --export=ALL
#SBATCH --array={barcode_IDs}

sample_id="$SLURM_ARRAY_TASK_ID"
absolute_path="path/to"
datafile=$(printf "barcode%02d" "${sample_id}")
env_path="${absolute_path}/miniforge3/envs/minimap2_env"

source "${absolute_path}/miniforge3/etc/profile.d/conda.sh"
conda activate "$env_path"

srun minimap2 \
    -t 10 \
    -ax splice \
    "${absolute_path}/{ref.fna}}" \
    "${absolute_path}/${datafile}.fastq.gz" > "${absolute_path}/${datafile}.sam"
#!/bin/bash

#SBATCH --job-name=sbatch_index
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
env_path="${absolute_path}/miniforge3/envs/samtools_env"

source "${absolute_path}/miniforge3/etc/profile.d/conda.sh"
conda activate "$env_path"

srun samtools index -b -@ 4 "${absolute_path}/${datafile}.sorted.bam" 
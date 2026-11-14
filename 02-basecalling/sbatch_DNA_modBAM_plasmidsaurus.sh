#!/bin/bash

#SBATCH --job-name=dorado
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=3
#SBATCH --mem-per-gpu=20G
#SBATCH --cpus-per-gpu=4
#SBATCH --time=14:00:00
#SBATCH --export=ALL

srun path/to/dorado-1.4.0-linux-x64/bin/dorado basecaller \
     path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
     path/to/pod5_file.pod5 \
     --modified-bases-models path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1,path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1 \
     --modified-bases-batchsize 128 \
     --min-qscore 9 \
     --batchsize 128 \
     > path/to/output.bam
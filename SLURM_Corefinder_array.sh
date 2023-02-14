#!/bin/bash 
#SBATCH --job-name=Corefinder
#SBATCH --ntasks=1
#SBATCH -c 64
#SBATCH -o /beegfs/projects/p450/out_files/Corefinder_%A_%a.out  
#SBATCH -t 1-00:00  
#SBATCH -p batch 
#SBATCH --array=1-80
#SBATCH -e /beegfs/projects/p450/error_filesCorefinder_error_%A_%a.txt



source ~/.bashrc
conda activate /beegfs/home/fbiermann/miniconda3_supernew/envs/frida


for file in /projects/p450/tryptorubin_cores/listpositives_threshold_01_split/"{$SLURM_ARRAY_TASK_ID}"_listpositiveNCBI_combined.fasta ; do python3 /projects/p450/tryptorubin_cores/Corefinder.py -i $file -o /projects/p450/tryptorubin_cores/corefinder_output_"{$SLURM_ARRAY_TASK_ID}"/ -e f.biermann@bio.uni-frankfurt.de

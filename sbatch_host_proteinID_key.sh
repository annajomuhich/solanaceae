#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/host_proteinID_key_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/host_proteinID_key_stderr-%j.txt
#SBATCH -J host_proteinID_key
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
INPUT_FILE=$1
COLUMN_NAME=$2
OUTPUT_FILE=$3

Rscript ~/solanaceae/host_proteinID_key.R \
"$INPUT_FILE" \
"$COLUMN_NAME" \
"$OUTPUT_FILE"

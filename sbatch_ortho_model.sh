#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/ortho_model_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/ortho_model_stderr-%j.txt
#SBATCH -J ortho_model
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
COUNTS_FILE=$1
OUTPUT_DIR=$2


Rscript ~/solanaceae/ortho_model.R \
"$COUNTS_FILE" \
"$OUTPUT_DIR"

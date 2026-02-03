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
COUNTS_FILE1=$1
COUNTS_FILE2=$2
ORTHO_FILE=$3
OUTPUT_DIR=$4


Rscript ~/solanaceae/ortho_model.R \
"$COUNTS_FILE1" "$COUNTS_FILE2" \
"$ORTHO_FILE" \
"$OUTPUT_DIR"

#!/bin/bash
#SBATCH -o /home/ajmuhich/slurm-log/blastp_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/blastp_stderr-%j.txt
#SBATCH -J blastp
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load conda
conda activate blast

if [ $# -lt 3 ]; then
  echo "Usage: sbatch sbatch_blastp_annot.sh <QUERY> <DB> <OUTPUT_DIR>"
  exit 1
fi

QUERY="$1"
DB="$2"
OUTPUT_CSV="$3"

RAW_OUT="raw_out_${SLURM_JOB_ID}.tsv"

blastp \
    -query $QUERY \
    -db $DB \
    -out "$RAW_OUT \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle"

module load R
Rscript ~/solanaceae/blastp_reformat.R $RAW_OUT $DB $OUTPUT_CSV

rm $RAW_OUT
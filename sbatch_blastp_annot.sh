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

mkdir -p $OUTPUT_DIR

blastp \
    -query $QUERY \
    -db $DB \
    -out raw_out.tsv \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle"

Rscript ~/solanaceae/blastp_reformat.R raw_out.tsv $DB $OUTPUT_CSV
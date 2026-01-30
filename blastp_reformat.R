### reformat blastp output
### January 2026 AJM

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
	stop("Usage: Rscript blastp_reformat.R <input_file> <db> <output_file>")
}

input_file <- args[1]
db <- args[2]
output_file <- args[3]

message("Loading input file ", input_file, "...")
df <- read.delim(input_file, header = F)
file.remove(input_file) #remove temporary input

colnames(df) <- c("qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore", "stitle")

db <- sub(".*/", "", db)
df <- df %>%
	mutate(description = sub("^[^ ]* ", "", stitle)) %>%
	select(qseqid, sseqid, pident, evalue, description)

names(df) <- ifelse(names(df) == "qseqid",
										"qseqid",
										paste0(db, "_", names(df)))

df %>% write.csv(output_file, row.names = F)

### Converting Protein IDs (orthofinder output) to Gene IDs (LOC IDs)
### Jan 2026 AJM

### --------------------- Parse arguments ------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
	stop("Usage: Rscript xp_to_loc.R <input_file> <column_name> <output_file>")
}

input_file <- args[1]
column_name <- args[2]
output_file <- args[3]

### --------------------- Load libraries ----------------

library(rentrez)
library(tidyverse)
library(readr)

### ------------------ Prepare input data --------------

df <- read_delim(input_file, stringsAsFactors = FALSE)

if (!column_name %in% colnames(df)) {
	stop("Column '", column_name, "' not found in input file")
}

xp_vec <- unique(
	trimws(
		unlist(
			strsplit(df[[column_name]], ",")
		)
	)
)

xp_vec <- xp_vec[1:10] #subset for testing

### ----------------- Define function ----------------------------

#define function to convert XP id to LOC ID
xp_to_loc_df <- function(xp_ids) {
	
	n <- length(xp_ids)
	
	data.frame(
		protein_id = xp_ids,
		loc_id = vapply(
			seq_along(xp_ids),
			function(i) {
				prot <- xp_ids[i]
				
				message("Getting LOC ID for ", prot, " (", i, " / ", n, ")")
				
				gene_link <- tryCatch(
					entrez_link(
						dbfrom = "protein",
						id = prot,
						db = "gene"
					)$links$protein_gene,
					error = function(e) NULL
				)
				
				if (is.null(gene_link) || length(gene_link) == 0)
					return(NA_character_)
				
				gene_sum <- entrez_summary(
					db = "gene",
					id = gene_link[1]
				)[[1]]
				
				loc <- gene_sum[[1]]
				
				if (!is.na(loc) && !grepl("^LOC", loc))
					loc <- paste0("LOC", loc)
				
				loc
			},
			FUN.VALUE = character(1)
		),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

### --------- Run and write output -------------------

loc_df <- xp_to_loc_df(xp_vec)

#create directory if it doesn't exist
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) {
	dir.create(out_dir, recursive = TRUE)
}

#write
loc_df %>% write.csv(output_file, row.names = F)




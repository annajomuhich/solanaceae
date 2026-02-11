### plotting mean-variance for each organism
### February 2026 AJM

library(tidyverse)

### ------------------------ Host --------------------------------

host1_counts_file <- "data/norm_counts/Pepper_Host_expression.csv"
host2_counts_file <- "data/norm_counts/Tomato_Host_expression.csv"
ortho_file <- "data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_1to1.csv"

#load count file for each host
message("Loading host1 counts file: ", host1_counts_file)
host1 <- read.csv(host1_counts_file) %>%
	dplyr::select(iso_name, tray, seq_batch, gene, cpm)
message("Loading host2 counts file: ", host2_counts_file)
host2 <- read.csv(host2_counts_file) %>%
	dplyr::select(iso_name, tray, seq_batch, gene, cpm)
message("Loading ortholog file: ", ortho_file)
ortho <- read.csv(ortho_file) %>%
	dplyr::select(!Orthogroup) %>%
	mutate(ortholog = paste0(Ca, "_", Sl))

#Join orthologs with data
host1 <- left_join(host1, ortho, join_by(gene == Ca)) %>%
	filter(!is.na(ortholog)) %>% #remove counts for non-orthologs
	dplyr::select(!Sl) %>%
	mutate(genotype = "pepper")
host2 <- left_join(host2, ortho, join_by(gene == Sl)) %>%
	filter(!is.na(ortholog)) %>%
	dplyr::select(!Ca) %>%
	mutate(genotype = "tomato")

colnames(host1) == colnames(host2)
df_long <- rbind(host1, host2) %>%
	dplyr::select(!gene)

#check for NA values in count data (needs to be FALSE)
any(is.na(df_long$cpm))

df_long <- df_long %>%
	mutate(sample_ID = paste(genotype, iso_name, tray, seq_batch, sep = "_"))

df_long <- df_long %>%
	dplyr::select(sample_ID, ortholog, cpm)

#keep only genes that have >=1 CPM in at least 20% of the samples.
removed_genes <- df_long %>%
	group_by(ortholog) %>%
	summarise(prop_expressed = mean(cpm >= 1), .groups = "drop") %>%
	filter(prop_expressed < 0.20) %>%
	pull(ortholog)

df_long <- df_long %>%
	filter(!ortholog %in% removed_genes)

#genes as rows, samples as columns
host_mat <- df_long %>%
	pivot_wider(names_from = sample_ID,
							values_from = cpm) %>%
	column_to_rownames("ortholog")


#Check for NAs 
anyNA(host_mat)
rowSums(is.na(host_mat))
#flag NAs
na_rows <- rowSums(is.na(host_mat)) > 0
#get list of genes to remove
removed_rows <- row.names(host_mat)[na_rows]
#remove genes
message(
	"Removing the following ", length(removed_rows),
	" genes (not present when infecting both hosts):\n",
	paste(removed_rows, collapse = "\n")
)
host_mat <- host_mat[!na_rows, ]
anyNA(host_mat)

host_summary <- data.frame(mean = rowMeans(host_mat),
													 var = apply(host_mat, 1, var))
host_summary %>%
	ggplot(aes(x = mean, y = var)) +
	geom_point()

### --------------------- Bcin -----------------------

host1_counts_file <- "data/norm_counts/Pepper_B.cinerea_expression_long.csv"
host2_counts_file <- "data/norm_counts/Tomato_B.cinerea_expression_long.csv"

#load count file for each host
message("Loading host1 counts file: ", host1_counts_file)
host1 <- read.csv(host1_counts_file)
message("Loading host2 counts file: ", host2_counts_file)
host2 <- read.csv(host2_counts_file)

host1 <- host1 %>% dplyr::select(genotype, iso_name, tray, seq_batch, gene, CPM)
host2 <- host2 %>% dplyr::select(genotype, iso_name, tray, seq_batch, gene, CPM)

df_long <- rbind(host1, host2)

df_long <- df_long %>%
	mutate(sample_ID = paste(genotype, iso_name, tray, seq_batch, sep = "_"))

df_long <- df_long %>%
	dplyr::select(sample_ID, gene, CPM)

#keep only genes that have >=1 CPM in at least 20% of the samples.
removed_genes <- df_long %>%
	group_by(gene) %>%
	summarise(prop_expressed = mean(CPM >= 1), .groups = "drop") %>%
	filter(prop_expressed < 0.20) %>%
	pull(gene)

#genes as rows, samples as columns
bcin_mat <- df_long %>%
	pivot_wider(names_from = sample_ID,
							values_from = CPM) %>%
	column_to_rownames("gene")

#Check for NAs 
anyNA(bcin_mat)
rowSums(is.na(bcin_mat))
#flag NAs
na_rows <- rowSums(is.na(bcin_mat)) > 0
#get list of genes to remove
removed_rows <- row.names(bcin_mat)[na_rows]
#remove genes
message(
	"Removing the following ", length(removed_rows),
	" genes (not present when infecting both bcins):\n",
	paste(removed_rows, collapse = "\n")
)
bcin_mat <- bcin_mat[!na_rows, ]
anyNA(bcin_mat)

bcin_summary <- data.frame(mean = rowMeans(bcin_mat),
													 var = apply(bcin_mat, 1, var))

### --------- Plot ----------------------------------------

host_summary %>%
	ggplot(aes(x = mean, y = var)) +
	geom_point() +
	geom_smooth(method = "loess") +
	theme(xlim(0, 8000),
				ylim(0, 5e+08)) +
	theme_minimal() +
	ylab("variance") +
	xlab("mean (cpm)")
ggsave("figures/transcript_mean-var/ortho_mean-var.png", height=5, width=5)

bcin_summary %>%
	filter(mean < 50000) %>%
	filter(var < 5.0e+08) %>%
	ggplot(aes(x = mean, y = var)) +
	geom_point() +
	geom_smooth(method = "loess") +
	theme(xlim(0, 8000),
				ylim(0, 5e+08)) +
	theme_minimal()+
	ylab("variance") +
	xlab("mean (cpm)")
ggsave("figures/transcript_mean-var/bcin_mean-var.png", height=5, width=5)

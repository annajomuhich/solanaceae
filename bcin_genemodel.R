##### Bcin gene model - Solanaceae
##### For comparison of Bcin gene expression across 2 hosts
##### January 2026 AJM

### Define paths from arguments =====================
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
host1_counts_file <- args[1]		# Counts file (host_norm_counts_expressed.csv)
host2_counts_file <- args[2]   # Counts file for other host (host_norm_counts_expressed.csv)
output_dir <- args[3]     # Output directory

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

### Load packages and data ===================
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)

#load count file for each host
message("Loading host1 counts file: ", host1_counts_file)
host1 <- read.csv(host1_counts_file)
message("Loading host2 counts file: ", host2_counts_file)
host2 <- read.csv(host2_counts_file)

### Reformat count data ======================

host1 <- host1 %>% dplyr::select(genotype, iso_name, tray, seq_batch, gene, CPM)
host2 <- host2 %>% dplyr::select(genotype, iso_name, tray, seq_batch, gene, CPM)

df <- rbind(host1, host2)

df <- df %>%
	pivot_wider(names_from = gene,
							values_from = CPM)

#convert categorical variables to factor
head(df)
df$genotype <- as.factor(df$genotype)
df$iso_name <- as.factor(df$iso_name)
df$tray <- as.factor(df$tray)
df$seq_batch <- as.factor(df$seq_batch)
head(df) 

### Remove genes not expressed in both hosts ======================

#Check for NAs 
anyNA(df)
colSums(is.na(df))
#flag NAs
na_cols <- colSums(is.na(df)) > 0
#get list of genes to remove
removed_cols <- names(df)[na_cols]
#remove genes
message(
	"Removing the following ", length(removed_cols),
	" genes (not present when infecting both hosts):\n",
	paste(removed_cols, collapse = "\n")
)
df <- df[, !na_cols]

### Define model function ======================

analyze_gene <- function(gene, df) {
	
	message("modeling ", gene)
	
	# model
	formula <- as.formula(paste(
		gene,
		"~ genotype + iso_name + genotype * iso_name + (1|tray) + (1|seq_batch)"
	))
	
	model <- glmmTMB(formula, data = df, family = nbinom2) %>%
		suppressMessages()
	
	## ---- ANOVA + variance ----
	anova <- car::Anova(model) %>%
		as.data.frame() %>%
		rownames_to_column("variable") %>%
		mutate(gene = gene)
	
	fixed_var <- diag(vcov(model)$cond) %>%
		as.data.frame() %>%
		rownames_to_column("term") %>%
		setNames(c("term", "variance"))
	
	var_sums <- fixed_var %>%
		summarise(
			genotype = sum(variance[grepl("^genotype", term) & !grepl(":", term)], na.rm = TRUE),
			iso_name = sum(variance[grepl("^iso_name", term)], na.rm = TRUE),
			`genotype:iso_name` = sum(variance[grepl("^genotype", term) & grepl(":", term)], na.rm = TRUE),
			intercept = sum(variance[grepl("Intercept", term)], na.rm = TRUE),
			tot_var = sum(genotype, iso_name, `genotype:iso_name`, intercept)
		) %>%
		mutate(across(-tot_var, ~ .x / tot_var)) %>%
		select(-tot_var) %>%
		pivot_longer(everything(), names_to = "variable", values_to = "variance")
	
	anova <- full_join(anova, var_sums, by = "variable") %>%
		select(gene, everything())
	
	## ---- EMMs ----
	emm_sum <- emmeans(model, ~ iso_name + genotype) %>%
		summary() %>%
		as.data.frame() %>%
		select(iso_name, genotype, emmean, SE) %>%
		mutate(gene = gene)
	
	## ---- DEG ----
	DEG <- emmeans(model, specs = "genotype") %>%
		contrast(method = "revpairwise") %>%
		summary() %>%
		mutate(
			log2FC = estimate / log(2),
			SE = SE / log(2),
			gene = gene
		) %>%
		select(gene, everything(), -estimate)
	
	list(
		anova = anova,
		emm = emm_sum,
		DEG = DEG
	)
}


### Run function across genes ===================================

#get list of genes to run
df_long <- df %>%
	pivot_longer(cols = starts_with("Bcin"),
							 names_to = "gene",
							 values_to = "CPM")
genes <- unique(df_long$gene)

#setup outputs
results <- list()
failed_genes <- character()

#run
for (gene in genes) {
	tryCatch(
		results[[gene]] <- analyze_gene(gene, df),
		error = function(e) {
			message("Error for gene ", gene, " — skipping")
			failed_genes <<- c(failed_genes, gene)
		}
	)
}

### Reformat and write outputs ===========================

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_all   <- bind_rows(lapply(results, `[[`, "emm"))
DEG_all   <- bind_rows(lapply(results, `[[`, "DEG"))

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
  p_values <- x$`Pr(>Chisq)`
  x$p_adj <- p.adjust(p_values, method = "BH")
  return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#Replace DEG p value with anova genotype p value
#remove DEG p value
DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "genotype") %>%
	dplyr::select(gene, p_adj)
#join with host anova
DEG_all <- left_join(DEG_all, anova_ps, by = "gene")

#Put failed genes in a df
failed_genes_df <- data.frame(failed_gene = failed_genes, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", output_dir)
dir.create(output_dir, recursive = T)
write.csv(emm_all, paste0(output_dir, "bcin_adjusted_emmeans.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "bcin_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "bcin_DEGs.csv"), row.names = F)
write.csv(failed_genes_df, paste0(output_dir, "failed_genes.csv"), row.names = FALSE)

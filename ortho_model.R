### Modeling all 1:1 orthologs
### February 2026 AJM

### Define paths from arguments =====================
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
ortho_counts_file <- args[1]   # Ortholog counts
output_dir <- args[2]     # Output directory

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

### Load packages and data ====================================

#install.packages('partR2')

library(tidyverse)
library(emmeans)
library(car)
#library(purrr)
#library(lme4)
library(glmmTMB)

#load count file for each host
message("Loading ortholog counts file: ", ortho_counts_file)
df_long <- read.csv(ortho_counts_file)

df <- df_long %>%
	pivot_wider(names_from = ortholog,
							values_from = cpm)
	
#convert categorical variables to factor
head(df)
df$genotype <- as.factor(df$genotype)
df$iso_name <- as.factor(df$iso_name)
df$tray <- as.factor(df$tray)
df$seq_batch <- as.factor(df$seq_batch)
df$infected <- as.factor(df$infected)
head(df) 

### Define model function ======================

analyze_gene <- function(gene, df) {
	
	message("modeling ", gene)
	
	# model
	formula <- as.formula(paste(
		gene,
		"~ genotype * iso_name +
       (1 | tray) + (1 | seq_batch)"
	))

	model <- glmmTMB(formula, data = df, family = gaussian) %>%
		suppressMessages()
	
	#collect convergence note
	conv_note <- model$fit$message
	
	## ---- ANOVA + variance ----
	anova <- car::Anova(model) %>%
		as.data.frame() %>%
		rownames_to_column("variable") %>%
		mutate(gene = gene) %>%
		mutate(convergence_note = conv_note)
	
	fixed_var <- diag(vcov(model)$cond) %>%
		as.data.frame() %>%
		rownames_to_column("term") %>%
		setNames(c("term", "variance"))
	
	var_sums <- fixed_var %>%
				summarise(
					genotype = sum(variance[grepl("^genotype", term) & !grepl(":", term)], na.rm = TRUE),
					#infected = sum(variance[grepl("^infected", term) & !grepl(":", term)], na.rm = TRUE),
					#`genotype:infected` = sum(variance[grepl("^genotype", term) & grepl(":", term)], na.rm = TRUE),
					#`infected:iso_name` = sum(variance[grepl("^infected", term) & grepl(":", term)], na.rm = TRUE),
					iso_name = sum(variance[grepl("^iso_name", term) & !grepl(":", term)], na.rm = TRUE),
					`genotype:iso_name` = sum(variance[grepl("^genotype", term) & grepl(":", term)], na.rm = TRUE),
		 			intercept = sum(variance[grepl("Intercept", term)], na.rm = TRUE),
					tot_var = sum(genotype, iso_name, `genotype:iso_name`, intercept)
				) %>%
		mutate(across(-tot_var, ~ .x / tot_var)) %>%
		dplyr::select(-tot_var) %>%
		pivot_longer(everything(), names_to = "variable", values_to = "variance")
	
	anova <- full_join(anova, var_sums, by = "variable") %>%
		dplyr::select(gene, everything())
	anova$gene <- gene
	
	## ---- EMMs ----
	emm_log <- emmeans(model, ~ iso_name + genotype) %>%
		summary() %>%
		as.data.frame() %>%
		#dplyr::select(iso_name, genotype, emmean, SE) %>%
		dplyr::rename(emmean_log = emmean) %>%
		mutate(gene = gene)
	emm_resp <- emmeans(model, ~ iso_name + genotype, type = "response") %>%
		summary() %>%
		as.data.frame() %>%
		#dplyr::select(iso_name, genotype, response, SE) %>%
		dplyr::rename(emmean_response = emmean) %>%
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
		dplyr::select(gene, everything(), -estimate)
	
	list(
		anova = anova,
		emm_log = emm_log,
		emm_resp = emm_resp,
		DEG = DEG
	)
}


### Run function across genes ===================================

#get list of genes to run
genes <- unique(df_long$ortholog)

#subset for testing
#genes <- genes[1100:1110]

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

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_log_all   <- bind_rows(lapply(results, `[[`, "emm_log"))
emm_resp_all   <- bind_rows(lapply(results, `[[`, "emm_resp"))
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

message("Writing results to output directory: ", output_dir)
dir.create(output_dir, recursive = T)
write.csv(emm_log_all, paste0(output_dir, "ortho_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "ortho_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "ortho_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "ortho_DEGs.csv"), row.names = F)
write.csv(failed_genes_df, paste0(output_dir, "failed_genes.csv"), row.names = FALSE)


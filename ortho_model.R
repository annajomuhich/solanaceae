### Modeling all 1:1 orthologs
### October 2025 AJM

#install.packages('partR2')

library(tidyverse)
library(emmeans)
library(car)
library(purrr)
library(lme4)

library(glmmTMB) #Had trouble pulling out variances properly with this package


### ---------------- Get Ca and Sl expression for all 1:1 orthologs ---------------------------------
#read in input files
cadf <- read.csv(file = "data/norm_counts/Pepper_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv") %>%
	dplyr::select(genotype, iso_name, tray, seq_batch, gene, cpm)
sldf <- read.csv(file = "data/norm_counts/Tomato_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv") %>%
	dplyr::select(genotype, iso_name, tray, seq_batch, gene, cpm)

#read in 1:1 orthos
ortho <- read.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_1to1.csv") %>%
	select(!Orthogroup) %>%
	mutate(ortholog = paste0(Ca, "_", Sl)) %>%
	distinct()

##THIS IS CREATING DUPLICATES SOMEHOW
df1 <- left_join(cadf, ortho, join_by(gene == Ca)) %>%
	filter(!is.na(ortholog)) %>% #remove counts for non-orthologs
	select(!Sl)
df2 <- left_join(sldf, ortho, join_by(gene == Sl)) %>%
	filter(!is.na(ortholog)) %>%
	select(!Ca)

colnames(df1) == colnames(df2)
df_long <- rbind(df1, df2) %>%
	select(!gene)

#check for NA values in count data (needs to be FALSE)
any(is.na(df_long$cpm))

#pivot wider
df <- df_long %>%
	pivot_wider(names_from = ortholog,
							values_from = cpm)

### ISSUE NOW IS THAT MY ORTHOLOG IDS SOME HAVE NA. Need to fix this

#convert categorical variables to factor
head(df)
df_long$genotype <- as.factor(df_long$genotype)
df_long$iso_name <- as.factor(df_long$iso_name)
df_long$infected <- as.factor(df_long$infected)
head(df_long)

### ---------------- Run GLMM on 1:1s ---------------------------------

#Want to check and see if we get similar results for p values with glmmTMB()
#so we can justify using lm() for the variance

### Define function to run anovas and get variance for all genes

run_anovas_glmm <- function(df_long, genes) {
	map_dfr(genes, function(gene) {
		message("\nModeling ", gene, "\n")
		# Filter data for this gene
		df_filt <- df_long %>% filter(ortholog == gene)
		# Skip if too few data points
		if (nrow(df_filt) < 5) return(NULL)
		# Fit model (suppress warnings but not errors)
		model <- tryCatch(
			suppressMessages(
				model <- glmmTMB(
					cpm ~ host_species + iso_name + host_species:iso_name,
					data = df_filt,
					family = nbinom2()))
			,
			error = function(e) {
				message("❌ Model failed for ", gene, ": ", e$message)
				return(NULL)
			}
		)
		
		# Try to extract ANOVA-like table
		anova_tbl <- tryCatch({
			aov <- car::Anova(model)                     # may fail for some merMod objects
			aov <- tibble::rownames_to_column(as.data.frame(aov), var = "variable")
			aov$ortholog <- gene
			aov
		}, error = function(e) {
			message("⚠️ car::Anova failed for ", gene, ". (", e$message, ")")
		})
	})
}

### Run on all 1:1s
#get list of genes to iterate through
genes <- unique(df_long$ortholog)
#run (prob take 1-2h)
anova_all <- run_anovas_glmm(df_long, genes)

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

#write out
anova_corrected %>% write.csv("data/host_ortho/anova_glmm.csv", row.names = F)

### GO enrichment analysis on Botrytis genes
### January 2026 AJM

library(tidyverse)
library(topGO)
library(GO.db)

#define function to setup geneID2GO using the annotation file
geneID2GO_setup <- function(df_annot) {
	df_annot <- df_annot %>%
		dplyr::select(gene, GO) %>%
		distinct(gene, .keep_all = TRUE) %>%
		filter(GO != "")
	# Rename gene column
	#names(df_annot)[names(df_annot) == "locusName"] <- "gene"
	# Build geneID2GO list directly in memory
	geneID2GO <- setNames(
		strsplit(df_annot$GO, "\\s*,\\s*|\\s+"), #split string by comma or space
		df_annot$gene
	)
	return(geneID2GO)
}

#define function to run topGO directly comparing two groups
run_topGO <- function(groupA, groupB, topNodes = 100) {
	#get entire list of genes to compare
	all_genes <- unique(c(groupA, groupB))
	#set up binary vector
	geneList <- factor(as.integer(all_genes %in% groupA))
	names(geneList) <- all_genes
	#run topGO
	GOdata <- new("topGOdata",
								ontology = "MF",  # or "MF" or "CC"
								allGenes = geneList,
								geneSel = function(p) p == 1, #only include genes in hostspec
								annot = annFUN.gene2GO,
								gene2GO = geneID2GO,
								nodeSize = 10)
	resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	# Analysis of Results
	results <- GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes)
	# get non-truncated GO terms	
	full_terms <- Term(GOTERM[results$GO.ID])
	results$Term <- full_terms
	# Correct for multiple testing
	p_values <- results$classicFisher
	results$p_adj <- p.adjust(p_values, method = "BH")
	return(results)
}

### Prepare gene universe (geneID2GO) -----------------------------------------------
#load annotation file
df_annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv")
# Reformat annotation
df_annot <- df_annot %>%
	distinct(X.Gene.ID., .keep_all = TRUE) %>% # Removing duplicates based on 'gene'
	dplyr::select(X.Gene.ID.,X.Computed.GO.Function.IDs.,X.Computed.GO.Process.IDs.) # Selecting only the 'gene' and 'GO' columns
# Rename the X.Gene.ID. column to gene
names(df_annot)[names(df_annot) == "X.Gene.ID."] <- "gene"
# Combine the two columns into one "GO" column, removing "N/A"
df_annot <- df_annot %>%
	mutate(
		GO = paste(
			ifelse(X.Computed.GO.Function.IDs. != "N/A", X.Computed.GO.Function.IDs., ""),
			ifelse(X.Computed.GO.Process.IDs. != "N/A", X.Computed.GO.Process.IDs., ""),
			sep = ",") %>%
			gsub("(^,|,$|,,)", "", .)) %>%  # Remove leading/trailing commas and duplicate commas
	dplyr::select(gene, GO)  # Keep only the Gene ID and the new GO column
#setup geneID2GO
geneID2GO <- geneID2GO_setup(df_annot)
#get all genes
all_genes <- df_annot %>% pull(gene)

### Run GO enrichment on Ca-specific genes vs all Botrytis genes -------------------------------------------

#get Ca-specific genes
deg <- read.csv("data/bcin_expr/DEG/bcin_DEG.csv")
deg_Ca <- deg %>% filter(significance == "Pepper Upregulated") %>% pull(gene)

result <- run_topGO(deg_Ca, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
#result %>% write.csv("data/bcin_expr/DEG/CaUp_GO.csv", row.names = F)

### Run GO enrichment on Sl-specific genes vs all Botrytis genes -------------------------------------------

#get Sl-specific genes
deg_Sl <- deg %>% filter(significance == "Tomato Upregulated") %>% pull(gene)

result <- run_topGO(deg_Sl, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
#result %>% write.csv("data/bcin_expr/DEG/SlUp_GO.csv", row.names = F)

### Run GO enrichment on 972 host-specific genes vs all Botrytis genes ----------------------------------------

#get host-specific only genes
host_spec <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept") %>%
	pivot_wider(names_from = "variable",
							values_from = "p_adj") %>%
	filter(genotype < 0.05 & `genotype:iso_name` > 0.05 & iso_name > 0.05) %>%
	pull(gene)

result <- run_topGO(host_spec, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
#result %>% write.csv("data/bcin_expr/host_iso_specific/hostonlyUp_GO.csv", row.names = F)

### Run GO enrichment on 798 iso-specific genes vs all Botrytis genes -----------------------------------------

#get iso-specific only genes
iso_spec <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept") %>%
	pivot_wider(names_from = "variable",
							values_from = "p_adj") %>%
	filter(genotype > 0.05 & `genotype:iso_name` > 0.05 & iso_name < 0.05) %>%
	pull(gene)

result <- run_topGO(iso_spec, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
#result %>% write.csv("data/bcin_expr/host_iso_specific/isoonlyUp_GO.csv", row.names = F)

### Run GO enrichment on all iso-specific genes vs all Botrytis genes -------------------------

#get all iso-specific genes
iso_spec <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept") %>%
	pivot_wider(names_from = "variable",
							values_from = "p_adj") %>%
	filter(iso_name < 0.05) %>%
	pull(gene)

#run
result <- run_topGO(iso_spec, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
result %>% write.csv("data/bcin_expr/host_iso_specific/isoallUp_GO.csv", row.names = F)

### Run GO enrichment on all host-specific genes vs all Botrytis genes -------------------------

#get all host-specific genes
host_spec <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept") %>%
	pivot_wider(names_from = "variable",
							values_from = "p_adj") %>%
	filter(genotype < 0.05) %>%
	pull(gene)

#run
result <- run_topGO(host_spec, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
result %>% write.csv("data/bcin_expr/host_iso_specific/hostallUp_GO.csv", row.names = F)

### Run GO enrichment on all interaction-specific genes vs all Botrytis genes -------------------------

#get all host-specific genes
interaction_spec <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept") %>%
	pivot_wider(names_from = "variable",
							values_from = "p_adj") %>%
	filter(`genotype:iso_name` < 0.05) %>%
	pull(gene)

#run
result <- run_topGO(interaction_spec, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
result %>% write.csv("data/bcin_expr/host_iso_specific/interactionallUp_GO.csv", row.names = F)

### Bcin genes differentially expressed across hosts
### January 2026 AJM

library(tidyverse)
library(ggrepel)
library(topGO)

### =================== Combine data - DEG and anova ====================================

deg <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_DEGs.csv")

#join with annotations
annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv")
annot <- annot[,c("X.Gene.ID.", "X.Gene.Name.or.Symbol.", "X.PFam.Description.")]
names(annot)[names(annot) == "X.Gene.Name.or.Symbol."] <- "gene_name"
names(annot)[names(annot) == "X.PFam.Description."] <- "gene_desc"
annot <- annot[!duplicated(annot$X.Gene.ID.),] #keep first entry of each gene
deg <- left_join(deg, annot, join_by(gene == X.Gene.ID.))

# Add -log10(p-value) and significance column
deg <- deg %>%
	mutate(
		neg_log10_p = -log10(p_adj),
		significance = case_when(
			p_adj < 0.05 & log2FC > 1 ~ "Tomato Upregulated",
			p_adj < 0.05 & log2FC < -1 ~ "Pepper Upregulated",
			TRUE ~ "Not Significant"
		)
	)

deg %>% write.csv("data/bcin_expr/DEG/bcin_DEG.csv", row.names = F)

### ========================= Plot DEGs ==========================================

# Volcano plot - no labels
deg %>%
ggplot(aes(x = log2FC, y = neg_log10_p, color = significance)) +
	geom_point(size = 1) +
	scale_color_manual(values = c("Tomato Upregulated" = "red",
																"Pepper Upregulated" = "navy",
																"Not Significant" = "dark grey")) +
	# geom_text_repel(aes(label = label),
	# 								color = "black",
	# 								fontface = "bold",
	# 								max.overlaps = 40,
	# 								#box.padding = 0.7,
	# 								#point.padding = 0.5,
	# 								size = 2) +
	labs(x = "Log2 Fold Change", y = "-Log10 P-Value") +
	theme_minimal()
	#ggtitle("Differentially expressed genes in Common Bean \ninfected with Botrytis cinerea")

# #manually labeled in excel
# deg_label <- read.csv("data/bcin_RNAseq/bcin_DEG/unadj_pval/bcin_DEG_labeled.csv") #this one has the unadj p values, can just take the label off it
# deg_label_column <- deg_label %>% dplyr::select(gene, label)
# deg_label <- left_join(deg, deg_label_column, by = 'gene')

# #decide which genes to label
# deg_label_cutoff <- deg %>%
# 	filter(significance != "Not Significant") %>%
# 	filter(neg_log10_p > 50 |
# 				 	log2FC > 4 |
# 				 	log2FC < -7) %>%
# 	filter(gene_name != "N/A" | gene_desc != "N/A") %>%
# 	mutate(label = if_else(gene_name == "N/A", NA, gene_name))
# 
# deg_label_cutoff %>% write.csv("data/bcin_RNAseq/bcin_DEG/bcin_DEG_labels.csv", row.names = F)
# #manually added remaining labels in excel then read back in
# labels <- read.csv("data/bcin_RNAseq/bcin_DEG/bcin_DEG_labels.csv") %>%
# 	dplyr::select(gene, label)
# 
# deg_label <- left_join(deg, labels, by = "gene")
# 
# # Volcano plot - with labels
# deg_label %>%
# 	filter(neg_log10_p > 1.3) %>% #clean up lower outliers
# 	mutate(label = if_else(gene == "Bcin15g03990", "oxidase", label)) %>% #add label for the big cowpea outlier based on blastp
# ggplot(aes(x = log2FC, y = neg_log10_p, color = significance)) +
# 	geom_point(size = 1) +
# 	scale_color_manual(values = c("Common Bean Upregulated" = "#88cdebff",
# 																"Cowpea Upregulated" = "#dd3497ff",
# 																"Not Significant" = "dark grey")) +
# 	geom_text_repel(aes(label = label),
# 									color = "black",
# 									fontface = "bold",
# 									max.overlaps = 40,
# 									#box.padding = 0.7,
# 									#point.padding = 0.5,
# 									size = 2) +
# 	labs(x = "Log2 Fold Change", y = "-Log10 P-Value") +
# 	theme_minimal()
# 	#ggtitle("Differentially expressed genes in Common Bean \ninfected with Botrytis cinerea")
# 
# ggsave("figures/bcin_transcriptome/bcin_DEG.png", height = 4, width = 9)

### ================ Format supplementary table ===================================

deg_table <- deg %>%
	filter(significance != "Not Significant")
deg_table <- deg_table %>%
	dplyr::select(gene, gene_name, gene_desc, significance, log2FC, neg_log10_p)
colnames(deg_table) <- c("Gene ID", "Gene Name", "Gene Description",
												 "Differential Expression", "log2FC", "Negative Log10 P Value")
deg_table %>% write.csv("supplementary_tables/csvs/Supplementary_Table_BcinDEGs.csv", row.names = F)

#### ============ SM-related GO enrichment of host specific genes - 5/29/25 =================

#Performing GO enrichment of host-specific vs non-host specific genes against a SM-specific gene universe
#I actually need not a list of GO terms but a list of genes, i.e. a subset of the annotation file

#load SM_GO terms
sm_go <- read.csv("~/UCDavis/Klieb_Lab/Projects/RNASeq/full_pipeline/remove_7isos/sm_go/SM_GO_annot.csv")

#remove excess columns
sm_go <- sm_go %>% dplyr::select(gene, GO)

#Convert this dataframe into a txt file so that it can be used for topGO
write.table(sm_go, file = "df_GO.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
#load as geneID2GO format for topGO
geneID2GO <- readMappings(file = "df_GO.txt") 
#remove the file, we no longer need it
file.remove(file = "df_GO.txt")

#get list of genes differentially expressed in either host specific direction
myInterestingGenes <- deg %>% filter(significance != "Not Significant") %>%
	pull(gene)

# Create new geneList 
geneNames <- names(geneID2GO)

gene_list <- factor(as.integer(geneNames %in% myInterestingGenes))
names(gene_list) <- geneNames

# topGO data object
GOdata <- new("topGOdata",
							ontology = "MF",
							allGenes = gene_list,
							annot = annFUN.gene2GO,
							gene2GO = geneID2GO)

# Enrichment Test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Analysis of Results
hostspec_results <- GenTable(GOdata, classicFisher = resultFisher)
hostspec_results$geneset <- "host_specific"

#now run for non- host specific
#get list of genes differentially expressed in either host specific direction
myInterestingGenes <- deg %>% filter(significance == "Not Significant") %>%
	pull(gene)

# Create new geneList 
geneNames <- names(geneID2GO)

gene_list <- factor(as.integer(geneNames %in% myInterestingGenes))
names(gene_list) <- geneNames

# topGO data object
GOdata <- new("topGOdata",
							ontology = "MF",
							allGenes = gene_list,
							annot = annFUN.gene2GO,
							gene2GO = geneID2GO)

# Enrichment Test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Analysis of Results
nonhostspec_results <- GenTable(GOdata, classicFisher = resultFisher)
nonhostspec_results$geneset <- "non_host_specific"

result <- rbind(hostspec_results, nonhostspec_results)

#save
result %>% write.csv("data/bcin_RNAseq/bcin_DEG/SM_GO/bcin_hostspec_SMGO_enrich.csv", row.names = F)

#### ------------- plotting SM gene classes for DEG groups -------------------

sm <- read.csv("../FabBc_Coex/output_preMay2025/select_networks/bcin_SM/bcinSMgenes_Suarezetal.csv")

unique(sm$class)

pks <- sm %>% filter(class == "PKS") %>% pull(gene)
pks_genes <- deg %>% filter(gene %in% pks)

nrps <- sm %>% filter(class == "NRPS") %>% pull(gene)
nrps_genes <- deg %>% filter(gene %in% nrps)

stc <- sm %>% filter(class == "STC") %>% pull(gene)
stc_genes <- deg %>% filter(gene %in% stc)

dtc <- sm %>% filter(class == "DTC") %>% pull(gene)
dtc_genes <- deg %>% filter(gene %in% dtc)

dmat <- sm %>% filter(class == "DMAT") %>% pull(gene)
dmat_genes <- deg %>% filter (gene %in% dmat)



#### ============= Network analysis ===================================

#Want to see if any genes in either host-specific group co-express in networks

#Define function to get list of clusters where ANY within a defined list of genes appear
extract_clusters_anygene <- function(net, genelist) {
	# Filter the dataframe to rows where 'gene' is in corr_neg_genes
	# Select only the 'gene' and 'clustfull_host' columns
	extracted_rows <- net %>%
		filter(gene %in% genelist) %>%
		select(gene, cluster_full)
	return(extracted_rows)
}

#Define function to get list of clusters where at least 3 genes in the list appear
extract_clusters_w3genes <- function(net, genelist) {
	# First, extract genes and clustfull_host
	initial_result <- net %>%
		filter(gene %in% genelist) %>%
		select(gene, cluster_full)
	# Count occurrences of each clustfull_host
	# Then filter to keep only those with count > 3
	filtered_result <- initial_result %>%
		group_by(cluster_full) %>%
		filter(n() > 3) %>%
		ungroup()
	return(filtered_result)
}

#Define function to get list of clusters where at least 5 genes in the list appear
extract_clusters_w5genes <- function(net, genelist) {
	# First, extract genes and clustfull_host
	initial_result <- net %>%
		filter(gene %in% genelist) %>%
		select(gene, cluster_full)
	# Count occurrences of each clustfull_host
	# Then filter to keep only those with count > 5
	filtered_result <- initial_result %>%
		group_by(cluster_full) %>%
		filter(n() > 5) %>%
		ungroup()
	return(filtered_result)
}

### Pv
#get list of genes upregulated
pv_upreg <- deg %>% filter(significance == "Common Bean Upregulated")
pv_upreg_list <- pv_upreg$gene

#load significant networks
pv_net <- read.csv("~/UCDavis/Klieb_Lab/Projects/RNASeq/full_pipeline/remove_7isos/reformatted_networks/Pv/all_clusters.csv")
pv_net <- pv_net %>% filter(P.value < 0.05)

#extract any clusters where these genes appear
result <- extract_clusters_anygene(pv_net, pv_upreg_list)
unique(result$cluster_full) %>% length()
#there are 254

result <- extract_clusters_w3genes(pv_net, pv_upreg_list)
unique(result$cluster_full) %>% length()
#there are 88

result <- extract_clusters_w5genes(pv_net, pv_upreg_list)
unique(result$cluster_full) %>% length()
#there are 61

#saving d5 networks with 5+ host specific genes
netlist <- unique(result$cluster_full)
filtered_nets <- pv_net %>% filter(cluster_full %in% netlist)
filtered_nets <- filtered_nets %>% filter(decay_rate == 5)
filtered_nets %>% write.csv("data/bcin_RNAseq/bcin_DEG/Pvnetworks_wHostSpecGenes.csv",
														row.names = F)

### Vu
#get list of genes upregulated
vu_upreg <- deg %>% filter(significance == "Cowpea Upregulated")
vu_upreg_list <- vu_upreg$gene

#load significant networks
vu_net <- read.csv("~/UCDavis/Klieb_Lab/Projects/RNASeq/full_pipeline/remove_7isos/reformatted_networks/Vu/all_clusters.csv")
vu_net <- vu_net %>% filter(P.value < 0.05)

#extract any clusters where these genes appear
result <- extract_clusters_anygene(vu_net, vu_upreg_list)
unique(result$cluster_full) %>% length()
#there are 240

result <- extract_clusters_w3genes(vu_net, vu_upreg_list)
unique(result$cluster_full) %>% length()
#there are 73

result <- extract_clusters_w5genes(vu_net, vu_upreg_list)
unique(result$cluster_full) %>% length()
#there are 35

#saving d5 networks with 5+ host specific genes
netlist <- unique(result$cluster_full)
filtered_nets <- vu_net %>% filter(cluster_full %in% netlist)
filtered_nets <- filtered_nets %>% filter(decay_rate == 5)
filtered_nets %>% write.csv("data/bcin_RNAseq/bcin_DEG/Vunetworks_wHostSpecGenes.csv",
														row.names = F)


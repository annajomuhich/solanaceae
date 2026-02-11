### Bcin PCA - comparison of transcriptomes on pepper, tomato, and lettuce 
### February 2026 AJM


### ------------- B05.10 only --------------------------------------------
#prepare pepper data
ca_bcin <- read.csv("data/norm_counts/Pepper_B.cinerea_expression_long.csv") %>%
	filter(iso_name == "B05_10") %>%
	mutate(host = "Pepper",
				 timepoint = 48,
				 count_log2 = log2(CPM + 1)) %>%
	dplyr::select(host, timepoint, gene, count_log2)

#prepare tomato data
sl_bcin <- read.csv("data/norm_counts/Tomato_B.cinerea_expression_long.csv") %>%
	filter(iso_name == "B05_10") %>%
	mutate(host = "Tomato",
				 timepoint = 48,
				 count_log2 = log2(CPM + 1)) %>%
	dplyr::select(host, timepoint, gene, count_log2)

#prepare lettuce data
df <- read.csv("../Lettuce_Denby/bot_gene_count_sizefactored_wideformat.csv") %>%
	dplyr::rename(gene = X) %>%
	pivot_longer(cols = !gene,
							 names_to = "sample_ID",
							 values_to = "count") %>%
	mutate(treatment = str_replace(sample_ID,
																 "(.*)_(.*)_(.*)",
																 "\\1"),
				 timepoint = str_replace(sample_ID,
				 												"(.*)_(.*)_(.*)",
				 												"\\2")) %>%
	filter(treatment == "Infected") %>%
	mutate(host = "Lettuce",
				 timepoint = str_sub(timepoint, 3, 4),
				 count_log2 = log2(count + 1)) %>%
	dplyr::select(host, timepoint, gene, count_log2)
	
df <- rbind(ca_bcin, sl_bcin, df)

#calculate means
means_pca <- df %>%
	group_by(host, timepoint, gene) %>%
	summarise(mean_count_log2 =  mean(count_log2)) %>%
	ungroup() %>%
	mutate(sample_ID = paste0(host, "_", timepoint))

#reformat means for PCA

#pivot wider
means_pca <- means_pca %>%
	pivot_wider(values_from = mean_count_log2,
							names_from = gene)
head(means_pca)

any(is.na(means_pca))

#remove NA genes
na_cols <- colSums(is.na(means_pca)) > 0
removed_cols <- names(means_pca)[na_cols]
message("removing ", length(removed_cols), " genes with missing data")
means_pca <- means_pca[, !na_cols]
any(is.na(means_pca))

#remove categorical variables
means_pca_matrix <- means_pca %>%
	dplyr::select(!c(timepoint, host, sample_ID))
head(means_pca_matrix)
#put in gene rownames
row.names(means_pca_matrix) <- means_pca$sample_ID

# Create a DGEList object
#dge <- DGEList(counts = means_pca_matrix)

pca_result <- prcomp(means_pca_matrix, scale. = FALSE)

#take a look
pca_summary <- summary(pca_result)
pca_summary
#extract PC1 and PC2 proportions of variance for axis labels
pca_summary <- as.data.frame(pca_summary$importance)
pc1 <- pca_summary["Proportion of Variance", "PC1"]
pc2 <- pca_summary["Proportion of Variance", "PC2"]
pc1 <- pc1 * 100
pc2 <- pc2 * 100
pc1 <- signif(pc1, 3) #reduce to 3 significant figures
pc2 <- signif(pc2, 3)

#convert to dataframe
pca_df <- as.data.frame(pca_result$x)
#initial plot
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

#append host geno and iso information
pca_df <- cbind(means_pca[,c(1,2)], pca_df)

#PC1 x PC2
pca_df %>% ggplot(aes(PC1, PC2)) +
	geom_point(aes(color = host), size = 3) +
	geom_text_repel(data=pca_df,
									aes(label = timepoint),
									position = position_jitter(width = 5, height = 5)) +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	theme_minimal() +
	scale_color_brewer(palette = "Dark2")
ggsave("figures/bcin_transcriptome/PCA/wLettuce/B0510_wLettuce.png", width = 6, height = 5)

### ------------- all isolates --------------------------------------------
#prepare pepper data
ca_bcin <- read.csv("data/norm_counts/Pepper_B.cinerea_expression_long.csv") %>%
	#filter(iso_name == "B05_10") %>%
	mutate(host = "Pepper",
				 timepoint = 48,
				 count_log2 = log2(CPM + 1)) %>%
	dplyr::select(host, iso_name, timepoint, gene, count_log2)

#prepare tomato data
sl_bcin <- read.csv("data/norm_counts/Tomato_B.cinerea_expression_long.csv") %>%
	#filter(iso_name == "B05_10") %>%
	mutate(host = "Tomato",
				 timepoint = 48,
				 count_log2 = log2(CPM + 1)) %>%
	dplyr::select(host, iso_name, timepoint, gene, count_log2)

#prepare lettuce data
df <- read.csv("../Lettuce_Denby/bot_gene_count_sizefactored_wideformat.csv") %>%
	dplyr::rename(gene = X) %>%
	pivot_longer(cols = !gene,
							 names_to = "sample_ID",
							 values_to = "count") %>%
	mutate(treatment = str_replace(sample_ID,
																 "(.*)_(.*)_(.*)",
																 "\\1"),
				 timepoint = str_replace(sample_ID,
				 												"(.*)_(.*)_(.*)",
				 												"\\2")) %>%
	filter(treatment == "Infected") %>%
	mutate(host = "Lettuce",
				 iso_name = "B05_10",
				 timepoint = str_sub(timepoint, 3, 4),
				 count_log2 = log2(count + 1)) %>%
	dplyr::select(host, iso_name, timepoint, gene, count_log2)

df <- rbind(ca_bcin, sl_bcin, df)

#calculate means
means_pca <- df %>%
	group_by(host, iso_name, timepoint, gene) %>%
	summarise(mean_count_log2 =  mean(count_log2)) %>%
	ungroup() %>%
	mutate(sample_ID = paste0(host, "_", iso_name, "_", timepoint))

#reformat means for PCA

#pivot wider
means_pca <- means_pca %>%
	pivot_wider(values_from = mean_count_log2,
							names_from = gene)
head(means_pca)

any(is.na(means_pca))

#remove NA genes
na_cols <- colSums(is.na(means_pca)) > 0
removed_cols <- names(means_pca)[na_cols]
message("removing ", length(removed_cols), " genes with missing data")
means_pca <- means_pca[, !na_cols]
any(is.na(means_pca))

#remove categorical variables
means_pca_matrix <- means_pca %>%
	dplyr::select(!c(timepoint, host, iso_name, sample_ID))
head(means_pca_matrix)
#put in gene rownames
row.names(means_pca_matrix) <- means_pca$sample_ID

# Create a DGEList object
#dge <- DGEList(counts = means_pca_matrix)

pca_result <- prcomp(means_pca_matrix, scale. = FALSE)

#take a look
pca_summary <- summary(pca_result)
pca_summary
#extract PC1 and PC2 proportions of variance for axis labels
pca_summary <- as.data.frame(pca_summary$importance)
pc1 <- pca_summary["Proportion of Variance", "PC1"]
pc2 <- pca_summary["Proportion of Variance", "PC2"]
pc1 <- pc1 * 100
pc2 <- pc2 * 100
pc1 <- signif(pc1, 3) #reduce to 3 significant figures
pc2 <- signif(pc2, 3)

#convert to dataframe
pca_df <- as.data.frame(pca_result$x)
#initial plot
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

#append host geno and iso information
pca_df <- cbind(means_pca[,c(1,2,3)], pca_df)

#PC1 x PC2
pca_df %>%
	ggplot(aes(PC1, PC2)) +
	geom_point(aes(color = host), size = 3) +
	geom_text_repel(
		data = dplyr::filter(pca_df, host == "Lettuce"),
		aes(label = timepoint),
		position = position_jitter(width = 5, height = 5)
	) +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	theme_minimal() +
	scale_color_brewer(palette = "Dark2")
ggsave("figures/bcin_transcriptome/PCA/wLettuce/allisos_wLettuce.png", width = 6, height = 5)


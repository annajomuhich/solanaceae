### Solanaceae
### Gene filtering of normalized counts - Botrytis
### February 2026 - AJM

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

#look at gene expression across all samples
df_long %>%
	group_by(gene) %>%
	summarise(prop_expressed = mean(CPM >= 1), .groups = "drop") %>%
	ggplot(aes(x = prop_expressed)) +
	geom_histogram(binwidth = 0.01) +
	xlab("proportion of samples with expression") +
	ylab("number of genes")

#keep only genes that have >=1 CPM in at least 20% of the samples.
removed_genes <- df_long %>%
	group_by(gene) %>%
	summarise(prop_expressed = mean(CPM >= 1), .groups = "drop") %>%
	filter(prop_expressed < 0.20) %>%
	pull(gene)

df_long$gene %>% unique() %>% length()
#filter
df_long_filtered <- df_long %>%
		filter(!gene %in% removed_genes)
df_long_filtered$gene %>% unique() %>% length()

df_long_filtered %>%
	group_by(gene) %>%
	summarise(prop_expressed = mean(CPM >= 1), .groups = "drop") %>%
	ggplot(aes(x = prop_expressed)) +
	geom_histogram(binwidth = 0.01) +
	xlab("proportion of samples with expression") +
	ylab("number of genes")

#split back into host files
ca <- df_long_filtered %>%
	filter(genotype == "Ca")
sl <- df_long_filtered %>%
	filter(genotype == "Sl")

#write them out
ca %>%
	write.csv("data/norm_counts/Pepper_B.cinerea_expression_long_filtered.csv", row.names = F)
sl %>%
	write.csv("data/norm_counts/Tomato_B.cinerea_expression_long_filtered.csv", row.names = F)

# ### Determine which genes to remove ======================================
# df <- read.csv("data/bcin_expr/bcin_genemodelNB_20260205/bcin_adjusted_emmeans.csv")
# canorm <- read.csv("data/norm_counts/Pepper_B.cinerea_expression_long.csv") %>%
# 	mutate(host = "Pepper") %>%
# 	dplyr::select(host, iso_name, rep, gene, CPM)
# slnorm <- read.csv("data/norm_counts/Tomato_B.cinerea_expression_long.csv") %>%
# 	mutate(host = "Tomato") %>%
# 	dplyr::select(host, iso_name, rep, gene, CPM)
# norm <- rbind(canorm, slnorm)
# 
# #Quick look at overall gene count
# ggplot(df, aes(x = emmean)) +
# 	geom_histogram(bins = 500) +
# 	ggtitle("gene count histogram - no transformation")
# #There's kind of a big tail.
# 
# #for genes in tail cluster, check mean count per group and percent zeros per group
# #first get list of genes in tail cluster
# df_sums <- df %>%
# 	group_by(gene) %>%
# 	summarise(sum_expr = sum(emmean))
# #histogram
# ggplot(df_sums, aes(x = sum_expr)) +
# 	geom_histogram(bins = 500) +
# 	ggtitle("gene sum histogram - no transformation")
# 
# #removing all genes with sums <-1500
# genes_to_remove <- df_sums %>%
# 	filter(sum_expr < -1000) %>%
# 	pull(gene)
# 

# norm %>%
# norm_wide <- norm %>% pivot_wider(names_from = host, values_from = CPM)
# norm_filt <- norm %>%
# 	filter(gene %in% genes_to_remove )
# 
# norm_filt_summary <- norm_filt %>%
# 	group_by(host, gene) %>%
# 	summarise(
# 		n_samples = n(),
# 		n_zero = sum(CPM == 0, na.rm = TRUE),
# 		pct_zero = 100 * n_zero / n_samples,
# 		mean_count = mean(CPM, na.rm = TRUE)
# 	) %>%
# 	ungroup()
# 
# #there are 36 of them
# df <- df %>%
# 	filter(!gene %in% genes_to_remove)
# 
# #looking at each genes totals again
# df_sums <- df %>%
# 	group_by(gene) %>%
# 	summarise(sum_expr = sum(emmean))
# #histogram
# ggplot(df_sums, aes(x = sum_expr)) +
# 	geom_histogram(bins = 500) +
# 	ggtitle("gene sum histogram - 36 gene tail removed")
# #this looks much better.
# ggplot(df, aes(x = emmean)) +
# 	geom_histogram(bins = 500) +
# 	ggtitle("gene count histogram - 36 gene tail removed")
# 
# ### Replace unfiltered datasets ===============================================
# 
# #write the filtered dataset out
# df %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_adjusted_emmeans.csv", row.names = F)
# 
# #filter anova
# anova <- read.csv("data/bcin_expr/bcin_genemodel_20260115/unfiltered/bcin_anova.csv")
# anova <- anova %>%
# 	filter(!gene %in% genes_to_remove)
# anova %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv", row.names = F)
# 
# #filter degs
# deg <- read.csv("data/bcin_expr/bcin_genemodel_20260115/unfiltered/bcin_DEGs.csv")
# deg <- deg %>%
# 	filter(!gene %in% genes_to_remove)
# deg %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_DEGs.csv", row.names = F)

### Solanaceae
### Gene filtering of normalized counts - Host orthologs
### February 2026 - AJM

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

#make infected column
df_long <- df_long %>%
	mutate(infected = if_else(iso_name == "Mock_48h", "no", "yes")) %>%
	dplyr::select(genotype, iso_name, infected, everything())

### Remove genes not expressed in both hosts ======================
df <- df_long %>%
	pivot_wider(names_from = ortholog,
							values_from = cpm)
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
df_long$ortholog %>% unique() %>% length()
df_long <- df_long %>%
	filter(!ortholog %in% removed_cols)
df_long$ortholog %>% unique() %>% length()

#look at gene expression across all samples
df_long %>%
	group_by(ortholog) %>%
	summarise(prop_expressed = mean(cpm >= 1), .groups = "drop") %>%
	ggplot(aes(x = prop_expressed)) +
	geom_histogram(binwidth = 0.01) +
	xlab("proportion of samples with expression") +
	ylab("number of genes")

#keep only orthologs that have >=1 CPM in at least 20% of the samples.
removed_orthologs <- df_long %>%
	group_by(ortholog) %>%
	summarise(prop_expressed = mean(cpm >= 1), .groups = "drop") %>%
	filter(prop_expressed < 0.20) %>%
	pull(ortholog)

df_long$ortholog %>% unique() %>% length()
#filter
df_long_filtered <- df_long %>%
		filter(!ortholog %in% removed_orthologs)
df_long_filtered$ortholog %>% unique() %>% length()

df_long_filtered %>%
	group_by(ortholog) %>%
	summarise(prop_expressed = mean(cpm >= 1), .groups = "drop") %>%
	ggplot(aes(x = prop_expressed)) +
	geom_histogram(binwidth = 0.01) +
	xlab("proportion of samples with expression") +
	ylab("number of orthologs")

# write out filtered ortholog expression
df_long_filtered %>%
	write.csv("data/norm_counts/PepperTomato_ortho_expression_long_filtered.csv", row.names = F)


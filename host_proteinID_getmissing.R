### Get LOC IDs for protein IDs missing from the existing keys
### January 2026 AJM

### --------------- Get list of XP IDs with missing LOC IDs for both hosts ------------------------------
ca <- read.csv("data/gene_descriptions/fromRitu_notusing/Pepper_Gene_Transcript_Protein_Mapping.csv")
sl <- read.csv("data/gene_descriptions/fromRitu_notusing/Tomato_Gene_Transcript_Protein_Mapping.csv")

ca_missing <- ca %>%
		filter(!startsWith(gene_id, "LOC"))
sl_missing <- sl %>%
	filter(!startsWith(gene_id, "LOC"))

xp_missing <- rbind(ca_missing, sl_missing)

xp_missing %>% write.csv("data/gene_descriptions/host_proteinID_getmissing/xp_missing.csv", row.names = F)

### -------------- Join missing data with existing key --------------------------------------------------
missing <- read.csv("data/gene_descriptions/host_proteinID_getmissing/CaSl_loc_missing.csv") %>%
	rename(new_gene_id = gene_id)

cadf <- left_join(ca, missing, by = "protein_id")
cadf <- cadf %>%
	dplyr::mutate(
		gene_id = if_else(startsWith(gene_id, "LOC"), gene_id, new_gene_id)
	)

sldf <- left_join(sl, missing, by = "protein_id")
sldf <- sldf %>%
	dplyr::mutate(
		gene_id = if_else(startsWith(gene_id, "LOC"), gene_id, new_gene_id)
	)

cadf <- cadf %>%
	select(gene_id, transcript_id, protein_id)
sldf <- sldf %>%
	select(gene_id, transcript_id, protein_id)

cadf %>% write.csv("data/gene_descriptions/pepper_proteinID_key.csv", row.names = F)
sldf %>% write.csv("data/gene_descriptions/tomato_proteinID_key.csv", row.names = F)

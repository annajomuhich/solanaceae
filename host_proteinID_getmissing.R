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

#xp_missing %>% write.csv("data/gene_descriptions/host_proteinID_getmissing/xp_missing.csv", row.names = F)

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

#cadf %>% write.csv("data/gene_descriptions/pepper_proteinID_key.csv", row.names = F)
#sldf %>% write.csv("data/gene_descriptions/tomato_proteinID_key.csv", row.names = F)

### --------- Get list of XP IDs missing from the key entirely -------------------------
#There are a bunch of these missing so when I try to join my ortholog data to LOC, I get NAs.
#Of all the unique values in the raw orthofinder output, how many of those are missing from my keys?
#Probably should have just started with this but whatever

#read in orthologs
ortho <- read_delim("data/ortho/orthofinder/Orthologues/Orthologues_GCF_002878395.1_UCD10Xv1.1_protein/GCF_002878395.1_UCD10Xv1.1_protein__v__GCF_036512215.1_SLM_r2.1_protein.tsv") %>%
	dplyr::rename(Ca = GCF_002878395.1_UCD10Xv1.1_protein,
								Sl = GCF_036512215.1_SLM_r2.1_protein)

#get list of unique protein IDs
unique_xp_ids <- unique(
	unlist(
		strsplit(
			na.omit(c(ortho$Ca, ortho$Sl)),
			",\\s*"
		)
	)
)

xp_ids_available <- rbind(cadf, sldf) %>%
	pull(protein_id) %>%
	unique()

xp_ids_missing <- setdiff(unique_xp_ids, xp_ids_available)

#xp_ids_missing %>% write.csv("data/gene_descriptions/host_proteinID_getmissing/xp_missing2.csv", row.names = F)

### ----------------- Join missing data with existing key ------------------------------------------------

missing <- read.csv("data/gene_descriptions/host_proteinID_getmissing/CaSl_loc_missing2.csv") %>%
	mutate(transcript_id = NA)

#Just going to bind them altogether. This means a few of the other host will be in here but isn't a big deal bc there are no replicates
cadf <- rbind(cadf, missing)
unique(cadf$protein_id) %>% length()

sldf <- rbind(sldf, missing)
unique(sldf$protein_id) %>% length()

cadf %>% write.csv("data/gene_descriptions/pepper_proteinID_key.csv", row.names = F)
sldf %>% write.csv("data/gene_descriptions/tomato_proteinID_key.csv", row.names = F)
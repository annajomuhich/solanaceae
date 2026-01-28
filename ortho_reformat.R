### Reformatting Orthofinder output
### January 2026 AJM

library(readr)
library(tidyr)
library(dplyr)

### ----- Reformat orthofinder --------------------

#read in orthologs
ortho <- read_delim("data/ortho/orthofinder/Orthologues/Orthologues_GCF_002878395.1_UCD10Xv1.1_protein/GCF_002878395.1_UCD10Xv1.1_protein__v__GCF_036512215.1_SLM_r2.1_protein.tsv")

#rename columns for simplicity
ortho <- ortho %>%
	dplyr::rename(Ca = GCF_002878395.1_UCD10Xv1.1_protein,
								Sl = GCF_036512215.1_SLM_r2.1_protein)

#remove NA values
ortho$Ca <- vapply(
	strsplit(ortho$Ca, ",\\s*"),
	\(x) paste(x[x != "NA"], collapse = ", "),
	character(1)
)
ortho$Sl <- vapply(
	strsplit(ortho$Sl, ",\\s*"),
	\(x) paste(x[x != "NA"], collapse = ", "),
	character(1)
)

ortho %>% filter(!startsWith(Sl, "LOC")) %>% pull(Sl) %>% unique()

#get protein ID to gene ID
ca_ids <- read.csv("data/gene_descriptions/Pepper_Gene_Transcript_Protein_Mapping.csv") %>%
	dplyr::select(!transcript_id)
sl_ids <- read.csv("data/gene_descriptions/Tomato_Gene_Transcript_Protein_Mapping.csv") %>%
	dplyr::select(!transcript_id)

ca <- ortho %>%
	mutate(row_id = row_number()) %>%
	separate_rows(Ca, sep = ",\\s*") %>%
	left_join(ca_ids, by = c("Ca" = "protein_id")) %>%
	group_by(row_id) %>%
	summarise(
		gene_id = paste(unique(gene_id), collapse = ", "),
		.groups = "drop") %>%
	rename(Ca = gene_id) %>%
	dplyr::select(!row_id)

sl <- ortho %>%
	mutate(row_id = row_number()) %>%
	separate_rows(Sl, sep = ",\\s*") %>%
	left_join(sl_ids, by = c("Sl" = "protein_id")) %>%
	group_by(row_id) %>%
	summarise(
		gene_id = paste(unique(gene_id), collapse = ", "),
		.groups = "drop") %>%
	rename(Sl = gene_id) %>%
	dplyr::select(!row_id)

ortho <- ortho %>% select(Orthogroup)
ortho <- cbind(ortho, ca, sl)

#there are many groups that have multiple orthologs in either species
#let's split these into diff dataframes to keep
#Retain just the 1 to 1 orthologs
ortho_1to1 <- ortho[!grepl(",", ortho$Ca) & !grepl(",", ortho$Sl), ]
ortho_1tomany <- ortho[!grepl(",", ortho$Ca) & grepl(",", ortho$Sl), ]
ortho_manyto1 <- ortho[grepl(",", ortho$Ca) & !grepl(",", ortho$Sl), ]
ortho_manytomany <- ortho[grepl(",", ortho$Ca) & grepl(",", ortho$Sl), ]

#There are some NAs in the 1to1 orthologs. Let's remove these



#write these out for further analyses
ortho %>% write.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_all.csv", row.names = F)
ortho_1to1 %>% write.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_1to1.csv", row.names=F)
ortho_1tomany %>% write.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_1tomany.csv", row.names=F)
ortho_manyto1 %>% write.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_manyto1.csv", row.names=F)
ortho_manytomany %>% write.csv("data/ortho/orthofinder/reformatted_orthologs/CaSl_ortho_manytomany.csv", row.names=F)

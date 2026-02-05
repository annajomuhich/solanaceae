### Update host expression files
### January 2026 AJM

#host files do not have entirely numerical gene IDs
#I'd like to change this before further analysis
#also generally reformat them to match

library(tidyverse)

### ------------------ Pepper -------------------------------
ca <- read.csv('data/norm_counts/fromRitu_notusing/Pepper_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv')

#pepper is entirely LOC IDs
ca %>%
	filter(!startsWith(gene, "LOC")) %>%
	pull(gene) %>%
	unique()

colnames(ca)

ca <- ca %>% 	select(sample_ID,
										 iso_number,
										 iso_name,
										 rep,
										 tray,
										 seq_batch,
										 gene,
										 cpm)

ca %>% write.csv("data/norm_counts/Pepper_Host_expression.csv", row.names = F)

### ------------------ Tomato -------------------------------
sl <- read.csv('data/norm_counts/fromRitu_notusing/Tomato_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv')

#tomato has a bunch of other stuff
sl %>%
	filter(!startsWith(gene, "LOC")) %>%
	pull(gene) %>%
	unique()

annot_sl <- read.csv("data/gene_descriptions/tomato_annotations.csv")

#need to manually change a bunch of mismatching names so they match what's in the annotation:
sl <- sl %>%
	mutate(gene = sub("LEAS1", "LeAS1", gene)) %>%
	mutate(gene = sub("AMT1;1", "AMT1-1", gene)) %>%
	mutate(gene = sub("LePK4", "PK4,PKv", gene)) %>%
	mutate(gene = sub("LeZIP2", "LeZIP2,SlbZIP1", gene)) %>%
	mutate(gene = sub("cRK_1", "CRK1", gene)) %>%
	mutate(gene = sub("Ls", "LS", gene)) %>%
	filter(gene != "LecRK_2") %>% #couldn't find a reliable match for this one
	filter(gene != "LecRK") #theres a ton of these and it could be any so just removing

#matches are across multiple columns in the annotation.
#Joining and coalescing
sl_joined <- sl %>%
	left_join(
		annot_sl %>% select(gene_ID, Name),
		by = c("gene" = "Name")
	) %>%
	left_join(
		annot_sl %>% select(gene_ID_syn = gene_ID, gene_synonym),
		by = c("gene" = "gene_synonym")
	) %>%
	left_join(
		annot_sl %>% select(gene_ID_loctag = gene_ID, locus_tag),
		by = c("gene" = "locus_tag")
	) %>%
	mutate(
		gene_ID = coalesce(gene_ID, gene_ID_syn, gene_ID_loctag),
		gene_ID = if_else(
			is.na(gene_ID) & grepl("^LOC", gene),
			gene,
			gene_ID
		)
	) %>%
	select(-c(gene_ID_syn, gene_ID_loctag)) %>%
	mutate(gene = sub("^Le", "", gene)) %>%
	left_join(
		annot_sl %>% select(gene_ID_name = gene_ID, Name),
		by = c("gene" = "Name")) %>%
	mutate(
		gene_ID = coalesce(gene_ID, gene_ID_name)
	) %>%
	select(-gene_ID_name)

#Check no gene_IDs that are NA
sl_joined %>%
	filter(is.na(gene_ID)) %>%
	pull(gene) %>% unique()

#remove "sum_" string from sample ID and reorder columns
sl_joined <- sl_joined %>%
	mutate(sample_ID = gsub("sum_", "", sample_ID)) %>%
	select(-gene) %>%
	rename(gene = gene_ID) %>%
	select(sample_ID,
				 iso_number,
				 iso_name,
				 rep,
				 tray,
				 seq_batch,
				 gene,
				 cpm)

sl_joined %>% write.csv("data/norm_counts/Tomato_Host_expression.csv", row.names = F)

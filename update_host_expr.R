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

ca <- ca %>% select(sample_ID,
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
annot_names <- annot_sl %>%	select(Name, gene_ID)

sl_join <- left_join(sl, annot_names, by = join_by("gene" == "Name"))
sl_join <- sl_join %>%
	select(!gene) %>%
	rename(gene = gene_ID)

colnames(sl_join)

sl_join <- sl_join %>% select(sample_ID,
															iso_number,
															iso_name,
															rep,
															tray,
															seq_batch,
															gene,
															cpm)

sl_join %>% write.csv("data/norm_counts/Tomato_Host_expression.csv", row.names = F)

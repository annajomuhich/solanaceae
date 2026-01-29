###### Dataset setup for mr2mods from emmeans
###### January 2026 AJM

library(tidyverse)
library(stringr)

###### Pepper - Make normalized.matrix ==================================
#load data
host <- read.csv("data/norm_counts/Pepper_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv") %>%
	select(gene, sample_ID, cpm)
host$sample_ID <- gsub("_Host", "", host$sample_ID)
bcin <- read.csv("data/norm_counts/Pepper_B.cinerea_expression_long.csv") %>%
	select(gene, expr_sample_ID, CPM) %>%
	rename(cpm = CPM,
				 sample_ID = expr_sample_ID)

#bind dataframes together if colnames match
if (all(colnames(host) == colnames(bcin))) {
	df <- rbind(host, bcin)
}

#remove tailing value in sample ID
#I am not sure what this signifies, would have to ask Ritu, but I'm guessing it denotes a backup sample or something
df$sample_ID <- gsub("^((?:[^_]*_){2}[^_]*).*", "\\1", df$sample_ID)

#pivot wide
#need sample ID as column names and genes as rownames for mr2mods!
df <- df %>%
	pivot_wider(names_from = sample_ID,
							values_from = cpm)

#look for NAs (needs to be FALSE)
any(is.na(df))

#convert gene column to rownames
df <- column_to_rownames(df, var = "gene")

#write matrix csv
df %>% write.table("data/gcn/input/CaBc_normalized.matrix", sep="\t", quote = FALSE)

###### Tomato - Make normalized.matrix ==================================

### start here - ran into a bunch of funkiness with the sample IDs/genes here

#load data
host <- read.csv("data/norm_counts/Tomato_Host_TMM_CPM_Merged_SampleID_48hpiMock.csv") %>%
	select(gene, sample_ID, cpm)
host$sample_ID <- gsub("_Host", "", host$sample_ID)
bcin <- read.csv("data/norm_counts/Tomato_B.cinerea_expression_long.csv") %>%
	select(gene, expr_sample_ID, CPM) %>%
	rename(cpm = CPM,
				 sample_ID = expr_sample_ID)

#bind dataframes together if colnames match
if (all(colnames(host) == colnames(bcin))) {
	df <- rbind(host, bcin)
}

#remove tailing value in sample ID
#I am not sure what this signifies, would have to ask Ritu, but I'm guessing it denotes a backup sample or something
df$sample_ID <- gsub("^((?:[^_]*_){2}[^_]*).*", "\\1", df$sample_ID)

#pivot wide
#need sample ID as column names and genes as rownames for mr2mods!
df <- df %>%
	pivot_wider(names_from = sample_ID,
							values_from = cpm)

#look for NAs (needs to be FALSE)
any(is.na(df))

#convert gene column to rownames
df <- column_to_rownames(df, var = "gene")

#write matrix csv
df %>% write.table("data/gcn/input/CaBc_normalized.matrix", sep="\t", quote = FALSE)
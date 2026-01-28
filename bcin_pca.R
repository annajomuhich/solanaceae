### Solanaceae
### Comparison of Botrytis reads with PCA
### January 2026 - AJM

library(tidyverse)
library(gridExtra)
library(edgeR)
library(readxl)
library(RColorBrewer)
library(ggrepel)

# ------------- PCA - 2026-01-21 --------------------
df <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_adjusted_emmeans.csv") %>%
	dplyr::select(!SE)

#Check for any NAs (needs to be FALSE)
any(is.na(df))

#pivot wide
df <- df %>%
	pivot_wider(names_from = gene,
							values_from = emmean)

#Check for any NAs (needs to be FALSE)
any(is.na(df))
colSums(is.na(df) > 0) #see by column


######### PCA

#get a numerical dataframe
dfpca <- df[,-c(1,2)]

#run PCA
pca_result <- prcomp(dfpca, scale. = FALSE)

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
pca_df <- cbind(df[,c(1,2)], pca_df)

pca_df <- pca_df %>%
	mutate(species = if_else(genotype == "Ca", "Pepper", "Tomato")) %>% 
	select(iso_name, species, genotype, everything())

#PC1 x PC2
pca_plot <- pca_df %>% ggplot(aes(PC1, PC2)) +
	geom_point(aes(color = species), size = 3) +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	theme_minimal() +
	scale_color_manual(values = c("navy","red"))

pca_plot

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_geno.png", width = 6, height = 5)

# ----------------- PCA colored by lesion size ------------------

#read in pheno data
pheno <- read.csv("data/Pheno/Alleudi_long_format_72iso.csv") %>%
	filter(Genotype == "LA3475" | Genotype == "Ca3") %>%
	select(Isolate_name, Taxon, LsMeans_mm2) %>%
	rename(iso_name = Isolate_name,
				 species = Taxon,
				 Lesion_size = LsMeans_mm2)

#join pheno data to PCA table
pca_pheno_df <- left_join(pca_df, pheno, by = c("iso_name", "species")) %>%
	select(iso_name, species, genotype, Lesion_size, everything())

pca_pheno_df$Lesion_size <- as.integer(pca_pheno_df$Lesion_size)

pca_pheno_plot <- pca_pheno_df %>% ggplot(aes(PC1, PC2, fill = Lesion_size)) +
	geom_point(size = 3, shape = 21, stroke = 0) +
	scale_fill_viridis_c(option = "C", direction = -1) +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	labs(fill = "Lesion Size") +
	theme_minimal()

pca_pheno_plot

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_Pheno.png", width = 6, height = 5)

####### PCA with both host species colored by lesion size ##########

pca_pheno_df %>%
	ggplot(aes(PC1, PC2, fill = Lesion_size, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	labs(fill = expression("Lesion Size (mm"^2*")")) +
	theme_minimal()

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_genopheno.png", width = 5, height = 4)

# #Colored by log10 lesion size
# pca_pheno_df %>%
# 	ggplot(aes(PC1, PC2, fill = Lesion_size, shape = species)) +
# 	geom_point(size = 3, stroke = 0.2, color = "black") +
# 	scale_shape_manual(values = c(21, 24)) +
# 	scale_fill_viridis_c(option = "D", direction = -1, trans = "log10") +
# #	ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
# 	xlab(paste0("PC1 - ", pc1, "%")) +
# 	ylab(paste0("PC2 - ", pc2, "%")) +
# 	labs(fill = expression("Lesion Size log10 (mm"^2*")")) +
# 	theme_minimal()

#ggsave("figures/bcin_transcriptome/bcin_PCA_log10.png", width = 5, height = 4)

##### 9/23/25 Statistical test - PCs correlation with lesion size --------------------------------------

model <- lm(Lesion_size ~ PC1 + PC2, data = pca_pheno_df)

anova(model)

#10/28/25 checking species too
model <- lm(Lesion_size ~ PC1 * species, data = pca_pheno_df)
anova(model)
model <- lm(Lesion_size ~ PC2 * species, data = pca_pheno_df)
anova(model)

####### Plotting Other PCs ############

#setup axis labels
pc3 <- pca_summary["Proportion of Variance", "PC3"]
pc4 <- pca_summary["Proportion of Variance", "PC4"]
pc5 <- pca_summary["Proportion of Variance", "PC5"]
pc3 <- pc3 * 100
pc4 <- pc4 * 100
pc5 <- pc5 * 100
pc3 <- signif(pc3, 3) #reduce to 3 significant figures
pc4 <- signif(pc4, 3)
pc5 <- signif(pc5, 3)

pca_pheno_df %>%
	ggplot(aes(PC1, PC3, fill = Lesion_size, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC3 - ", pc3, "%")) +
	labs(fill = expression("Lesion Size (mm"^2*")")) +
	theme_minimal()

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_PC3.png", width = 5, height = 4)

pca_pheno_df %>%
	ggplot(aes(PC1, PC4, fill = Lesion_size, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC4 - ", pc4, "%")) +
	labs(fill = expression("Lesion Size (mm"^2*")")) +
	theme_minimal()

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_PC4.png", width = 5, height = 4)

pca_pheno_df %>%
	ggplot(aes(PC1, PC5, fill = Lesion_size, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC5 - ", pc4, "%")) +
	labs(fill = expression("Lesion Size (mm"^2*")")) +
	theme_minimal()

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_PC5.png", width = 5, height = 4)



###### Colored by Botrytis read count ##########

#Read in readcount data
peppercount <- read.csv("data/bcin_expr/transcript_abundance/Pepper_Bcin_Percent_Across_Isolates.csv")
tomatocount <- read.csv("data/bcin_expr/transcript_abundance/Tomato_Bcin_Percent_Across_Isolates.csv")
peppercount$species <- "Pepper"
tomatocount$species <- "Tomato"
peppercount <- peppercount %>%
	rename(iso_name = Species_Isolate)
peppercount <- peppercount %>% dplyr::select(species, iso_name, Total_Reads_Bcin, Percent_Bcin)
tomatocount <- tomatocount %>% dplyr::select(species, iso_name, Total_Reads_Bcin, Percent_Bcin)
count <- rbind(peppercount, tomatocount)
#rename columns to match downstream code
count <- count %>% dplyr::rename(bcin_count = "Total_Reads_Bcin",
																 bcin_pct = "Percent_Bcin")


#join pheno data to PCA table
pca_count_df <- left_join(pca_df, count, by = c("iso_name", "species")) %>%
	select(iso_name, species, genotype, bcin_count, bcin_pct, everything())

pca_count_df %>%
	ggplot(aes(PC1, PC2, fill = bcin_count, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	labs(fill = "Botrytis \ntranscript \nabundance") +
	theme_minimal()

ggsave("figures/bcin_transcriptome/PCA/bcin_PCA_rdcount.png", width = 5, height = 4)

#check botrytis read percentage
pca_count_df %>%
	ggplot(aes(PC1, PC2, fill = bcin_pct, shape = species)) +
	geom_point(size = 2.5, stroke = 0.2, color = "black") +
	scale_shape_manual(values = c(21, 24)) +
	scale_fill_viridis_c(option = "D", direction = -1) +
	#ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	labs(fill = "Botrytis read percentage") +
	theme_minimal()

# #label points with iso name
# pca_count_df %>%
# 	ggplot(aes(PC1, PC2, fill = bcin_count, shape = species)) +
# 	geom_point(size = 2.5, stroke = 0.2, color = "black") +
# 	geom_text_repel(aes(label = ifelse(PC1 < -50, iso_name, NA_character_)),
# 									size = 3, max.overlaps = 100, na.rm = TRUE) +
# 	scale_shape_manual(values = c(21, 24)) +
# 	scale_fill_viridis_c(option = "D", direction = -1) +
# 	# ggtitle("PCA of Botrytis cinerea reads on\nFabales species by lesion size") +
# 	xlab(paste0("PC1 - ", pc1, "%")) +
# 	ylab(paste0("PC2 - ", pc2, "%")) +
# 	labs(fill = "Botrytis read count") +
# 	theme_minimal()

# #ggsave("figures/bcin_transcriptome/bcin_PCA_rdcount_lowlabeled.png", width = 5, height = 4)

### 9/23/25 statistical test - PCs correlation with total read count -------------------------------------3

model <- lm(bcin_count ~ PC1 + PC2, data = pca_count_df)

anova(model)

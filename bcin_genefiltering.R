### Solanaceae
### Post-model gene filtering
### January 2026 - AJM


### Determine which genes to remove ======================================
df <- read.csv("data/bcin_expr/bcin_genemodel_20260115/unfiltered/bcin_adjusted_emmeans.csv")

#Quick look at overall gene count
ggplot(df, aes(x = emmean)) +
	geom_histogram(bins = 500) +
	ggtitle("gene count histogram - no transformation")
#There's kind of a big tail.

#looking at each genes totals.
df_sums <- df %>%
	group_by(gene) %>%
	summarise(sum_expr = sum(emmean))
#histogram
ggplot(df_sums, aes(x = sum_expr)) +
	geom_histogram(bins = 500) +
	ggtitle("gene sum histogram - no transformation")
#let's remove the big tail from this.

#removing all genes with sums <-1500
genes_to_remove <- df_sums %>%
	filter(sum_expr < -1500) %>%
	pull(gene)
#there are 36 of them
df <- df %>%
	filter(!gene %in% genes_to_remove)

#looking at each genes totals again
df_sums <- df %>%
	group_by(gene) %>%
	summarise(sum_expr = sum(emmean))
#histogram
ggplot(df_sums, aes(x = sum_expr)) +
	geom_histogram(bins = 500) +
	ggtitle("gene sum histogram - 36 gene tail removed")
#this looks much better.
ggplot(df, aes(x = emmean)) +
	geom_histogram(bins = 500) +
	ggtitle("gene count histogram - 36 gene tail removed")

### Replace unfiltered datasets ===============================================

#write the filtered dataset out
df %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_adjusted_emmeans.csv", row.names = F)

#filter anova
anova <- read.csv("data/bcin_expr/bcin_genemodel_20260115/unfiltered/bcin_anova.csv")
anova <- anova %>%
	filter(!gene %in% genes_to_remove)
anova %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv", row.names = F)

#filter degs
deg <- read.csv("data/bcin_expr/bcin_genemodel_20260115/unfiltered/bcin_DEGs.csv")
deg <- deg %>%
	filter(!gene %in% genes_to_remove)
deg %>% write.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_DEGs.csv", row.names = F)

###### Bcin gene p values - Upset plot
###### AJM January 2026

library(tidyverse)
library(ComplexUpset)

##### ------------- Set up data -------------------------
df <- read.csv("data/bcin_expr/bcin_genemodel_20260115/bcin_anova.csv") %>%
	dplyr::select(gene, variable, p_adj) %>%
	filter(variable != "intercept")

any(is.na(df))
colSums(is.na(df) > 0) #see by column


df_wide <- pivot_wider(df, names_from = "variable", values_from = "p_adj")

# convert p-values to significance (TRUE/FALSE)
sig_df <- df_wide %>%
	mutate(
		`Host Species` = genotype < 0.05,
		`Isolate` = iso_name < 0.05,
		`Host Species X Isolate` = `genotype:iso_name` < 0.05
	) %>%
	select(gene, `Host Species`, `Isolate`, `Host Species X Isolate`)

sig_df <- sig_df %>% rename("Host Species\nX Isolate" = `Host Species X Isolate`)

### Plot ============================================================
upset(
	sig_df,
	intersect = c("Host Species", "Isolate", "Host Species\nX Isolate"),
	annotations = list(),
	base_annotations = list(
		'Significant Botrytis Genes' = intersection_size(text = list(size = 3))
	),
	width_ratio = c(0, 1),        # removes left panel entirely
	set_sizes = FALSE             # <- this removes the "Set Size" title
) +
	theme_minimal() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x  = element_blank(),
		axis.ticks.x = element_blank()
	)
ggsave("figures/bcin_transcriptome/bcin_anova_upset.png", width = 4, height = 4)

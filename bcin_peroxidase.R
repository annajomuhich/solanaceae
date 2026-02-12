### Investigating host-specific peroxidase expression in Botrytis
### February 2026 AJM

GO <- read.csv("data/bcin_expr/DEG/NB/CaUp_GO.csv")
deg <- read.csv("data/bcin_expr/DEG/NB/bcin_DEG.csv")
annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv")
emmeans <- read.csv("data/bcin_expr/bcin_genemodelNB_20260205/bcin_adjusted_emmeans.csv")
anova <- read.csv("data/bcin_expr/bcin_genemodelNB_20260205/bcin_anova.csv")

emmeans %>%
	filter(genotype == "Ca") %>%
	ggplot(aes(x = emmean)) +
	geom_histogram(binwidth = 0.1)

emmeans %>%
	filter(genotype == "Sl") %>%
	ggplot(aes(x = emmean)) +
	geom_histogram(binwidth = 0.1)

GO$Term
GO[,1:2]

#strong ROS association GO terms
strong <- c("GO:0016491","GO:0016705", "GO:0016684", "GO:0016209", "GO:0004601")

peroxidase <- "GO:0004601"

strong_genes <- annot %>%
	filter(
		if_any(
			c(
				X.Computed.GO.Function.IDs.,
				X.Computed.GO.Process.IDs.),
			~str_detect(., str_c(peroxidase, collapse = "|"))
		)
	)

strong_genes <- strong_genes %>%
	dplyr::select(X.Gene.ID.,
								X.Computed.GO.Function.IDs.,
								X.Computed.GO.Process.IDs.,
								X.Gene.Name.or.Symbol., X.PFam.Description.
								)

strong_geneIDs <- strong_genes %>% pull(X.Gene.ID.)

converged <- anova %>%
	filter(convergence_note == "relative convergence (4)") %>%
	pull(gene) %>%
	unique()

emmeans %>%
	filter(gene %in% converged) %>%
	filter(emmean > -10) %>%
	dplyr::select(!SE) %>%
	filter(gene %in% strong_geneIDs) %>%
	pivot_wider(names_from = genotype,
							values_from = emmean) %>%
	ggplot(aes(x = Ca, y = Sl)) +
	geom_point() +
	ylim(-2, 9) +
	xlim(-2, 9)

# strong_genes <- annot %>%
# 	filter(
# 		if_any(
# 			c(X.Computed.GO.Component.IDs.,
# 				X.Computed.GO.Function.IDs.,
# 				X.Computed.GO.Process.IDs.,
# 				X.Curated.GO.Component.IDs.,
# 				X.Curated.GO.Function.IDs.,
# 				X.Curated.GO.Process.IDs.),
# 			~str_detect(., str_c(strong, collapse = "|"))
# 		)
# 	)
# 
# strong_genes <- strong_genes %>%
# 	dplyr::select(X.Gene.ID.,
# 								X.Computed.GO.Function.IDs.,
# 								X.Computed.GO.Process.IDs.,
# 								X.Curated.GO.Component.IDs.,
# 								X.Curated.GO.Function.IDs.,
# 								X.Curated.GO.Process.IDs.)


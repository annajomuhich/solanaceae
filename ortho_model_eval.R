### Evaluating ortholog model outputs
### February 2026 AJM

nb <- read.csv("data/ortho/ortho_model_20260203/ortho_anova.csv")
gaus <- read.csv("data/ortho/ortho_modelGAUSSIAN_20260203/ortho_anova.csv")

table(nb$convergence_note)
table(gaus$convergence_note)

gaus %>%
	filter(variable != "intercept") %>%
	ggplot(aes(x = p_adj)) +
	geom_histogram(binwidth = .01) +
	facet_wrap("variable")
ggsave("figures/ortho_model_eval/gaussian_pvals.png", height = 3, width = 7.5)

nb %>%
	filter(variable != "intercept") %>%
	ggplot(aes(x = p_adj)) +
	geom_histogram(binwidth = .01) +
	facet_wrap("variable")
ggsave("figures/ortho_model_eval/nb_pvals.png", height = 3, width = 7.5)

read.csv("data/ortho/ortho_model_20260203/failed_genes.csv") %>% nrow()
read.csv("data/ortho/ortho_modelGAUSSIAN_20260203/failed_genes.csv") %>% nrow()

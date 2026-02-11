### Evaluating ortholog model outputs
### February 2026 AJM

nb <- read.csv("data/bcin_expr/bcin_genemodelNB_20260205/bcin_anova.csv")
gaus <- read.csv("data/bcin_expr/bcin_genemodelGAUSSIAN_20260205/bcin_anova.csv")

table(nb$convergence_note)
table(gaus$convergence_note)

gaus %>%
	filter(variable != "intercept") %>%
	ggplot(aes(x = p_adj)) +
	geom_histogram(binwidth = .01) +
	facet_wrap("variable")
ggsave("figures/bcin_model_eval/gaussian_pvals.png", height = 3, width = 7.5)

nb %>%
	filter(variable != "intercept") %>%
	ggplot(aes(x = p_adj)) +
	geom_histogram(binwidth = .01) +
	facet_wrap("variable")
ggsave("figures/bcin_model_eval/nb_pvals.png", height = 3, width = 7.5)

read.csv("data/bcin_expr/bcin_genemodelNB_20260205/failed_genes.csv") %>% nrow()
read.csv("data/bcin_expr/bcin_genemodelGAUSSIAN_20260205/failed_genes.csv") %>% nrow()

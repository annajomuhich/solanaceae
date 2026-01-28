### Presence absence matrix of Botrytis GO terms host specific and iso specific
### January 2026 AJM

go_host <- read.csv("data/bcin_expr/host_iso_specific/hostallUp_GO.csv") %>%
	pull(Term)
go_iso <- read.csv("data/bcin_expr/host_iso_specific/isoallUp_GO.csv") %>%
	pull(Term)
go_interaction <- read.csv("data/bcin_expr/host_iso_specific/interactionallUp_GO.csv") %>%
	pull(Term)
all_go <- sort(unique(c(go_host, go_iso, go_interaction)))

go_pa <- data.frame(GO = all_go,
										Host.specific = all_go %in% go_host,
										Isolate.specific = all_go %in% go_iso,
										Interaction.specific = all_go %in% go_interaction)

go_long <- pivot_longer(
	go_pa,
	cols = -GO,
	names_to = "Variable",
	values_to = "Present"
)

go_long <- go_long %>%
	mutate(Present = if_else(Present == TRUE, "yes", "no")) %>%
	mutate(Variable = if_else(Variable == "Host.specific", "Host specific", 
														if_else(Variable == "Isolate.specific", "Isolate specific", "Interaction specific")))

ggplot(go_long, aes(x = Variable, y = GO, fill = Present)) +
	geom_tile(color = "grey80") +
	scale_fill_manual(
		values = c("no" = "white", "yes" = "dark orange"),
		name = "Upregulated"
	) +
	theme_minimal() +
	theme(
		axis.text.y = element_text(size = 7),
		panel.grid = element_blank(),
		legend.position = "none",
		axis.text.x = element_text(angle = 45, hjust = 1)
	) +
	scale_y_discrete(limits = rev) +
	ylab("") +
	xlab("")

ggsave("figures/bcin_transcriptome/bcin_allvar_GO.png", height = 5, width = 6)

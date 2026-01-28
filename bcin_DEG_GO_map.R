### Presence absence matrix of Botrytis GO terms upregulated on either host
### January 2026 AJM

df <- read.csv("data/bcin_expr/DEG/GO_integrated.csv")

go_pepper <- df %>%
	filter(host_upregulated == "Pepper") %>%
	pull(Term)
go_tomato <- df %>%
	filter(host_upregulated == "Tomato") %>%
	pull(Term)
all_go <- sort(unique(df$Term))

go_pa <- data.frame(GO = all_go,
										Pepper = all_go %in% go_pepper,
										Tomato = all_go %in% go_tomato)

go_long <- pivot_longer(
	go_pa,
	cols = -GO,
	names_to = "Host",
	values_to = "Present"
)

go_long <- go_long %>%
	mutate(Present = if_else(Present == TRUE, "yes", "no"))

ggplot(go_long, aes(x = Host, y = GO, fill = Present)) +
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

ggsave("figures/bcin_transcriptome/bcin_DEG_GO.png", height = 5, width = 6.75)

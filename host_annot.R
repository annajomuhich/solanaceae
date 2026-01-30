### Improving host annotation files and fixing related count inputs
### January 2026 AJM

### ------------------ Pepper -------------------------------------------------
ca <- read.csv("data/gene_descriptions/fromRitu_notusing/pepper_cleaned_gene_annotations.csv")

#issues for pepper:
#Only have written descriptions for ~300 genes
#No GO terms
#no LOC ID for a few genes

#Make LOC ID column
ca$gene_ID <- ca$Dbxref
ca$gene_ID <- gsub("GeneID:", "", ca$gene_ID)

#Check these are all numerical values as expected (TRUE)
all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", ca$gene_ID))

ca <- ca %>%
	mutate(gene_ID = paste0("LOC", gene_ID))

#These columns are exactly the same so we can get rid of one
all(ca$Name == ca$gene)
ca <- ca %>% select(!gene)

### ------------------ Tomato -------------------------------------------------
sl <- read.csv("data/gene_descriptions/fromRitu_notusing/tomato_cleaned_gene_annotations.csv")


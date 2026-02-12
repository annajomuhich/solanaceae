### Improving host annotation files and fixing related count inputs
### January 2026 AJM

library(dplyr)

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
ca <- ca %>% dplyr::select(!gene)

#There is no variation in these columns so can remove
ca <- ca %>% dplyr::select(!c("score", "phase", "gbkey"))

#I am thinking these are not very informative so cleaning these out too
ca <- ca %>% dplyr::select(!c("partial", "start_range", "end_range", "exception"))

#Changing column to numerical ID in case I need it
ca$ID_number <- ca$Dbxref
ca$ID_number <- gsub("GeneID:", "", ca$ID_number)

#read in GO annot
go <- read.delim(
		"data/gene_descriptions/GO/pepper/GCF_002878395.1_UCD10Xv1.1_gene_ontology.gaf",
		header = F,
		sep = "\t",
		comment.char = "!",
		stringsAsFactors = FALSE
	)
colnames(go) <- c("DB",
									"GeneID",
									"Symbol",
									"Qualifier",
									"GO_ID",
									"Reference",
									"Evidence_Code",
									"With_From",
									"Aspect",
									"Gene_Name",
									"Gene_Synonym",
									"Type",
									"Taxon",
									"Date",
									"Assigned_By",
									"Annot_Ext",
									"Gene_Product_Form_ID")
#get a column for all GO terms
go_by_gene <- go %>%
	dplyr::group_by(GeneID) %>%
	dplyr::summarise(
		GO_all = paste(unique(GO_ID), collapse = ", "),
		.groups = "drop"
	)
#get separate columns for each GO aspect too
go_split <- go %>%
	dplyr::group_by(GeneID, Aspect) %>%
	dplyr::summarise(
		GO_terms = paste(unique(GO_ID), collapse = ", "),
		.groups = "drop"
	) %>%
	tidyr::pivot_wider(
		names_from = Aspect,
		values_from = GO_terms,
		names_prefix = "GO_"
	)
go_by_gene <- left_join(go_by_gene, go_split, by = "GeneID")
#join GO info with annotation
go_by_gene$GeneID <- as.character(go_by_gene$GeneID)
ca <- left_join(ca, go_by_gene, by = join_by("ID_number" == "GeneID"))

#adding protein IDs. Removing additional protein isoforms because I don't think they have a lot of different annotation info
prot <- read.csv("data/gene_descriptions/pepper_proteinID_key.csv") %>%
	dplyr::select(!transcript_id) %>%
	dplyr::rename(protein_ID = protein_id) %>%
	distinct(gene_id, .keep_all = T) #just keeping the first occurrence of each gene ID/protein isoform
ca <- left_join(ca, prot, by = join_by("gene_ID" == "gene_id"))

#next: add blastp descriptions
#some blastp results haev multiple matches. Keep only the one with the highest percent identity
cavsl <- read.csv("data/gene_descriptions/blastp/Pepper_v_Tomato_blastp.csv") %>%
		group_by(qseqid) %>%                 # group by qseqid
		slice_max(tomato_pident, n = 1) %>%  # keep row(s) with max
		ungroup() %>%
		distinct(qseqid, .keep_all = T) #this is to get rid of any remaining ones that might have the same % identity
cavat <- read.csv("data/gene_descriptions/blastp/Pepper_v_Arabidopsis_blastp.csv") %>%
	group_by(qseqid) %>%                 # group by qseqid
	slice_max(arabidopsis_pident, n = 1) %>%  # keep row(s) with max
	ungroup() %>%
	distinct(qseqid, .keep_all = T) %>% #this is to get rid of any remaining ones that might have the same % identity
	mutate(
		arabidopsis_symbol = arabidopsis_description %>% #separate the symbol and description string
			str_extract("\\|[^|]*\\|") %>%
			str_remove_all("\\|") %>%
			str_remove("^\\s*Symbols:\\s*") %>%
			str_trim() %>%
			na_if(""),
		
		arabidopsis_description = arabidopsis_description %>%
			str_remove("\\|[^|]*\\|\\s*") %>%
			str_trim()
	)
blastp <- full_join(cavsl, cavat, by = "qseqid")

ca <- left_join(ca, blastp, by = join_by(protein_ID == qseqid))

#writing for now
ca %>% write.csv("data/gene_descriptions/pepper_annotations.csv", row.names = F)

### ------------------ Tomato -------------------------------------------------
rm(list=ls())

sl <- read.csv("data/gene_descriptions/fromRitu_notusing/tomato_cleaned_gene_annotations.csv") %>%
	distinct(Dbxref, .keep_all=TRUE) #get rid of duplicate

#Make LOC ID column
sl$gene_ID <- sl$Dbxref
sl$gene_ID <- gsub("GeneID:", "", sl$gene_ID)

#there's a bunch that end in a weird string so removing these
sl <- sl %>%
	mutate(gene_ID = sub(",.*", "", gene_ID))

#Check these are all numerical values as expected (TRUE)
all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", sl$gene_ID))

sl <- sl %>%
	mutate(gene_ID = paste0("LOC", gene_ID))

#These columns are exactly the same so we can get rid of one
all(sl$Name == sl$gene)
sl <- sl %>% dplyr::select(!gene)

#There is no variation in these columns so can remove
sl <- sl %>% dplyr::select(!c("score", "phase", "gbkey"))

#I am thinking these are not very informative so cleaning these out too
sl <- sl %>% dplyr::select(!c("partial", "start_range", "end_range", "exception"))

#Changing column to numerical ID in case I need it
sl$ID_number <- sl$Dbxref
sl$ID_number <- gsub("GeneID:", "", sl$ID_number)
sl <- sl %>%
	mutate(ID_number = sub(",.*", "", ID_number))
all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", sl$ID_number))

#read in GO annot
go <- read.delim(
	"data/gene_descriptions/GO/tomato/GCF_036512215.1-RS_2024_09_gene_ontology.gaf.gz",
	header = F,
	sep = "\t",
	comment.char = "!",
	stringsAsFactors = FALSE
)
colnames(go) <- c("DB",
									"GeneID",
									"Symbol",
									"Qualifier",
									"GO_ID",
									"Reference",
									"Evidence_Code",
									"With_From",
									"Aspect",
									"Gene_Name",
									"Gene_Synonym",
									"Type",
									"Taxon",
									"Date",
									"Assigned_By",
									"Annot_Ext",
									"Gene_Product_Form_ID")
#get a column for all GO terms
go_by_gene <- go %>%
	dplyr::group_by(GeneID) %>%
	dplyr::summarise(
		GO_all = paste(unique(GO_ID), collapse = ", "),
		.groups = "drop"
	)
#get separate columns for each GO aspect too
go_split <- go %>%
	dplyr::group_by(GeneID, Aspect) %>%
	dplyr::summarise(
		GO_terms = paste(unique(GO_ID), collapse = ", "),
		.groups = "drop"
	) %>%
	tidyr::pivot_wider(
		names_from = Aspect,
		values_from = GO_terms,
		names_prefix = "GO_"
	)
go_by_gene <- left_join(go_by_gene, go_split, by = "GeneID")
#join GO info with annotation
go_by_gene$GeneID <- as.character(go_by_gene$GeneID)
sl <- left_join(sl, go_by_gene, by = join_by("ID_number" == "GeneID"))

#adding protein IDs. Removing additional protein isoforms because I don't think they have a lot of different annotation info
prot <- read.csv("data/gene_descriptions/tomato_proteinID_key.csv") %>%
	dplyr::select(!transcript_id) %>%
	dplyr::rename(protein_ID = protein_id) %>%
	distinct(gene_id, .keep_all = T) #just keeping the first occurrence of each gene ID/protein isoform

sl <- left_join(sl, prot, by = join_by("gene_ID" == "gene_id"))

#Add arabidopsis descriptions
#Some arabidopsis blastp results have multiple matches. Keep only the one with the highest percent identity
slvat <- read.csv("data/gene_descriptions/blastp/Tomato_v_Arabidopsis_blastp.csv") %>%
	group_by(qseqid) %>%                 # group by qseqid
	slice_max(arabidopsis_pident, n = 1) %>%  # keep row(s) with max arabidopsis_pident
	ungroup() %>%
	distinct(qseqid, .keep_all = T) %>% #this is to get rid of any remaining ones that might have the same % identity
	mutate(
		arabidopsis_symbol = arabidopsis_description %>% #separate the symbol and description string
			str_extract("\\|[^|]*\\|") %>%
			str_remove_all("\\|") %>%
			str_remove("^\\s*Symbols:\\s*") %>%
			str_trim() %>%
			na_if(""),
		
		arabidopsis_description = arabidopsis_description %>%
			str_remove("\\|[^|]*\\|\\s*") %>%
			str_trim()
	)

sl <- left_join(sl, slvat, by = join_by(protein_ID == qseqid))

#writing for now
sl %>% write.csv("data/gene_descriptions/tomato_annotations.csv", row.names = F)

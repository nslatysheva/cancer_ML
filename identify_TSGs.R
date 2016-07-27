# Identify TF genes and TF families
require(dplyr)
require(stringr)
require(tidyr)

# Get data
tsg <- read.table("~/Projects/fusion_ML/data/gene_identity_data/TSGs/tsgs_from_tsgene.txt", header=TRUE, sep='\t')
downreg_tsg <- read.table("~/Projects/fusion_ML/data/gene_identity_data/TSGs/downregulated_tsgs.txt", header=TRUE, sep='\t')
nrow(tsg); nrow(downreg_tsg); head(downreg_tsg)

# combine tsg tables
tsg <- merge(tsg, select(downreg_tsg, -gene_symbol), by="entrez", all.x=TRUE); head(tsg); nrow(tsg)

# Convert to ensg
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)
tsg_genes <- merge(tsg, entrez_ensg, on="entrez", all.x = TRUE); head(tsg_genes, 500); nrow(tsg_genes)

# Anything we can do to assign ENSGs to the remaining genes? Down from ~210 to 72 NAs, not much else to do
filter(tsg_genes, is.na(ensg)) %>% select(entrez)

# Any original entries now missing? All ok.
n_distinct(tsg$gene_symbol); n_distinct(tsg_genes$gene_symbol)
tsg$gene_symbol[!tsg$gene_symbol %in% unique(tsg_genes$gene_symbol)]

# Any original entries duplicated? This is fine.
tsg_genes %>% group_by(gene_symbol) %>%
  summarise(
    n_ensg = n_distinct(ensg)
  ) %>% arrange(-n_ensg)

filter(tsg_genes, gene_symbol == "MDC1")

head(tsg_genes)

# reformat and write out gene list
tsg_list_n_cancer_count <- 
  tsg_genes %>%
  select(ensg, gene_type, n_cancers) %>%
  replace_na(list(n_cancers=0)) %>%
  na.omit() 

# dummy variables for type of gene
tsg_type_dummy <- cbind(tsg_summary, model.matrix(~ gene_type - 1, tsg_summary))
head(tsg_type_dummy)

tsg_summary <- 
  tsg_type_dummy %>% 
  select(ensg, n_cancers, `gene_typeprotein-coding`, gene_typencRNA) %>%
  rename(num_cancers_downreg_TSG = n_cancers) %>%
  mutate(is_TSG = 1) %>% select(ensg, is_TSG, n_cancers_downreg_TSG)

# # check counts of dummy variables
# tsg_summary %>% 
#   select(`gene_typeprotein-coding`, gene_typencRNA) %>%
#   mutate_each(funs(as.factor)) %>% summary

# write out tsg table
nrow(tsg_summary); nrow(distinct(tsg_summary))
head(tsg_summary)

write.table(select(tsg_summary, ensg, is_TSG, n_cancers_downreg_TSG), 
            "~/Projects/fusion_ML/features/tsg_genes.csv", sep=",", quote = FALSE, row.names = FALSE)

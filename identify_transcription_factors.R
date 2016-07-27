# Identify TF genes and TF families
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)

# Get tf data
tf <- read.table("~/Projects/fusion_ML/data/gene_identity_data/transcriptional_control/trrust_rawdata.txt", header=TRUE, sep='\t'); head(tf)

# Get tg gene symbol to entrez labels
tf_entrez <- read.table("~/Projects/fusion_ML/data/gene_identity_data/transcriptional_control/trrust_gene_symbol_entrez.csv", header=TRUE, sep=','); head(tf_entrez)

# Convert gene symbol to entrez
tf <- merge(tf, tf_entrez, on="gene_symbol"); head(tf)
tf <- select(tf, gene_symbol, entrez, target, interaction_type) 
head(tf,20); nrow(tf)

# Convert to ensg
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)
tf_genes <- merge(tf, entrez_ensg, on="entrez", all.x = TRUE); head(tf_genes, 500); nrow(tf_genes); head(tf_genes)

# Any original entries now missing? All ok.
tf$gene_symbol[!tf$gene_symbol %in% unique(tf_genes$gene_symbol)]

# Any original entries duplicated? This is fine.
n_distinct(tf$gene_symbol); n_distinct(tf_genes$ensg)

head(tf_genes)

tf_genes %>% group_by(gene_symbol) %>%
  summarise(
    n_ensg = n_distinct(ensg)
  ) %>% arrange(-n_ensg)

filter(tf_genes, gene_symbol == "ZNRD1")

# Anything we can do to assign ENSGs to the remaining genes? All sorted now.
table(filter(tf_genes, is.na(ensg)) %>% select(entrez))

# write out tf gene list
tf_genes$is_TF <- 1; head(tf_genes)
write.table(na.omit(distinct(select(tf_genes, ensg, is_TF))), 
            "~/Projects/fusion_ML/features/tf_genes.csv", sep=",", quote = FALSE, row.names = FALSE)

# get total number of interactions per gene
n_interactions <- tf_genes %>%
  mutate(
    interaction_id = paste(ensg, interaction_type, target, sep="_")
  ) %>% 
  group_by(ensg) %>%
  summarise(
    num_TF_interactions= n_distinct(interaction_id)
  )

### type of interactions per gene
head(tf_genes,50); nrow(tf_genes)

n_tf_interaction_type <-
  tf_genes %>%
  mutate(
    interaction_id = paste(ensg, interaction_type, target, sep="_")
  ) %>% 
  group_by(ensg, interaction_type) %>%
  summarise(
    interaction_count = n_distinct(interaction_id)
  )

head(n_tf_interaction_type, 30); nrow(n_tf_interaction_type)

# Get counts of activation and repression interactions by gene
act_rep_counts <-
  n_tf_interaction_type %>%
  filter(interaction_type != "Unknown") %>%
  spread(key=interaction_type, value=interaction_count) %>%
  replace_na(list(Activation=0, Repression=0))

head(act_rep_counts)

# write out tf table
head(n_interactions); head(n_tf_interaction_type); head(act_rep_counts)
nrow(n_interactions); nrow(n_tf_interaction_type); nrow(act_rep_counts)

tf_counts <-
  merge(n_interactions, act_rep_counts, by="ensg", all.x = TRUE) %>%
  replace_na(list(Activation=0, Repression=0)) %>%
  rename(num_TF_activation = Activation,
         num_TF_repression = Repression)

head(tf_counts)

write.table(tf_counts, "~/Projects/fusion_ML/features/tf_classes.csv", sep=",", quote = FALSE, row.names = FALSE)

# Are two TF tables needed?
head(na.omit(distinct(select(tf_genes, ensg, is_tf)))); nrow(na.omit(distinct(select(tf_genes, ensg, is_tf))))
head(tf_counts); nrow(tf_counts)








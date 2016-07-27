# GO multifunctionality

require(dplyr)
require(tidyr)

# Get data
go <- read.table("~/Projects/fusion_ML/data/gene_identity_data/GO/GO_terms_by_ensg.txt", header=TRUE, sep='\t')
go_slim <- read.table("~/Projects/fusion_ML/data/gene_identity_data/GO/GOSlim_terms_by_ensg.txt", header=TRUE, sep='\t')

head(go_slim)

# GO term count
go_clean <- 
  go %>% group_by(Ensembl.Gene.ID) %>%
  summarise(
    num_GO_terms = n_distinct(GO.Term.Accession))

# GO slim term count 
go_slim_clean <- 
  go_slim %>% rename(GOSlim_accession = GOSlim.GOA.Accession.s.) %>%
  group_by(Ensembl.Gene.ID) %>%
  summarise(
    num_GOSlim_terms = n_distinct(GOSlim_accession))

# merge the two datasets
head(go_clean); head(go_slim_clean); nrow(go_clean); nrow(go_slim_clean)
go <- full_join(go_clean, go_slim_clean, by = "Ensembl.Gene.ID") %>% rename(ensg = Ensembl.Gene.ID)
go[is.na(go)] <- 0; head(go)

# There are actually more GoSlim terms :P
sum(go$num_GO_terms)
sum(go$num_GOSlim_terms)

# write out features
write.table(go, "~/Projects/fusion_ML/features/GO_counts.csv", sep=",", quote = FALSE, row.names = FALSE)

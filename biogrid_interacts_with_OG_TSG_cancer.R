# biogrid interactions with cancer genes, OGs, TSGs

require(dplyr)
require(igraph)

# load gene sets
oncogene <- read.table("~/Projects/fusion_ML/features/oncogene_list.csv", header=TRUE, sep=','); head(oncogene)
tsg <- read.table("~/Projects/fusion_ML/features/tsg_genes.csv", header=TRUE, sep=','); head(tsg)
cancer_genes <- read.table("~/Projects/fusion_ML/features/cancer_genes.csv", header=TRUE, sep=','); head(cancer_genes)

# load biogrid
d <- read.csv("~/Projects/fusion_ML/data/network_data/biogrid_human_phys.csv", sep=","); head(d)
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)

d <- merge(merge(d, entrez_ensg, by.x="Entrez.Gene.Interactor.A", by.y="entrez", all.x = TRUE),
           entrez_ensg, by.x="Entrez.Gene.Interactor.B", by.y="entrez", all.x = TRUE) %>%
  rename(entrez_a = Entrez.Gene.Interactor.A,
         entrez_b = Entrez.Gene.Interactor.B,
         ensg_a = ensg.x,
         ensg_b = ensg.y) %>% 
  select(entrez_a, entrez_b, ensg_a, ensg_b); head(d)

# omit any remaining rows with NA values
d <- na.omit(d); d <- d %>% select(ensg_a, ensg_b) %>% distinct(); head(d)

# create graph
g <- graph.data.frame(d, directed = FALSE)

# get genes interacting with oncogene
interacts_with_OG_1 <- d %>% filter(d$ensg_b %in% oncogene$ensg) %>% rename(ensg=ensg_a, OG=ensg_b); head(interacts_with_OG_1)
interacts_with_OG_2 <- d %>% filter(d$ensg_a %in% oncogene$ensg) %>% rename(ensg=ensg_b, OG=ensg_a) %>% select(ensg, OG); head(interacts_with_OG_2)
interacts_with_OG <- arrange(rbind(interacts_with_OG_1, interacts_with_OG_2), ensg, OG); head(interacts_with_OG,100)

interacts_with_OG_clean <- 
  interacts_with_OG %>%
  group_by(ensg) %>% 
  summarise(num_OG_interactions = n_distinct(OG)) %>% 
  mutate(interacts_with_OG = 1)

write.table(interacts_with_OG_clean, 
            "~/Projects/fusion_ML/features/interacts_with_OG.csv", sep=",", quote = FALSE, row.names = FALSE)

# get genes interacting with TSGs
interacts_with_TSG_1 <- d %>% filter(d$ensg_b %in% tsg$ensg) %>% rename(ensg=ensg_a, TSG=ensg_b); head(interacts_with_TSG_1)
interacts_with_TSG_2 <- d %>% filter(d$ensg_a %in% tsg$ensg) %>% rename(ensg=ensg_b, TSG=ensg_a) %>% select(ensg, TSG); head(interacts_with_TSG_2)
interacts_with_TSG <- arrange(rbind(interacts_with_TSG_1, interacts_with_TSG_2), ensg, TSG); head(interacts_with_TSG,100)

interacts_with_TSG_clean <- 
  interacts_with_TSG %>%
  group_by(ensg) %>% 
  summarise(num_TSG_interactions = n_distinct(TSG)) %>% 
  mutate(interacts_with_TSG = 1)

write.table(interacts_with_TSG_clean, 
            "~/Projects/fusion_ML/features/interacts_with_TSG.csv", sep=",", quote = FALSE, row.names = FALSE)

# get genes interacting with cancer genes
interacts_with_cancer_1 <- d %>% filter(d$ensg_b %in% cancer_genes$ensg) %>% rename(ensg=ensg_a, cancer=ensg_b); head(interacts_with_cancer_1)
interacts_with_cancer_2 <- d %>% filter(d$ensg_a %in% cancer_genes$ensg) %>% rename(ensg=ensg_b, cancer=ensg_a) %>% select(ensg, cancer); head(interacts_with_cancer_2)
interacts_with_cancer <- arrange(rbind(interacts_with_cancer_1, interacts_with_cancer_2), ensg, cancer); head(interacts_with_cancer,100)

interacts_with_cancer_clean <- 
  interacts_with_cancer %>%
  group_by(ensg) %>% 
  summarise(num_cancer_interactions = n_distinct(cancer)) %>% 
  mutate(interacts_with_cancer = 1)

write.table(interacts_with_cancer_clean, 
            "~/Projects/fusion_ML/features/interacts_with_cancer.csv", sep=",", quote = FALSE, row.names = FALSE)


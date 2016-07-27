# Processing biogrid data

require(dplyr)

biogrid <- read.csv("~/Projects/fusion_ML/data/network_data/BIOGRID-ALL-3.4.135.tab2.txt", sep="\t")
head(biogrid); nrow(biogrid)

# get human interactions
human <- filter(biogrid, Organism.Interactor.A == 9606 & Organism.Interactor.B == 9606)
head(human); nrow(human)

# How many genes are covered now?
n_distinct(c(human$Entrez.Gene.Interactor.A, human$Entrez.Gene.Interactor.B))

# get only interactions identified by physical experimental systems
human_phys <- filter(human, Experimental.System.Type == "physical")
n_distinct(c(human_phys$Entrez.Gene.Interactor.A, human_phys$Entrez.Gene.Interactor.B))

# In the human set, what experimental systems used?
human_phys %>%
  group_by(Experimental.System, Throughput) %>%
  summarise(count=n()) %>%
  arrange(-count)

# If I take only low throughput, how many data points left?
human_phys_low <- filter(human, Throughput == "Low Throughput")
n_distinct(c(human_phys_low$Entrez.Gene.Interactor.A, human_phys_low$Entrez.Gene.Interactor.B))

# Let's use both low and high throughput data because yolo
human_phys <- human_phys %>% select(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B) %>% distinct()
n_distinct(c(human_phys$Entrez.Gene.Interactor.A, human_phys$Entrez.Gene.Interactor.B))
head(human_phys); nrow(human_phys)

write.table(human_phys, "~/Projects/fusion_ML/data/network_data/biogrid_human_phys.csv", sep=",", quote = FALSE, row.names = FALSE)


### OLD: had intended to use BioGRID multi-validated physical PPIs, but the dataset is too small
### Only 8072 genes (Entrez) have interaction data available

biogrid_phys <- read.csv("~/Projects/fusion_ML/data/network_data/BIOGRID-MV-Physical-3.4.135.tab2.txt", sep="\t")
head(biogrid_phys); nrow(biogrid_phys)

human <- filter(biogrid_phys, Organism.Interactor.A == 9606 & Organism.Interactor.B == 9606)
head(human); nrow(human)

# How many genes are covered now?
n_distinct(c(biogrid_phys$Entrez.Gene.Interactor.A, biogrid_phys$Entrez.Gene.Interactor.B))

# How many organisms + interactions in dataset?
biogrid_phys %>%
  group_by(Organism.Interactor.A) %>%
  summarise(count=n()) %>%
  arrange(-count)

# In the human set, what experimental systems used?
human %>%
  group_by(Experimental.System, Experimental.System.Type) %>%
  summarise(count=n()) %>%
  arrange(-count)

# Looks good. Deduplicate and export with Entrez symbols
human <- human %>% select(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B) %>% distinct()
n_distinct(c(human$Entrez.Gene.Interactor.A, human$Entrez.Gene.Interactor.B))

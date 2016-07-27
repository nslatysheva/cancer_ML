###### Confirming PPI trends using Wang data
###### 

require(igraph)
require(dplyr)

# Fetch biogrid binary interaction data as entrez symbols
d <- read.csv("~/Projects/fusion_ML/data/network_data/biogrid_human_phys.csv", sep=","); head(d)

# convert Entrez to ENSG
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)

d <- 
  merge(merge(d, entrez_ensg, by.x="Entrez.Gene.Interactor.A", by.y="entrez", all.x = TRUE),
           entrez_ensg, by.x="Entrez.Gene.Interactor.B", by.y="entrez", all.x = TRUE) %>%
  rename(entrez_a = Entrez.Gene.Interactor.A,
         entrez_b = Entrez.Gene.Interactor.B,
         ensg_a = ensg.x,
         ensg_b = ensg.y) %>% 
  select(entrez_a, entrez_b, ensg_a, ensg_b); head(d)

# Anything we can do to assign ENSGs to the remaining genes? Down from 743 to 435 NAs
unmapped_entrez <- unique(c((filter(d, is.na(ensg_a)) %>% select(entrez_a))$entrez_a, 
                            (filter(d, is.na(ensg_b)) %>% select(entrez_b))$entrez_b)) %>% length()
unmapped_entrez %>% length()
write.table(unmapped_entrez, "~/Projects/fusion_ML/data/network_data/unmapped_entrez.csv",  sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# omit any remaining rows with NA values
d <- na.omit(d); d <- d %>% select(ensg_a, ensg_b) %>% distinct(); head(d)

# convert to undirected graph object
g <- graph.data.frame(d, directed=FALSE)
V(g); E(g)

### degree centrality
degree_score <- as.data.frame(degree(g))
degree_score <- cbind(rownames(degree_score), degree_score)
colnames(degree_score) <- c("ensg","degree_centrality"); head(degree_score)

### betweenness centrality 
betweenness_score <- as.data.frame(betweenness(g, normalize = TRUE))
betweenness_score <- cbind(rownames(betweenness_score), betweenness_score)
colnames(betweenness_score) <- c("ensg","betweenness_centrality"); head(betweenness_score)

### closeness centrality
closeness_score <- as.data.frame(closeness(g, normalized = TRUE))
closeness_score <- cbind(rownames(closeness_score), closeness_score)
colnames(closeness_score) <- c("ensg","closeness_centrality"); head(closeness_score)

# ### authority centrality
# authority_score <- as.data.frame(authority.score(g))
# authority_score <- select(authority_score, 1)
# authority_score <- cbind(rownames(authority_score), authority_score); 
# colnames(authority_score) <- c("ensg","authority_centrality"); head(authority_score)
# 
# ### hub score centrality
# hub_score <- as.data.frame(hub.score(g))
# hub_score <- select(hub_score, 1)
# hub_score <- cbind(rownames(hub_score), hub_score)
# colnames(hub_score) <- c("ensg","hub_score"); head(hub_score)

### eigenvector centrality
eigenvector_centralities <- as.data.frame(eigen_centrality(g))
eigenvector_centralities <- select(eigenvector_centralities, 1)
eigenvector_centralities <- cbind(rownames(eigenvector_centralities), eigenvector_centralities)
colnames(eigenvector_centralities) <- c("ensg","eigenvector_centrality"); head(eigenvector_centralities)

# write out centrality score table
centrality <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                     list(degree_score, 
                          betweenness_score, 
                          closeness_score,
                          # authority_score,
                          # hub_score,
                          eigenvector_centralities))



write.table(centrality, "~/Projects/fusion_ML/features/centrality_biogrid.csv", sep=",", quote = FALSE, row.names = FALSE)


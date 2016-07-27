# pca and clustering
# both on fusion parents and fusion partnerings
require(dplyr)
require(tidyr)
require(FactoMineR)
require(ggplot2)

# themeset
theme_set(theme_bw(base_size = 20))

# load features
feature_space <- read.csv("~/Projects/fusion_ML/features/feature_spaces/feature_space_imputation_completed.csv", header=TRUE)
head(feature_space); dim(feature_space)

######### 1) fusion parents ######### 
##narrowing down features
# only parent genes
fusions <- feature_space %>% filter(is_parent==1); nrow(fusions)

# get rid of other parent labels
fusions <- fusions %>% select(-c(is_parent, is_5_prime_parent, is_3_prime_parent, is_recurrent_parent, num_fusions))

# how many n are set to 0
is_zero <- function(data) { sum(data == 0) }
counts_zero <- as.data.frame(fusions %>% summarize_each(funs(is_zero))/nrow(fusions)); counts_zero
percent_zero <- counts_zero %>% gather(key = var, value = percent_zero) %>% arrange(-percent_zero); percent_zero
vars_mostly_zero <- percent_zero %>% filter(percent_zero > 0.90) %>% select(var); vars_mostly_zero

# filter columns so that only 90 or more are filled
to_cluster <- fusions[, !(names(fusions) %in% vars_mostly_zero$var)]

# remove outliers (once these are determined! recursive)
outlier_rows; length(outlier_rows)
to_cluster <- to_cluster[rownames(to_cluster)[!(rownames(to_cluster) %in% outlier_rows)], ]
dim(to_cluster)

# as matrix and scale
m <- as.matrix(to_cluster[, 2:ncol(to_cluster)])
m2 <- scale(m, center = TRUE, scale = TRUE)

# factominer
res.pca <- PCA(m2, scale.unit=FALSE, ncp=10, graph = FALSE)

# plot eigenvalues
eigenvalues <- as.data.frame(cbind(rownames(res.pca$eig), seq(1:nrow(res.pca$eig)), res.pca$eig)) %>% select(-1); head(eigenvalues)
colnames(eigenvalues) <- c("component", "eigenvalue", "percentage_of_variance", "var explained")
plot_eigenvalues <- eigenvalues %>% select(-percentage_of_variance) %>% gather(key = variable, value = values, -component)
  
# plot
p <- ggplot(data = plot_eigenvalues, aes(x=component, y=values, colour=variable)) + 
  geom_line(size=0.5) + facet_grid(variable ~ ., scales = "free") + theme(legend.position="none") +
  scale_colour_manual(values = c("blue", "#cccccc")) + ylab("") + xlab(""); p
ggsave(p, filename = "~/Projects/fusion_ML/analyses/PCA/parent_eigenvalues_variance_explained.pdf", width = 4, height = 4)

######### Examine loadings on principal components
loadings <- as.data.frame(cbind(rownames(res.pca$var$cor), res.pca$var$cor))
loadings <- loadings %>% mutate_each(funs(as.character), Dim.1:Dim.6) %>% mutate_each(funs(as.numeric), Dim.1:Dim.6)
head(loadings)
# First component
first <- loadings %>% filter(abs(as.numeric(Dim.1)) >= 0.5) %>% arrange(-Dim.1) %>% select(V1, Dim.1); first
# Second component
second <- loadings %>% filter(abs(as.numeric(Dim.2)) >= 0.5) %>% arrange(-Dim.2) %>% select(V1, Dim.2); second
# Third component
third <- loadings %>% filter(abs(as.numeric(Dim.3)) >= 0.5) %>% arrange(-Dim.3) %>% select(V1, Dim.3); third
# Fourth component
fourth <- loadings %>% filter(abs(as.numeric(Dim.4)) >= 0.5) %>% arrange(-Dim.4) %>% select(V1, Dim.4); fourth
# Fifth component
fifth <- loadings %>% filter(abs(as.numeric(Dim.5)) >= 0.4) %>% arrange(-Dim.5) %>% select(V1, Dim.5); fifth
# Sixth component
sixth <- loadings %>% filter(abs(as.numeric(Dim.6)) >= 0.4) %>% arrange(-Dim.6) %>% select(V1, Dim.6); sixth

# plot loadings
first_plot <- first %>% top_n(n = 10, wt = Dim.1) %>% ggplot(aes(x=reorder(V1, Dim.1), y=Dim.1)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); first_plot
second_plot <- second %>% top_n(n = 10, wt = Dim.2) %>% ggplot(aes(x=reorder(V1, Dim.2), y=Dim.2)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); second_plot
third_plot <- third %>% top_n(n = 10, wt = Dim.3) %>% ggplot(aes(x=reorder(V1, Dim.3), y=Dim.3)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); third_plot
fourth_plot <- fourth %>% top_n(n = 10, wt = Dim.4) %>% ggplot(aes(x=reorder(V1, Dim.4), y=Dim.4)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); fourth_plot
fifth_plot <- fifth %>% top_n(n = 10, wt = Dim.5) %>% ggplot(aes(x=reorder(V1, Dim.5), y=Dim.5)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); fifth_plot
sixth_plot <- sixth %>% top_n(n = 10, wt = Dim.6) %>% ggplot(aes(x=reorder(V1, Dim.6), y=Dim.6)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); sixth_plot

ggsave(first_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_first_component_loadings.pdf", height=7, width=5.5)
ggsave(second_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_second_component_loadings.pdf", height=3.3, width=5)
ggsave(third_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_third_component_loadings.pdf", height=2.5, width=4)
ggsave(fourth_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_fourth_component_loadings.pdf", height=1.5, width=5.8)
ggsave(fifth_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_fifth_component_loadings.pdf", height=2.5, width=6)
ggsave(sixth_plot, file = "~/Projects/fusion_ML/analyses/PCA/PARENTS_sixth_component_loadings.pdf", height=2.5, width=5)


######### Clustering on principal components
hcpc <- HCPC(res.pca, nb.clust=10, conso=0, min=3, max=10, graph = FALSE)
hcpc <- HCPC(res.pca, nb.clust=-1, conso=0, min=3, max=10, graph = FALSE)

# examine clusters
hcpc$data.clust
# Description of the clusters by the variables
hcpc$desc.var
# Description of the clusters by the factors (axes)
hcpc$desc.axes

# Nicer tree plot
plot(hcpc, choice ="tree", draw.tree = TRUE,
     ind.names = FALSE, title="", tree.barplot = FALSE, rect = TRUE)

# Inertia gain bar plot
plot(hcpc, choice ="bar", title="")

# Nicer factor map
plot(hcpc, choice ="map", draw.tree = FALSE,
     ind.names = FALSE, centers.plot = TRUE, title="")

# examine clusters
# closest individuals to barycentre of cluster
cluster_counts <- hcpc$data.clust %>% group_by(clust) %>% summarise(count =n()); cluster_counts
cluster_counts$clust <- as.numeric(cluster_counts$clust)
too_small_clusters <- (cluster_counts %>% filter(count <10))$clust; too_small_clusters
# exclude data points in these tiny clusters
outlier_rows <- rownames(hcpc$data.clust[hcpc$data.clust$clust %in% too_small_clusters,]); outlier_rows

# plot cluster size
cluster_counts$cluster <- paste("cluster ", cluster_counts$clust, sep="")
p <- ggplot(data = cluster_counts, aes(x=reorder(cluster, -clust), y=count)) + geom_bar(stat="identity", aes(fill=cluster)) +
  coord_flip() + xlab("") + scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) +
  theme(legend.position="none"); p
ggsave(p, file="~/Projects/fusion_ML/analyses/PCA/PARENTS_cluster_counts.pdf", height=2.5, width=3.5)


########## Plot features
source(file = "~/Projects/fusion_ML/scripts/analysis/multiplot_function.R") # multiplot ggplots

# get back original dataset before scaling
clustered_data <- cbind(to_cluster, hcpc$data.clust %>% select(clust))
# write out dataset
write.csv(clustered_data, "~/Projects/fusion_ML/Tables/Table S6 - PCA and clustering tables/parent_protein_clustering.csv", quote = FALSE, row.names = FALSE)

dim(clustered_data); head(clustered_data)
head(first); head(second); head(third); 

p1 <- ggplot(data = clustered_data, aes(x=clust, y=num_cancer_interactions)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,65)); p1

p2 <- ggplot(data = clustered_data, aes(x=clust, y=num_ANCHOR_LMs)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,150)); p2

p3 <- ggplot(data = clustered_data, aes(x=clust, y=avg_ensp_length)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,1700)); p3

head(fourth); head(fifth); head(sixth)
clustered_data$betweenness_centrality
p4 <- ggplot(data = clustered_data, aes(x=clust, y=betweenness_centrality)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = FALSE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,0.0025)); p4

p5 <- ggplot(data = clustered_data, aes(x=clust, y=density_INstruct_domains)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = FALSE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,0.02)); p5

proportions <- clustered_data %>% group_by(clust) %>% summarise(proportion = sum(loc_nucleus)/n())
p6 <- ggplot(data = proportions, aes(x=clust, y=proportion)) +
  geom_bar(stat="identity", aes(fill = clust), colour="#404040", size=0.9, width=0.5) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F")) + xlab("") + ylab("% loc_nucleus") +
  guides(fill=FALSE)  + coord_cartesian(ylim = c(0,0.85)); p6

# p6 <- ggplot(data = clustered_data, aes(x=clust, y=prot_half_life)) +
# geom_violin(aes(fill=is_recurrent_parent), alpha=0.2, color="lightgrey") +
#   geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#404040", size=1.2, width=0.7) +
#   scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
#   guides(fill=FALSE) + coord_cartesian(ylim = c(0,00)); p6

# good size is 7.6x12, multiplot_PARENTS_features
multiplot(p1, p4, p2, p5, p3, p6, cols=3)

# Description of the clusters by the variables
# cluster 1
head(data.frame(rownames(hcpc$desc.var$quanti$`1`), 
                hcpc$desc.var$quanti$`1`)  %>% arrange(p.value), 20)
# cluster 2
head(data.frame(rownames(hcpc$desc.var$quanti$`2`), 
                hcpc$desc.var$quanti$`2`)  %>% arrange(p.value), 20)
# cluster 3
head(data.frame(rownames(hcpc$desc.var$quanti$`3`), 
                hcpc$desc.var$quanti$`3`)  %>% arrange(p.value), 10)

colnames(clustered_data)
ggplot(data = clustered_data, aes(x=clust, y=num_OG_interactions)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = FALSE, aes(fill = clust), color="#404040", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE)



###### Do fusions occur most often between parents in the same class?
# read in fusion list and ensg conversion
fusions <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/klijn_gene_fusions.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
fusions <- fusions %>% select(id, entrez_id_3, entrez_id_5) %>% distinct; head(fusions); nrow(fusions)
# convert fuson entrez to ensg
ensg_entrez <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
# label conversion
fusions_ensg <- merge(merge(fusions, ensg_entrez, by.x = "entrez_id_5", by.y = "entrez", all.x = TRUE), 
                      ensg_entrez, by.x = "entrez_id_3", by.y = "entrez", all.x = TRUE) %>% 
  rename(ensg_3 = ensg.y, ensg_5 = ensg.x, fusion_id = id) %>% distinct() 
head(fusions_ensg)

# get parents with cluster labels
parents_with_clust <- clustered_data %>% select(ensg, clust); head(parents_with_clust)

# overlap onto fusion event
fusions_with_clust <- merge(merge(fusions_ensg, parents_with_clust, by.x = "ensg_5", by.y = "ensg", all.x = TRUE), 
                            parents_with_clust, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>% 
  rename(clust_3 = clust.y, clust_5 = clust.x) %>% 
  distinct %>%
  select(fusion_id, ensg_5, ensg_3, clust_5, clust_3) %>% 
  filter(!is.na(clust_3) & !is.na(clust_5))

head(fusions_with_clust)

# summarise frequencies
n_fusions_in_dataset =  fusions_with_clust %>% select(fusion_id) %>% nrow
fusion_freqs_by_cluster <- fusions_with_clust %>% 
  group_by(clust_5, clust_3) %>% 
  summarise(n_fusions = n_distinct(fusion_id)) %>% 
  ungroup %>%
  mutate(percent = 100*n_fusions/n_fusions_in_dataset); fusion_freqs_by_cluster

# test
chi <- chisq.test(matrix(fusion_freqs_by_cluster$n_fusions, nrow = 3, byrow = TRUE)); chi
chi$expected

# 1/2 fusion should count as a 2/1 fusion, etc
fusion_freqs_by_cluster_dedup <- fusion_freqs_by_cluster %>% 
  mutate_each(funs(as.character), clust_5, clust_3) %>%
  mutate(cluster_A = pmin(clust_5, clust_3), 
         cluster_B = pmax(clust_5, clust_3)) %>%
  group_by(cluster_A, cluster_B) %>%
  summarise(total_n_fusions = sum(n_fusions)) %>%
  ungroup %>%
  mutate(percent = 100*total_n_fusions/n_fusions_in_dataset)

fusion_freqs_by_cluster_dedup

# just summarise which clusters act as 5' and 3' partners
head(fusions_with_clust %>% arrange(fusion_id))
five_summ <- summary(as.factor(fusions_with_clust$clust_5)); five_summ
three_summ <- summary(as.factor(fusions_with_clust$clust_3)); three_summ
chi <- chisq.test(matrix(c(five_summ, three_summ), nrow=2, byrow=TRUE)); chi
chi$expected

# summarise same vs. different cluster fusion frequencies
percent_same_cluster_observed <- fusion_freqs_by_cluster %>% filter(clust_5 == clust_3) %>% 
  summarise(percent_same_cluster = sum(percent)); percent_same_cluster_observed
fusion_freqs_by_cluster

# get expected proportion of same cluster fusions
get_random_proportions_clusters <- function(n_rand_fusions = 10000){
  cluster_proportions <- matrix(NA, nrow = n_rand_fusions, ncol = 10)
  
  for(i in 1:n_rand_fusions){
    # get random fusoions
    random_5_clust <- sample(1:nrow(parents_with_clust), 1000, replace = TRUE); parents_with_clust$clust[random_5_clust]
    random_3_clust <- sample(1:nrow(fusions_with_clust), 1000, replace = TRUE); parents_with_clust$clust[random_3_clust]
    random_fusion_clusters <- data.frame(parents_with_clust$clust[random_5_clust], 
                                         parents_with_clust$clust[random_3_clust]); names(random_fusion_clusters) <- c("clust_5", "clust_3")
    # get all proportions
    random_fusion_clusters_summary <- random_fusion_clusters %>% 
      group_by(clust_5, clust_3) %>% 
      summarise(n_fusions = n()) %>% 
      ungroup %>%
      mutate(percent = 100*n_fusions/1000); random_fusion_clusters_summary
    
    # count proportion same fusion
    percent_same_cluster <- random_fusion_clusters_summary %>% filter(clust_5 == clust_3) %>% 
      summarise(percent_same_cluster = sum(percent)); percent_same_cluster
    
    # output
    cluster_proportions[i,] <- c(random_fusion_clusters_summary$percent, percent_same_cluster$percent_same_cluster)
    
  }
  # format names
  cols <- c(paste("prop_", random_fusion_clusters_summary$clust_5, "_", random_fusion_clusters_summary$clust_3, sep=""),
            "prop_same_clust"); cols
  cluster_proportions <- data.frame(cluster_proportions); names(cluster_proportions) <- cols
  return(cluster_proportions)
}

random_proportions <- get_random_proportions_clusters(n_rand_fusions = 10000); head(random_proportions)
summary(random_proportions$prop_same_clust)
summary(random_proportions$prop_1_1)
str(random_proportions)


# set up dataframe for plotting (dark grey in observed)
head(random_proportions)

# plot
# same clust
p <- ggplot(data = random_proportions, aes(prop_same_clust)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = percent_same_cluster_observed$percent_same_cluster, size=2, color="gold") +
  annotate("text", x=37, y=150,label=length(random_proportions$prop_same_clust[random_proportions$prop_same_clust<percent_same_cluster_observed$percent_same_cluster])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_same_cluster.pdf", height = 4, width = 6)     

# other plots
observed = (fusion_freqs_by_cluster %>% filter(clust_5==1 & clust_3==1))$percent
p <- ggplot(data = random_proportions, aes(prop_1_1)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") + 
  annotate("text", x=24.6, y=20,label=length(random_proportions$prop_1_1[random_proportions$prop_1_1<observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_1_1.pdf",height = 4, width = 6)    

observed = (fusion_freqs_by_cluster %>% filter(clust_5==1 & clust_3==2))$percent; observed
p <- ggplot(data = random_proportions, aes(prop_1_2)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") +
  annotate("text", x=13, y=30,label=length(random_proportions$prop_1_2[random_proportions$prop_1_2<observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_1_2.pdf", height = 4, width = 6)  

observed = (fusion_freqs_by_cluster %>% filter(clust_5==1 & clust_3==3))$percent
p <- ggplot(data = random_proportions, aes(prop_1_3)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") +
  annotate("text", x=11, y=300,label=length(random_proportions$prop_1_3[random_proportions$prop_1_3>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_1_3.pdf", height = 4, width = 6)  

# clust 2s
observed = (fusion_freqs_by_cluster %>% filter(clust_5==2 & clust_3==1))$percent
p <- ggplot(data = random_proportions, aes(prop_2_1)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") +
  annotate("text", x=18.5, y=50,label=length(random_proportions$prop_2_1[random_proportions$prop_2_1>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_2_1.pdf", height = 4, width = 6)  

observed = (fusion_freqs_by_cluster %>% filter(clust_5==2 & clust_3==2))$percent
p <- ggplot(data = random_proportions, aes(prop_2_2)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") +
  annotate("text", x=10.5, y=300,label=length(random_proportions$prop_2_2[random_proportions$prop_2_2>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_2_2.pdf", height = 4, width = 5)  

observed = (fusion_freqs_by_cluster %>% filter(clust_5==2 & clust_3==3))$percent
p <- ggplot(data = random_proportions, aes(prop_2_3)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") + 
  annotate("text", x=6.5, y=150,label=length(random_proportions$prop_2_3[random_proportions$prop_2_3>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_2_3.pdf",height = 4, width = 6)  

# clust 3s
observed = (fusion_freqs_by_cluster %>% filter(clust_5==3 & clust_3==1))$percent
p <- ggplot(data = random_proportions, aes(prop_3_1)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") +
  annotate("text", x=10.5, y=300,label=length(random_proportions$prop_3_1[random_proportions$prop_3_1>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_3_1.pdf", height = 4, width = 6)  

observed = (fusion_freqs_by_cluster %>% filter(clust_5==3 & clust_3==2))$percent
p <- ggplot(data = random_proportions, aes(prop_3_2)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") + 
  annotate("text", x=7, y=100,label=length(random_proportions$prop_3_2[random_proportions$prop_3_2>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_3_2.pdf", height = 4, width = 6)    

observed = (fusion_freqs_by_cluster %>% filter(clust_5==3 & clust_3==3))$percent
p <- ggplot(data = random_proportions, aes(prop_3_3)) + 
  geom_histogram(binwidth = 0.1, fill="grey") + xlab("") + 
  geom_vline(xintercept = observed, size=2, color="gold") + 
  annotate("text", x=5.5, y=40,label=length(random_proportions$prop_3_3[random_proportions$prop_3_3>observed])/10000, size=5, color="darkgrey"); p
ggsave(p, filename = "~/Projects/fusion_ML/figures/Figure 5 - Clustering fusion events/random_proportion_3_3.pdf", height = 4, width = 6)  


##### get paragon examples outs
# paragons, most typical individuals in each cluster
paragons1_rows <- data.frame(hcpc$desc.ind$para$`1`); paragons1_rows
paragons2_rows <- data.frame(hcpc$desc.ind$para$`2`); paragons2_rows
paragons3_rows <- data.frame(hcpc$desc.ind$para$`3`); paragons3_rows

dim(hcpc$data.clust); dim(m2)
head(clustered_data)
# get paragons
paragons1 <- clustered_data[(rownames(clustered_data) %in% rownames(paragons1_rows)), c(1:2, 51)]; paragons1; paragons1_rows
paragons2 <- clustered_data[(rownames(clustered_data) %in% rownames(paragons2_rows)), c(1:2, 51)]; paragons2; paragons2_rows
paragons3 <- clustered_data[(rownames(clustered_data) %in% rownames(paragons3_rows)), c(1:2, 51)]; paragons3; paragons3_rows

dim(to_cluster)
dim(clustered_data)

# ENSG00000187164, ENSG00000160007, ENSG00000077684, ENSG00000047188, ENSG00000161847
# ENSG00000119335, ENSG00000117676, ENSG00000185651, ENSG00000126581, ENSG00000134058



######### 2) fusion pairings ######### 
fusion_pairs <- read.csv("~/Projects/fusion_ML/features/feature_spaces/fusion_events_with_features_original.csv", header=TRUE)
head(fusion_pairs)
names(fusion_pairs)

# how many n are set to 0
is_zero <- function(data) { sum(data == 0) }
counts_zero <- as.data.frame(fusion_pairs %>% summarize_each(funs(is_zero))/nrow(fusion_pairs)); counts_zero
percent_zero <- counts_zero %>% gather(key = var, value = percent_zero) %>% arrange(-percent_zero); percent_zero
vars_mostly_zero <- percent_zero %>% filter(percent_zero > 0.90) %>% select(var); vars_mostly_zero

# filter columns so that only relatively filled ones remain
to_cluster <- fusion_pairs[, !(names(fusion_pairs) %in% vars_mostly_zero$var)]; head(fusion_pairs)

# remove outliers (once these are determined! recursive)
outlier_rows; length(outlier_rows)
to_cluster <- to_cluster[rownames(to_cluster)[!(rownames(to_cluster) %in% outlier_rows)], ]
dim(to_cluster)

# create matrix, scale
m <- as.matrix(to_cluster %>% select(-c(ensg_3, ensg_5, fusion_id, cancer_type)))
m2 <- scale(m, center = TRUE, scale = TRUE)

# factominer
res.pca <- PCA(m2, scale.unit=FALSE, ncp=10, graph = FALSE)

# plot eigenvalues
eigenvalues <- as.data.frame(cbind(rownames(res.pca$eig), seq(1:nrow(res.pca$eig)), res.pca$eig)) %>% select(-1)
colnames(eigenvalues) <- c("component", "eigenvalue", "percentage_of_variance", "var explained")
plot_eigenvalues <- eigenvalues %>% select(-percentage_of_variance) %>% gather(key = variable, value = values, -component)

# plot
p <- ggplot(data = plot_eigenvalues, aes(x=component, y=values, colour=variable)) + 
  geom_line(size=1) + facet_grid(variable ~ ., scales = "free") + theme(legend.position="none") +
  scale_colour_manual(values = c("blue", "#cccccc")) + ylab("") + xlab(""); p
ggsave(p, filename = "~/Projects/fusion_ML/analyses/PCA/fusion_pairs_eigenvalues_variance_explained.pdf", width = 4, height = 4)

# Examine loadings on principal components
loadings <- as.data.frame(cbind(rownames(res.pca$var$cor), res.pca$var$cor))
loadings <- loadings %>% mutate_each(funs(as.character), Dim.1:Dim.6) %>% mutate_each(funs(as.numeric), Dim.1:Dim.6)

# First component
first <- loadings %>% filter(abs(as.numeric(Dim.1)) >= 0.5) %>% arrange(-Dim.1) %>% select(V1, Dim.1); first
# Second component
second <- loadings %>% filter(abs(as.numeric(Dim.2)) >= 0.5) %>% arrange(-Dim.2) %>% select(V1, Dim.2); second
# Third component
third <- loadings %>% filter(abs(as.numeric(Dim.3)) >= 0.5) %>% arrange(-Dim.3) %>% select(V1, Dim.3); third
# Fourth component
fourth <- loadings %>% filter(abs(as.numeric(Dim.4)) >= 0.5) %>% arrange(-Dim.4) %>% select(V1, Dim.4); fourth
# Fifth component
fifth <- loadings %>% filter(abs(as.numeric(Dim.5)) >= 0.4) %>% arrange(-Dim.5) %>% select(V1, Dim.5); fifth
# Sixth component
sixth <- loadings %>% filter(abs(as.numeric(Dim.6)) >= 0.4) %>% arrange(-Dim.6) %>% select(V1, Dim.6); sixth

# plot
first_plot <- first %>% top_n(n = 10, wt = Dim.1) %>% ggplot(aes(x=reorder(V1, Dim.1), y=Dim.1)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); first_plot
second_plot <- second %>% top_n(n = 10, wt = Dim.2) %>% ggplot(aes(x=reorder(V1, Dim.2), y=Dim.2)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); second_plot
third_plot <- third %>% top_n(n = 10, wt = Dim.3) %>% ggplot(aes(x=reorder(V1, Dim.3), y=Dim.3)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); third_plot
fourth_plot <- fourth %>% top_n(n = 10, wt = Dim.4) %>% ggplot(aes(x=reorder(V1, Dim.4), y=Dim.4)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); fourth_plot
fifth_plot <- fifth %>% top_n(n = 10, wt = Dim.5) %>% ggplot(aes(x=reorder(V1, Dim.5), y=Dim.5)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); fifth_plot
sixth_plot <- sixth %>% top_n(n = 10, wt = Dim.6) %>% ggplot(aes(x=reorder(V1, Dim.6), y=Dim.6)) + geom_bar(stat="identity", fill="lightgrey", colour = "#404040", width = 0.8) + coord_flip() + labs(x="", y=""); sixth_plot

ggsave(first_plot, file = "~/Projects/fusion_ML/analyses/PCA/first_component_loadings.pdf", height=6.2, width=5.5)
ggsave(second_plot, file = "~/Projects/fusion_ML/analyses/PCA/second_component_loadings.pdf", height=3.5, width=5.5)
ggsave(third_plot, file = "~/Projects/fusion_ML/analyses/PCA/third_component_loadings.pdf", height=3, width=5.5)
ggsave(fourth_plot, file = "~/Projects/fusion_ML/analyses/PCA/fourth_component_loadings.pdf", height=2, width=5.5)
ggsave(fifth_plot, file = "~/Projects/fusion_ML/analyses/PCA/fifth_component_loadings.pdf", height=2.5, width=5.8)
ggsave(sixth_plot, file = "~/Projects/fusion_ML/analyses/PCA/sixth_component_loadings.pdf", height=2.5, width=5.5)

# graph other stuff
# res.pca <- PCA(m2, scale.unit=TRUE, ncp=10, graph = TRUE)

# Clustering on principal components
hcpc_fusion_pairs <- HCPC(res.pca, nb.clust=10, conso=0, min=3, max=10, graph = FALSE)
hcpc_fusion_pairs <- HCPC(res.pca, nb.clust=-1, conso=0, min=3, max=10, graph = FALSE)

# Nicer tree plot, 6x6
plot(hcpc_fusion_pairs, choice ="tree", draw.tree = TRUE,
     ind.names = FALSE, title="", tree.barplot = FALSE, rect = TRUE)

# Inertia gain bar plot, 15x8
plot(hcpc_fusion_pairs, choice ="bar", title="")

# Nicer factor map
plot(hcpc_fusion_pairs, choice ="map", draw.tree = FALSE,
     ind.names = FALSE, centers.plot = TRUE, title="")

# examine clusters
hcpc_fusion_pairs$data.clust

# examine clusters
# closest individuals to barycentre of cluster
cluster_counts <- hcpc_fusion_pairs$data.clust %>% group_by(clust) %>% summarise(count =n()); cluster_counts
cluster_counts$clust <- as.numeric(cluster_counts$clust)
too_small_clusters <- (cluster_counts %>% filter(count <10))$clust; too_small_clusters
# exclude data points in these tiny clusters
outlier_rows <- rownames(hcpc_fusion_pairs$data.clust[hcpc_fusion_pairs$data.clust$clust %in% too_small_clusters,]); outlier_rows

# plot cluster size
cluster_counts$cluster <- paste("cluster ", cluster_counts$clust, sep="")
p <- ggplot(data = cluster_counts, aes(x=reorder(cluster, -clust), y=count)) + geom_bar(stat="identity", aes(fill=cluster)) +
  coord_flip() + xlab("") + scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) +
  theme(legend.position="none"); p
ggsave(p, file="~/Projects/fusion_ML/analyses/PCA/PARENTS_cluster_counts.pdf", height=2.5, width=3.5)



# Description of the clusters by the variables
# cluster 1
head(data.frame(rownames(hcpc_fusion_pairs$desc.var$quanti$`1`), 
                hcpc_fusion_pairs$desc.var$quanti$`1`)  %>% arrange(p.value), 10)
# cluster 2
head(data.frame(rownames(hcpc_fusion_pairs$desc.var$quanti$`2`), 
                hcpc_fusion_pairs$desc.var$quanti$`2`)  %>% arrange(p.value), 20)
# cluster 3
head(data.frame(rownames(hcpc_fusion_pairs$desc.var$quanti$`3`), 
                hcpc_fusion_pairs$desc.var$quanti$`3`)  %>% arrange(p.value), 40)

# Description of the clusters by the factors (axes)
hcpc_fusion_pairs$desc.axes

# make some distribution plots of differneces between the three clusters of parents
source(file = "~/Projects/fusion_ML/scripts/analysis/multiplot_function.R") # multiplot ggplots

# get back original dataset before scaling
clustered_data <- cbind(to_cluster, hcpc_fusion_pairs$data.clust %>% select(clust))
# print out clustered fusion event data
write.csv(clustered_data, "~/Projects/fusion_ML/Tables/Table S6 - PCA and clustering tables/fusion_event_clustering.csv", quote = FALSE, row.names = FALSE)

dim(clustered_data); dim(to_cluster)
head(first); head(second); head(third); 

# PC1
p1 <- ggplot(data = clustered_data, aes(x=clust, y=num_OG_interactions_3_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = FALSE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,16)); p1

p2 <- ggplot(data = clustered_data, aes(x=clust, y=density_PTMs_3_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,0.2)); p2

# PC 2
p3 <- ggplot(data = clustered_data, aes(x=clust, y=num_cancer_interactions_5_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,70)); p3

p4 <- ggplot(data = clustered_data, aes(x=clust, y=degree_centrality_5_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,410)); p4

# PC3
p5 <- ggplot(data = clustered_data, aes(x=clust, y=num_ANCHOR_LMs_5_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,130)); p5

# other
head(fourth); head(fifth); head(sixth)

p6 <- ggplot(data = clustered_data, aes(x=clust, y=avg_ensp_length_5_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=0.9, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,1600)); p6

p7 <- ggplot(data = clustered_data, aes(x=clust, y=mean_expression_5_prime)) +
  geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=1.2, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,12)); p7

p8 <- ggplot(data = clustered_data, aes(x=clust, y=mean_expression_3_prime)) +
geom_violin(aes(fill=clust), alpha=0.2, color="lightgrey") +
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = clust), color="#808080", size=1.2, width=0.7) +
  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F", "#3A55A2", "#7DD0DE", "#B7559E")) + xlab("") +
  guides(fill=FALSE) + coord_cartesian(ylim = c(0,12)); p8

# proportions <- clustered_data %>% group_by(clust) %>% summarise(proportion = sum(loc_nucleus)/n())
# p6 <- ggplot(data = proportions, aes(x=clust, y=proportion)) +
#   geom_bar(stat="identity", aes(fill = clust), colour="#808080", size=0.9, width=0.5) +
#   scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F")) + xlab("") + ylab("% loc_nucleus") +
#   guides(fill=FALSE)  + coord_cartesian(ylim = c(0,0.85)); p6


# good size is 8X13, multiplot_FUSION_EVENTS
multiplot(p1, p5, p2, p6, p3, p7, p4, p8, cols=4)


##### get paragon examples out
# paragons, most typical individuals in each cluster
paragons1_rows <- data.frame(hcpc_fusion_pairs$desc.ind$para$`1`); paragons1_rows
paragons2_rows <- data.frame(hcpc_fusion_pairs$desc.ind$para$`2`); paragons2_rows
paragons3_rows <- data.frame(hcpc_fusion_pairs$desc.ind$para$`3`); paragons3_rows

dim(hcpc_fusion_pairs$data.clust); dim(m2)
head(clustered_data); dim(clustered_data)

# get gene symbols onto clustered data
gene_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_gene_symbol_conversion.txt", header=TRUE, sep="\t")
clustered_data_genes <- merge(merge(clustered_data, gene_ensg, by.x = "ensg_5", by.y = "ensg", all.x = TRUE),
                              gene_ensg, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>%
  rename(gene_5 = name.x, 
         gene_3 = name.y) %>%
  select(ensg_5, gene_5, ensg_3, gene_3, fusion_id:clust)

# get paragons
head(clustered_data_genes)
paragons1 <- clustered_data_genes[(rownames(clustered_data_genes) %in% rownames(paragons1_rows)), 
                                  c(1:4, ncol(clustered_data_genes))]; paragons1; paragons1_rows
paragons2 <- clustered_data_genes[(rownames(clustered_data_genes) %in% rownames(paragons2_rows)), 
                                  c(1:4, ncol(clustered_data_genes))]; paragons2; paragons2_rows
paragons3 <- clustered_data_genes[(rownames(clustered_data_genes) %in% rownames(paragons3_rows)), 
                                  c(1:4, ncol(clustered_data_genes))]; paragons3; paragons3_rows

# cancer types and clusters - spread equally?
head(clustered_data_genes)
too_rare <- clustered_data_genes %>% group_by(cancer_type) %>% summarise(cancer_n = n()) %>% filter(cancer_n < 50)

clust_by_cancer <- clustered_data_genes %>% filter(!(cancer_type %in% too_rare$cancer_type)) %>%
  group_by(cancer_type, clust) %>% summarise(n = n()) %>%
  inner_join(clustered_data_genes %>% group_by(cancer_type) %>% summarise(cancer_n = n())) %>%
  mutate(percent = 100*n/cancer_n, 
         percent_clust_1 = ifelse(clust==1, percent, 0),
         percent_clust_2 = ifelse(clust==2, percent, 0),
         percent_clust_3 = ifelse(clust==3, percent, 0))

p <- ggplot(clust_by_cancer, aes(reorder(cancer_type, percent_clust_3), percent, fill=clust)) + 
  geom_bar(stat = "identity", width=0.9) +  scale_fill_manual(values=c('#525050', '#EB242E', "#47B44F")) +
  xlab("cancer type") + ylab("% of fusions") + guides(fill=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)); p
ggsave(p, file="~/Projects/fusion_ML/analyses/PCA/cancer_type_clusters.pdf", height=5, width=3.5)


# quick test
m <- matrix(as.numeric(as.character(clust_by_cancer$n)), ncol=3, byrow = TRUE); m
m <- data.frame(m); names(m) <- c("C1", "C2", "C3"); m <- m %>% filter(C1 >=20 & C2 >= 20 & C3 >= 20); m
chi <- chisq.test(m); chi
chi$observed; chi$expected
fisher.test(m)
chi$observed; chi$expected
# quick looks
clust_by_cancer %>% arrange()

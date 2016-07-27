# Pathway data
require(dplyr)
require(tidyr)


##### GET REACTOME DATA

# get pathways
pathways <- read.table("~/Projects/fusion_ML/data/pathway_data/REACTOME_from_Biomart_clean.txt", sep="\t", header=TRUE, stringsAsFactors = TRUE)
head(pathways); nrow(pathways); str(pathways)

# get pathway annotations
pathway_names <- read.table("~/Projects/fusion_ML/data/pathway_data/REACTOME_pathway_mapping_human.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)
summary(as.factor(pathway_names$Species)); pathway_names <- pathway_names %>% select(-Species) %>% distinct; head(pathway_names)
# check for no multiple annotations
multiple_descriptions <- pathway_names %>% group_by(Reactome_ID) %>% summarise(n_descriptions = n_distinct(Description)) %>% 
  arrange(-n_descriptions) %>% filter(n_descriptions > 1); multiple_descriptions
pathway_names %>% filter(Reactome_ID %in% multiple_descriptions$Reactome_ID) %>% arrange(Reactome_ID)

# merge pathway annotations onto pathways
pathways_with_annot <- merge(pathways, pathway_names, by="Reactome_ID", all.x=TRUE) %>% select(ensg, Reactome_ID, Description) %>% distinct
head(pathways_with_annot)
nrow(pathways_with_annot); nrow(na.omit(pathways_with_annot))
write.csv(pathways_with_annot, "~/Projects/fusion_ML/data/pathway_data/pathways_with_annot.csv", quote = FALSE, row.names = FALSE)

# improve pathway annotations
unannotated_REACTOME_IDs <- pathways_with_annot %>% filter(is.na(Description)) %>% select(Reactome_ID) %>% distinct; unannotated_REACTOME_IDs
unannotated_REACTOME_IDs$Reactome_ID
write.csv(unannotated_REACTOME_IDs, "~/Projects/fusion_ML/data/pathway_data/unannotated_REACTOME_IDs.csv", quote = FALSE, row.names = FALSE)
pathways_with_annot <- na.omit(pathways_with_annot)

# how many ensg have pathways?
unique(pathways_with_annot$ensg) %>% length

# how many pathways per ensg?
pathway_counts <- pathways_with_annot %>% group_by(ensg) %>% 
  summarise(n_pathways = n_distinct(Reactome_ID)) %>% arrange(-n_pathways)
summary(pathway_counts$n_pathways)

#### LOAD FUSION EVENTS WITH ENSGs

# read in fusion list and ensg conversion
fusions <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/klijn_gene_fusions.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
fusions <- fusions %>% select(id, entrez_id_3, entrez_id_5) %>% distinct; head(fusions); nrow(fusions)
# convert fuson entrez to ensg
ensg_entrez <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
# label conversion
fusions_ensg <- merge(merge(fusions, ensg_entrez, by.x = "entrez_id_5", by.y = "entrez", all.x = TRUE), 
                      ensg_entrez, by.x = "entrez_id_3", by.y = "entrez", all.x = TRUE) %>% rename(ensg_3 = ensg.y, ensg_5 = ensg.x, fusion_id = id) %>% distinct(); head(fusions_ensg)

### WHICH REACTOME PATHWAYS ARE FUSED?

fusions_with_pathways <- merge(merge(fusions_ensg, pathways_with_annot, by.x = "ensg_5", by.y = "ensg", all.x = TRUE), 
                               pathways_with_annot, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>% 
  rename(pathway_3 = Reactome_ID.y, pathway_5 = Reactome_ID.x, description_5 = Description.x, description_3 = Description.y) %>% 
  distinct %>%
  select(fusion_id, ensg_5, pathway_5, description_5, ensg_3, pathway_3, description_3) %>% 
  filter(!is.na(pathway_5) & !is.na(pathway_3))

# how many fusions are left?
length(unique(fusions_ensg$fusion_id))
length(unique(fusions_with_pathways$fusion_id))

# largest possible number of fusions between pathways is number of fusions
head(fusions_with_pathways); str(fusions_with_pathways)
max_possible_n_fusions <- length(unique(fusions_with_pathways$fusion_id)); max_possible_n_fusions

n_fusions <- fusions_with_pathways %>% group_by(pathway_5, description_5, pathway_3, description_3) %>% 
  summarise(times_fused = n_distinct(fusion_id),
            times_fused_percent = 100*times_fused/max_possible_n_fusions) %>% 
  ungroup %>% arrange(-times_fused)

# verify things
fusions_with_pathways %>% filter(description_5 == "Signal Transduction", description_3 == "Signal Transduction") %>% arrange(fusion_id)
fusions_with_pathways %>% filter(description_5 == "Immune System", description_3 == "Signal Transduction") %>% arrange(fusion_id)
fusions_with_pathways %>% filter(description_5 == "Signal Transduction", description_3 == "Immune System") %>% arrange(fusion_id)

# have to account for A,B order
n_fusions_dedup <- n_fusions %>% mutate_each(funs(as.character), pathway_5, pathway_3) %>%
  mutate(pathway_1 = pmin(pathway_5, pathway_3), 
         pathway_2 = pmax(pathway_5, pathway_3)) %>%
  group_by(pathway_1, pathway_2) %>%
  summarise(total_times_fused=sum(times_fused),
            total_percent = sum(times_fused_percent)) %>%
  ungroup %>%
  arrange(-total_times_fused); head(n_fusions_dedup)


#### GET FREQUENCIES OF PATHWAY FUSIONS BY RANDOM

# get all genes
all_ensg <- read.table("~/Projects/fusion_ML/features/ensg_basic_structural_info.csv", sep=",", header=TRUE) %>% select(ensg) %>% distinct; head(all_ensg)

# make 100,000 random fusions
create_random_fusions <- function(gene_set = all_ensg, n_random_fusions = 10000){
  
  # output matrix
  random_pairs <- matrix(NA, nrow=n_random_fusions, ncol=3)
  # generate random pairs
  for(i in 1:n_random_fusions){
    random_pair <-  c(i,
                      as.character(gene_set[sample(1:nrow(gene_set), 1),]), 
                      as.character(gene_set[sample(1:nrow(gene_set), 1),]))
    random_pairs[i,] <- random_pair
  }
  # clean and output
  random_pairs <- as.data.frame(random_pairs)
  colnames(random_pairs) <- c("fusion_id", "ensg_5", "ensg_3")
  return(random_pairs)
}
random_pairs <- create_random_fusions(n_random_fusions = 70000); head(random_pairs)

# get pathways of these things
random_fusions_with_pathways <- merge(merge(random_pairs, pathways_with_annot, by.x = "ensg_5", by.y = "ensg", all.x = TRUE), 
                               pathways_with_annot, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>% 
  rename(pathway_3 = Reactome_ID.y, pathway_5 = Reactome_ID.x, description_5 = Description.x, description_3 = Description.y) %>% 
  distinct() %>%
  select(fusion_id, ensg_5, pathway_5, description_5, ensg_3, pathway_3, description_3) %>% na.omit

head(random_fusions_with_pathways); nrow(random_fusions_with_pathways)

# how many fusions are left?
length(unique(fusions_with_pathways$fusion_id))
length(unique(random_fusions_with_pathways$fusion_id))

# get just 10,000
ten_k_fusions_indices <- sample(1:length(unique(random_fusions_with_pathways$fusion_id)), size = 10000, replace = FALSE); length(ten_k_fusions_indices)
ten_k_fusions <- unique(random_fusions_with_pathways$fusion_id)[ten_k_fusions_indices]
random_fusions_with_pathways <- random_fusions_with_pathways %>% filter(fusion_id %in% ten_k_fusions)
length(unique(random_fusions_with_pathways$fusion_id))

# count by pathway joins
max_f <- length(unique(random_fusions_with_pathways$fusion_id)); max_f

n_fusions_random <- random_fusions_with_pathways %>% group_by(pathway_5, pathway_3) %>% 
  summarise(times_fused_random = as.integer(n_distinct(fusion_id))) %>%
  mutate(times_fused_percent_random = 100*times_fused_random/max_f) %>%
  ungroup %>% arrange(-times_fused_random)

head(n_fusions_random)

# have to account for A,B order
n_fusions_random_dedup <- n_fusions_random %>% mutate_each(funs(as.character), pathway_5, pathway_3) %>%
  mutate(pathway_1 = pmin(pathway_5, pathway_3), 
         pathway_2 = pmax(pathway_5, pathway_3)) %>%
  group_by(pathway_1, pathway_2) %>%
  summarise(total_times_fused_random =sum(times_fused_random),
            total_fused_percent_random = sum(times_fused_percent_random)) %>%
  ungroup %>%
  arrange(-total_times_fused_random); head(n_fusions_random_dedup)

#### merge frequencies and find enrichments
head(n_fusions_dedup)
head(n_fusions_random_dedup)
nrow(n_fusions_dedup); nrow(n_fusions_random_dedup)

n_with_random <- merge(n_fusions_random_dedup, n_fusions_dedup, by=c("pathway_1", "pathway_2"), all.x = TRUE) %>%
  mutate(total_times_fused = ifelse(is.na(total_times_fused), 0, total_times_fused),
         total_percent = ifelse(is.na(total_percent), 0, total_percent),
         enrichment = total_percent/total_fused_percent_random) %>%
  arrange(-enrichment) 

# add pathway labels and clean
enrichments <- merge(merge(n_with_random, pathway_names, by.x = "pathway_1", by.y = "Reactome_ID"), 
      pathway_names, by.x = "pathway_2", by.y = "Reactome_ID") %>%
  rename(description_1 = Description.x,
         description_2 = Description.y) %>%
  select(pathway_1, description_1, pathway_2, description_2, 
         total_times_fused, total_percent, total_times_fused_random, total_fused_percent_random, enrichment) %>%
  arrange(-enrichment)

enrichments_10_or_more_fusions <- enrichments %>% filter(total_times_fused >= 10)
write.csv(enrichments_10, "~/Projects/fusion_ML/analyses/pathway_joining/enrichment_of_pathway_joining_10_or_more_fusions.csv", quote = FALSE, row.names = FALSE)

# by enrichment
summary(enrichments$enrichment)
write.csv(enrichments, "~/Projects/fusion_ML/analyses/pathway_joining/enrichment_of_pathway_joining.csv", quote = FALSE, row.names = FALSE)

# examine depletions
enrichments_2_or_more_fusions <- enrichments %>% filter(total_times_fused >= 2)

# Metastatic fusions only
fusions <- read.table("~/Projects/fusion_ML/features/fusion_events_with_ensg_cancer_type_is_metastatic.txt", header=TRUE, sep="\t")
head(fusions); nrow(fusions)

metastatic <- fusions %>% filter(is_metastatic == 1)

### WHICH REACTOME PATHWAYS ARE FUSED?

fusions_with_pathways <- merge(merge(metastatic, pathways_with_annot, by.x = "ensg_5", by.y = "ensg", all.x = TRUE), 
                               pathways_with_annot, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>% 
  rename(pathway_3 = Reactome_ID.y, pathway_5 = Reactome_ID.x, description_5 = Description.x, description_3 = Description.y) %>% 
  distinct %>%
  select(fusion_id, ensg_5, pathway_5, description_5, ensg_3, pathway_3, description_3) %>% 
  filter(!is.na(pathway_5) & !is.na(pathway_3))

# how many fusions are left?
length(unique(fusions_with_pathways$fusion_id))

# largest possible number of fusions between pathways is number of fusions
max_possible_n_fusions <- length(unique(fusions_with_pathways$fusion_id)); max_possible_n_fusions

n_fusions <- fusions_with_pathways %>% group_by(pathway_5, description_5, pathway_3, description_3) %>% 
  summarise(times_fused = n_distinct(fusion_id),
            times_fused_percent = 100*times_fused/max_possible_n_fusions) %>% 
  ungroup %>% arrange(-times_fused)

# have to account for A,B order
n_fusions_dedup <- n_fusions %>% mutate_each(funs(as.character), pathway_5, pathway_3) %>%
  mutate(pathway_1 = pmin(pathway_5, pathway_3), 
         pathway_2 = pmax(pathway_5, pathway_3)) %>%
  group_by(pathway_1, pathway_2) %>%
  summarise(total_times_fused=sum(times_fused),
            total_percent = sum(times_fused_percent)) %>%
  ungroup %>%
  arrange(-total_times_fused); head(n_fusions_dedup)

n_with_random <- merge(n_fusions_random_dedup, n_fusions_dedup, by=c("pathway_1", "pathway_2"), all.x = TRUE) %>%
  mutate(total_times_fused = ifelse(is.na(total_times_fused), 0, total_times_fused),
         total_percent = ifelse(is.na(total_percent), 0, total_percent),
         enrichment = total_percent/total_fused_percent_random) %>%
  arrange(-enrichment) 

# add pathway labels and clean
enrichments <- merge(merge(n_with_random, pathway_names, by.x = "pathway_1", by.y = "Reactome_ID"), 
                     pathway_names, by.x = "pathway_2", by.y = "Reactome_ID") %>%
  rename(description_1 = Description.x,
         description_2 = Description.y) %>%
  select(pathway_1, description_1, pathway_2, description_2, 
         total_times_fused, total_percent, total_times_fused_random, total_fused_percent_random, enrichment) %>%
  arrange(-enrichment)

metastasis_3_or_more <- enrichments %>% filter(total_times_fused >= 3)

write.csv(enrichments, "~/Projects/fusion_ML/analyses/pathway_joining/metastasis_enrichment_of_pathway_joining.csv", quote = FALSE, row.names = FALSE)

# do metastatic parents have more pathways annotated?
primary_fusions <- fusions %>% filter(is_metastatic == 0)
primary_parents <- unique(c(as.character(primary_fusions$ensg_5), as.character(primary_fusions$ensg_3)))

metastatic_fusions <- fusions %>% filter(is_metastatic == 1)
metastatic_parents <- unique(c(as.character(metastatic_fusions$ensg_5), as.character(metastatic_fusions$ensg_3)))

# subset pathway counts values
head(pathway_counts)
metastatic_parents_n_pathways <- pathway_counts %>% filter(ensg %in% metastatic_parents)
primary_parents_n_pathways <- pathway_counts %>% filter(!(ensg %in% metastatic_parents) & ensg %in% primary_parents)

# summary
summary(primary_parents_n_pathways$n_pathways)
summary(metastatic_parents_n_pathways$n_pathways)
wilcox.test(metastatic_parents_n_pathways$n_pathways, primary_parents_n_pathways$n_pathways)


####### Get pathway joinings for fusion events in different cancer types

# get pathway info
fusions_with_pathways <- merge(merge(fusions, pathways_with_annot, by.x = "ensg_5", by.y = "ensg", all.x = TRUE), 
                               pathways_with_annot, by.x = "ensg_3", by.y = "ensg", all.x = TRUE) %>% 
  rename(pathway_3 = Reactome_ID.y, pathway_5 = Reactome_ID.x, description_5 = Description.x, description_3 = Description.y) %>% 
  distinct %>%
  select(fusion_id, ensg_5, pathway_5, description_5, ensg_3, pathway_3, description_3) %>% 
  filter(!is.na(pathway_5) & !is.na(pathway_3))

# get cancer type information
fusions_with_type <- merge(fusions_with_pathways, fusions %>% select(fusion_id, cancer_type), by = "fusion_id")
head(fusions_with_type)
head(fusions_with_pathways)

# largest possible number of fusions between pathways is number of fusions
n_by_cancer <- fusions_with_type %>% group_by(cancer_type) %>% summarise(total_fusion_count = n_distinct(fusion_id))
fusions_with_type %>% select(fusion_id) %>% distinct

n_fusions <- fusions_with_type %>% group_by(cancer_type, pathway_5, description_5, pathway_3, description_3) %>% 
  summarise(times_fused = n_distinct(fusion_id)) %>% 
  ungroup %>% arrange(-times_fused)

n_fusions <- merge(n_fusions, n_by_cancer, by="cancer_type")
n_fusions <- n_fusions %>%  mutate(percent = 100*times_fused/total_fusion_count)

# have to account for A,B order
n_fusions_dedup <- n_fusions %>% mutate_each(funs(as.character), pathway_5, pathway_3) %>%
  mutate(pathway_1 = pmin(pathway_5, pathway_3), 
         pathway_2 = pmax(pathway_5, pathway_3)) %>%
  group_by(cancer_type, pathway_1, pathway_2) %>%
  summarise(total_times_fused = sum(times_fused),
            total_percent = sum(percent)) %>%
  ungroup %>%
  arrange(-total_times_fused); head(n_fusions_dedup)

# with random
n_with_random <- merge(n_fusions_random_dedup, n_fusions_dedup, by=c("pathway_1", "pathway_2"), all.x = TRUE) %>%
  mutate(total_times_fused = ifelse(is.na(total_times_fused), 0, total_times_fused),
         total_percent = ifelse(is.na(total_percent), 0, total_percent),
         enrichment = total_percent/total_fused_percent_random) %>%
  arrange(-enrichment); head(n_with_random)

# add pathway labels and clean
enrichments <- merge(merge(n_with_random, pathway_names, by.x = "pathway_1", by.y = "Reactome_ID"), 
                     pathway_names, by.x = "pathway_2", by.y = "Reactome_ID") %>%
  rename(description_1 = Description.x,
         description_2 = Description.y) %>%
  select(cancer_type, pathway_1, description_1, pathway_2, description_2, 
         total_times_fused, total_percent, total_times_fused_random, total_fused_percent_random, enrichment) %>%
  arrange(cancer_type, -enrichment)

head(enrichments)

cancer_type_5_or_more <- enrichments %>% filter(total_times_fused >= 5)


# # create mostly sparse matrix of pathway participation
# pathway_participation_matrix <- pathways_genes %>% 
#   mutate(participation = 1) %>% 
#   spread(key = Reactome_ID, value = participation, fill = 0)
# 
# # print out this sparse matrix
# write.table(pathway_participation_matrix, "~/Projects/fusion_ML/features/REACTOME_pathway_participation_matrix.csv", quote = FALSE, row.names = FALSE)
# 
# dim(pathway_participation_matrix)

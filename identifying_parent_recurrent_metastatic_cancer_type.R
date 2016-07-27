# Fusion event data
require(dplyr)
require(tidyr)
require(ggplot2)

############ Summarize parent and fusion event data
fusions <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/klijn_gene_fusions.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE); head(fusions)
head(fusions)

# convert fuson entrez to ensg
ensg_entrez <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
head(ensg_entrez)

# label conversion
fusions_ensg <- merge(merge(fusions, ensg_entrez, by.x = "entrez_id_5", by.y = "entrez", all.x = TRUE), 
                      ensg_entrez, by.x = "entrez_id_3", by.y = "entrez", all.x = TRUE) %>% rename(ensg_3 = ensg.y, ensg_5 = ensg.x) %>% distinct()   
head(fusions_ensg)

# Anything we can do to assign ENSGs to the remaining entrez?
missing_5 <- filter(fusions_ensg, is.na(ensg_5)) %>% select(entrez_id_5)
missing_3 <- filter(fusions_ensg, is.na(ensg_3)) %>% select(entrez_id_3)
missing <- sort(unique(c(missing_5$entrez_id_5, missing_3$entrez_id_3))); missing
write.csv(missing, "~/Projects/fusion_ML/analyses/missing_fusion_parent_entrez_ensg.txt", quote = FALSE, row.names = FALSE)

# sort out 5' and 3' business
five <- data.frame(fusions_ensg$ensg_5, 
                   rep(1, length(fusions_ensg$ensg_5)),
                   rep(0, length(fusions_ensg$ensg_3))); names(five) <- c("ensg", "is_5_prime_parent", "is_3_prime_parent")

three <- data.frame(fusions_ensg$ensg_3, 
                   rep(1, length(fusions_ensg$ensg_3)),
                   rep(0, length(fusions_ensg$ensg_5))); names(three) <- c("ensg", "is_3_prime_parent", "is_5_prime_parent")

parent_genes <- rbind(five, three) %>% na.omit %>% distinct %>% arrange(ensg); head(parent_genes); nrow(parent_genes)

# write out clean set of parent_genes (remove those that act as both 3' and 5' parents)
clean_parent_genes <- parent_genes %>% select(ensg) %>% mutate(is_parent = 1) %>% distinct %>% na.omit

write.table(clean_parent_genes,
            "~/Projects/fusion_ML/features/parent_genes.csv", sep=",", quote = FALSE, row.names = FALSE)


############# write out 5' and 3' gene list ############# 
both_5_and_3_parents <- parent_genes %>% group_by(ensg) %>% summarise(count = n()) %>% filter(count==2)
parent_genes_5_3_prime <- rbind(parent_genes %>% filter(!ensg %in% both_5_and_3_parents$ensg),
                                parent_genes %>% filter(ensg %in% both_5_and_3_parents$ensg) %>% mutate(is_5_prime_parent=1, is_3_prime_parent=1) %>% arrange(ensg) %>% distinct) 

parent_genes_5_3_prime_clean <- parent_genes_5_3_prime %>% na.omit %>% arrange(ensg)
head(parent_genes_5_3_prime_clean); nrow(parent_genes_5_3_prime_clean)

write.table(parent_genes_5_3_prime_clean, "~/Projects/fusion_ML/features/five_three_prime_parents.csv", sep=",", quote = FALSE, row.names = FALSE)

############# recurrent parent genes ############# 
parent_genes <- rbind(five, three) %>% na.omit
recurrent_parents <- parent_genes %>% select(ensg) %>% na.omit %>% group_by(ensg) %>% 
  summarise(num_fusions = n()) %>% arrange(-num_fusions) 

# count n fusions
recurrent_parents <- recurrent_parents %>% 
  mutate(is_recurrent_parent = ifelse(num_fusions>=2, 1, 0)) %>% 
  distinct  %>% select(ensg, is_recurrent_parent, num_fusions)

write.table(recurrent_parents,
            "~/Projects/fusion_ML/features/recurrent_parent_genes.csv", sep=",", quote = FALSE, row.names = FALSE)


#############       FUSION EVENTS       ##############
#############   metastasis + cancer type #############  
# clean to ensure we only keep fusion events where we know both ENSGs
fusions_ensg <- fusions_ensg %>% filter(!is.na(ensg_5) & !is.na(ensg_3))

# merge cell line data
cell_line_data <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/klijn_cell_lines_simple.csv", header=TRUE, sep=",", na.strings = "<NA>", stringsAsFactors = FALSE); head(fusions)
head(cell_line_data)

# join on cell line information
fusions_with_cell_line <- merge(fusions_ensg, cell_line_data, by="cell_line", all.x = TRUE) %>%
  rename(fusion_id = id.x)  %>% 
  select(fusion_id, gene_5, gene_3, ensg_5, ensg_3, cell_line, tissue, metastatic_tissue, primary_tissue) %>%
  mutate(is_metastatic = ifelse(metastatic_tissue=="NA" | is.na(metastatic_tissue) | metastatic_tissue=="NULL", 0, 1)) %>% 
  rename(cancer_type = tissue)
head(fusions_with_cell_line)

# count fusion IDs, reflects multiple Ensg
fusions_with_cell_line %>% select(fusion_id) %>% nrow
fusions_with_cell_line %>% select(fusion_id) %>% distinct %>% nrow

# check metastatic labels are ok
nrow(fusions_with_cell_line %>% filter(is_metastatic == 1) %>% select(fusion_id) %>% distinct)
nrow(fusions_with_cell_line %>% filter(is_metastatic == 0) %>% select(fusion_id) %>% distinct)

# counts, percentages
nrow(fusions_with_cell_line %>% filter(is_metastatic == 1) %>% select(fusion_id) %>% distinct)/length(unique(fusions_with_cell_line$fusion_id))
nrow(fusions_with_cell_line %>% filter(is_metastatic == 0) %>% select(fusion_id) %>% distinct)/length(unique(fusions_with_cell_line$fusion_id))

# how many fusion events don't have cell line info?
fusions_with_cell_line %>% filter(is.na(cell_line))

# write out fusion event data
head(fusions_with_cell_line)
write.table(fusions_with_cell_line %>% arrange(fusion_id),
            "~/Projects/fusion_ML/features/fusion_events_with_ensg_cancer_type_is_metastatic.txt", sep="\t", quote = FALSE, row.names = FALSE)

### Plot number of cell lines per tissue supergroup
head(fusions_with_cell_line)
# group and count
n_fusions_by_cancer <- fusions_with_cell_line %>% 
  select(fusion_id, cancer_type) %>% distinct %>%
  group_by(cancer_type) %>% 
  summarise(fusion_count = n()) %>% 
  arrange(-fusion_count) %>%
  mutate(y_position = 1) %>%
  mutate(new_tissue_label = paste(fusion_count, "_", cancer_type, sep=""))

head(n_fusions_by_cancer)

# do bar plot
theme_set(theme_bw(base_size = 20))

p <- ggplot(n_fusions_by_cancer, 
            aes(x=reorder(cancer_type, -fusion_count), y=fusion_count, fill=reorder(new_tissue_label, -fusion_count))) + 
  geom_bar(width = 0.5, stat = "identity", show.legend = FALSE) + xlab("cancer type") + ylab("number of fusion events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p

ggsave(plot = p, filename = "n_fusions_by_cancer_type.pdf", path = "~/Projects/fusion_ML/figures/", width = 6, height = 4)

# 
# ############### Summarize tissue migrations with circos plot 
# # http://mkweb.bcgsc.ca/tableviewer/
# klijn %>% group_by(primary_tissue) %>% summarise(n = n())
# klijn %>% group_by(metastatic_tissue) %>% summarise(n = n())
# 
# combo <- klijn %>% 
#   group_by(primary_tissue, metastatic_tissue) %>%
#   summarise(n = n()) %>%
#   na.omit() %>% filter(metastatic_tissue != "NULL"); combo
# 
# circos_data <- combo %>% spread(metastatic_tissue, n, fill=0)
# write.table(circos_data, "~/Projects/fusion_ML/figures/circos_data.txt", sep="\t", quote=FALSE, row.names=FALSE)

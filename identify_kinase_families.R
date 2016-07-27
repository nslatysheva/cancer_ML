# Identify kinase genes and kinase families

require(dplyr)
require(stringr)

# Get kinase data
kinase <- read.table("~/Projects/fusion_ML/data/gene_identity_data/kinases/kinases.csv", header=TRUE, sep=',')
head(kinase); nrow(kinase)

# Trim whitespace surrounding uniprot accessions
kinase$uniprot <- str_trim(kinase$uniprot, side = c("both", "left", "right"))

# Convert to ensg
uniprot_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_uniprot_conversion_global.txt", header=TRUE, sep="\t"); head(uniprot_ensg)
kin_genes <- merge(kinase, uniprot_ensg, on="uniprot", all.x = TRUE); head(kin_genes,500); nrow(kin_genes)

# Anything we can do to assign ENSGs to the remaining uniprots?
filter(kin_genes, is.na(ensg)) %>% select(uniprot)

# What rows are repeated?
# I think this is fine as is
kin_genes %>%
  group_by(uniprot) %>%
  summarise(
    count = n()) %>%
  arrange(-count)

filter(kin_genes, uniprot=="Q08345")

# write out kinase gene list
kin_genes$is_kinase <- 1
write.table(na.omit(select(kin_genes, ensg, is_kinase)), 
            "~/Projects/fusion_ML/features/kinase_genes.csv", sep=",", quote = FALSE, row.names = FALSE)

# create dummy variables denoting presence of each kinase class
kin_classes <- na.omit(select(kin_genes, ensg, family)); head(kin_classes); dim(kin_classes)
expanded_kin_classes <- cbind(kin_classes, model.matrix( ~ family - 1, data=kin_classes))

# clean up column names
expanded_kin_classes_clean <- select(expanded_kin_classes, -family) 
names(expanded_kin_classes_clean)[2:length(names(expanded_kin_classes_clean))] = names(expanded_kin_classes_clean)[2:length(names(expanded_kin_classes_clean))] %>% 
  gsub(pattern = "family", replacement = "kinase_")
head(expanded_kin_classes_clean)

write.table(expanded_kin_classes_clean, "~/Projects/fusion_ML/features/kinase_classes.csv", sep=",", quote = FALSE, row.names = FALSE)


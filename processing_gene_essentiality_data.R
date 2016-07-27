# Processing gene essentiality data

require(dplyr)
require(tidyr)

essential <- read.table("~/Projects/fusion_ML/data/gene_identity_data/essentiality/gene_essentiality.txt", sep="\t", header = TRUE)
essential_human <- essential %>% filter(sciName == "Homo sapiens")

essential_clean <- 
  essential_human %>% filter(essential == "Y") %>%
  mutate(is_essential = 1) %>% 
  select(locus, is_essential) %>% na.omit() %>% distinct() %>% 
  rename(ensg=locus)

head(essential_clean)

# write out feature set
write.csv(essential_clean, "~/Projects/fusion_ML/features/essential_genes.csv", quote = FALSE, row.names = FALSE)



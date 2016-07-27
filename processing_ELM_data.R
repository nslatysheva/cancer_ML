# Processing linear motif data
require(dplyr)
require(tidyr)

# Get data
lm <- read.table("~/Projects/fusion_ML/data/structural_data/linear_motifs/human_elm_instances.csv", header=TRUE, sep=',')
head(lm); nrow(lm); n_distinct(lm$Accession)

# generate unique structural accession
lm <-
  lm %>% mutate(LM_length = (End-Start)+1,
                unique_lm_accession = paste(Primary_Acc, Start, End, sep = "_")) %>%
  select(unique_lm_accession, ELMType, Primary_Acc, LM_length)

# Fetch number of lms per proteins and number of aas covered
lm_count_coverage <-
  lm %>%
  group_by(Primary_Acc) %>% 
  summarise(
    num_lms = n_distinct(unique_lm_accession)
 #   aa_covered_lms = sum(LM_length)
  ) %>%
  rename(uniprot=Primary_Acc) %>% arrange(-num_lms)

head(lm_count_coverage); nrow(lm_count_coverage)

# convert to ensg
uniprot_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_uniprot_conversion_global.txt", header=TRUE, sep="\t"); head(uniprot_ensg)
lm_count_coverage_ensg <- merge(lm_count_coverage, uniprot_ensg, on="uniprot", all.x = TRUE); head(lm_count_coverage_ensg)

# try to fix stranded uniprots
unmapped_uniprot_ELM_mapping <- lm_count_coverage_ensg %>% filter(is.na(ensg)) %>% select(uniprot)
unmapped_uniprot_ELM_mapping %>% filter(grepl("-", uniprot))

write.table(unmapped_uniprot_ELM_mapping, "~/Projects/fusion_ML/data/structural_data/linear_motifs/unmapped_uniprot_ELM_mapping.csv", 
            sep=",", quote = FALSE, col.names = FALSE, row.names = FALSE)

# output LM results
lm_features_clean <- lm_count_coverage_ensg %>% na.omit() %>% select(ensg, num_lms)
write.table(lm_features_clean, "~/Projects/fusion_ML/features/lm_features_ELM.csv", 
            sep=",", quote = FALSE, row.names = FALSE)


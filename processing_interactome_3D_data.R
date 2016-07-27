# Getting proportion of protein involved in interactions
require(dplyr)
require(tidyr)
require(igraph)

si <- read.csv("~/Projects/fusion_ML/data/structural_data/structural_interacting_segments_interactome3D/structural_interactions_interactome3D.txt", sep="\t")
head(si); nrow(si)

# get just the uniprots and sequence cutoffs
si <- si %>% select(PROT1, PROT2, SEQ_BEGIN1, SEQ_END1, SEQ_BEGIN2, SEQ_END2) %>% distinct()
head(si); nrow(si)

# we want to do this gene by gene, reformat df
prot1 <- si %>% select(PROT1, SEQ_BEGIN1, SEQ_END1) %>% rename(uniprot = PROT1, begin = SEQ_BEGIN1, end = SEQ_END1)
prot2 <- si %>% select(PROT2, SEQ_BEGIN2, SEQ_END2) %>% rename(uniprot = PROT2, begin = SEQ_BEGIN2, end = SEQ_END2)
prot <- rbind(prot1, prot2) %>% distinct(); head(prot)

# there are indeed plenty of overlapping interacting regions
# for starters, I will just count the longest stretch per protein
prot_summary <- 
  prot %>% mutate(stretch = end - begin) %>%
  group_by(uniprot) %>%
  summarise(
    longest_stretch = max(stretch),
    number_of_stretches = n()
  )

head(prot_summary)

# Need to convert this to Ensembl ids
uniprot_ensembl <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_uniprot_conversion_global.txt", sep="\t", header=TRUE)
prot_ensg <- merge(prot_summary, uniprot_ensembl, on="uniprot", all.x=TRUE)

# try to reduce number of uniprots not mapped. Reduced 202 down to 145.
unmapped_uniprot <- prot_ensg %>% filter(is.na(ensg)) %>% select(uniprot) %>% distinct(); unmapped_uniprot
write.table(unmapped_uniprot, "~/Projects/fusion_ML/data/structural_data/unmapped_uniprot.csv", sep=",", row.names = FALSE, quote = FALSE, col.names = FALSE)

# get just ensg and write out
interacting_tracts <- prot_ensg %>% select(ensg, longest_stretch, number_of_stretches) %>% na.omit() %>% distinct

# find where duplicates occur
duplicates <- interacting_tracts %>% group_by(ensg) %>% summarise(count=n()) %>% arrange(-count) %>% filter(count > 1)
duplicates %>% select(ensg) %>% distinct %>% nrow
interacting_tracts %>% filter(ensg %in% duplicates$ensg) %>% arrange(ensg)

# deduplicate, take max values
interacting_tracts_clean <- interacting_tracts %>% group_by(ensg) %>% summarise(longest_stretch = max(longest_stretch),
                                                                                number_of_stretches = max(number_of_stretches))

interacting_tracts_clean %>% group_by(ensg) %>% summarise(count=n()) %>% arrange(-count) %>% filter(count > 1)

# examine what happened to dup values
interacting_tracts_clean %>% filter(ensg %in% duplicates$ensg) %>% arrange(ensg)


write.table(interacting_tracts_clean, "~/Projects/fusion_ML/features/interacting_tracts.csv", sep=",", quote = FALSE, row.names = FALSE)

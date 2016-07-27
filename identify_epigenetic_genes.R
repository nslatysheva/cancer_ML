# Identify TF genes and TF families
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)

#### Get epigenetic data
epi <- read.csv("~/Projects/fusion_ML/data/gene_identity_data/epigenetic/EpiGenes_main_clean.csv", header=TRUE, stringsAsFactors = FALSE); head(epi)
epi <- epi %>% select(Id, UniProt_ID, Function) %>% filter(Function != "#")
head(epi)

# regroup epi data into fewer categories 
epi %>% group_by(Function) %>% summarise(count=n()) %>% arrange(-count) %>% nrow

epi_new <- epi %>% 
  mutate(is_chromatin_remodelling = ifelse(grepl(x = Function, pattern="Chromatin remodelling"), 1, 0)) %>%
  mutate(is_histone_modification = ifelse(grepl(x = Function, pattern="Histone modification"), 1, 0)) %>% 
  mutate(is_other_epigenetic = ifelse(is_chromatin_remodelling==0 & is_histone_modification==0, 1, 0)) %>%
  mutate_each(funs(as.factor))

#### epigenetic complexes
epi_complexes <- read.csv("~/Projects/fusion_ML/data/gene_identity_data/epigenetic/EpiGenes_complexes_clean.csv", header=TRUE); head(epi_complexes)

# break up by comma OR the bar character (have to escape out the bar since it's normally logical OR)
epi_complexes_split <- epi_complexes %>% select (Id, UniProt_ID) %>% 
  separate(into = c(seq(1:37)), col = UniProt_ID, sep = ",|\\|") %>% 
  gather(key=complex_member, value=is_present, -Id, na.rm = TRUE) %>% 
  arrange(Id) %>% rename(UniProt_ID = is_present)

# accessions are all still messed up
# get rid of spaces, question marks, parentheses and plusses
epi_complexes_split$UniProt_ID <- gsub(pattern = "\\s*|\\?|\\(|\\)|\\+", replacement = "", epi_complexes_split$UniProt_ID)

# these IDs are a bit non-strandard, ensembl biomart won't convert
# going to use http://www.uniprot.org/uploadlists/
# which IDs need mapping

unmapped_IDs <- rbind(epi_new %>% select(UniProt_ID),
                      epi_complexes_split %>% select(UniProt_ID)) %>% distinct

write.table(unmapped_IDs, "~/Projects/fusion_ML/data/gene_identity_data/epigenetic/uniprot_IDs_to_map.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)

# got the conversions finally
ensg_uniprot_species <- read.table("~/Projects/fusion_ML/data/gene_identity_data/epigenetic/uniprot_species_ensg_mapping", sep="\t", header=TRUE, stringsAsFactors = FALSE); head(ensg_uniprot_species)

# get the epigenetic gene table out
epi_ensg <- inner_join(ensg_uniprot_species, epi_new, by="UniProt_ID"); head(epi_ensg)
write.table(epi_ensg %>% select(ensg, chromatin_remodelling:other_epigenetic) %>% distinct, 
            "~/Projects/fusion_ML/features/epigenetic_gene.csv", sep=",", quote = FALSE, row.names = FALSE)

# get the epigentic complex table out
epi_complexes_ensg <- inner_join(ensg_uniprot_species, epi_complexes_split, by="UniProt_ID")
write.table(epi_complexes_ensg %>% select(ensg) %>% mutate(is_epigenetic_complex = 1) %>% distinct, 
            "~/Projects/fusion_ML/features/epigenetic_complexes.csv", sep=",", quote = FALSE, row.names = FALSE)

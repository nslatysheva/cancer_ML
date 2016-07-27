# subcellular localization mapping
require(dplyr)
require(stringr)
require(tidyr)

#### Get experimental and knowledge data
knowledge <- read.table("~/Projects/fusion_ML/data/gene_identity_data/localization/human_compartment_knowledge_full.tsv", header=TRUE, stringsAsFactors = FALSE, sep="\t"); head(knowledge)
exp <- read.table("~/Projects/fusion_ML/data/gene_identity_data/localization/human_compartment_experiments_full.tsv", header=TRUE, stringsAsFactors = FALSE, sep="\t"); head(exp)

# combine
local <- rbind(knowledge %>% select(ensp, gene, localization),
               exp %>% select(ensp, gene, localization)) %>% distinct

# how many data points, for how many proteins?
nrow(local)
local %>% select(ensp) %>% distinct %>% nrow

# how many groups of localizations?
groups <- local %>% group_by(localization) %>% summarise(count = n()) %>% arrange(-count)

# categories to keep 
categories <- c("Cytoplasm", "Nucleus", "Plasma membrane", "Organelle lumen", "Vesicle", 
                "Endomembrane system", "Extracellular region", "Cytoskeleton", "Nucleolus", "Cell junction")

local_10 <- local %>% filter(localization %in% categories)
local_10 %>% group_by(localization) %>% summarise(count = n()) %>% arrange(-count)

# great! Now convert to ensemble gene IDs
ensp_ensg <- read.csv("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_ensp_conversion_protein_coding.csv", header=TRUE)

local_10_ensg <- inner_join(local_10, ensp_ensg, by="ensp")

write.table(local_10_ensg %>% filter(is.na(ensg)) %>% select(ensp) %>% distinct,
            "~/Projects/fusion_ML/data/gene_identity_data/localization/unmapped_ensp.csv",
            row.names = FALSE, quote = FALSE, sep=",")

local_10_ensg %>% group_by(ensg) %>% summarise(p = n_distinct(ensp)) %>% arrange(-p)

# Build matrix of things that are in this category
ensg_local <- local_10_ensg %>% select(ensg, localization) %>% distinct
ensg_localization_matrix <- ensg_local %>% mutate(is_in_location = 1) %>% spread(key = localization, value = is_in_location, fill = 0)

# Make sure column names are beautiful
clean_colnames <- names(ensg_localization_matrix) %>% tolower %>% gsub(pattern=" ", replacement="_")
names(ensg_localization_matrix)[2:ncol(ensg_localization_matrix)] <- paste("loc_", clean_colnames[2:length(clean_colnames)], sep="")

# write out features
write.csv(ensg_localization_matrix, "~/Projects/fusion_ML/features/subcellular_localization.csv", quote = FALSE, row.names = FALSE)

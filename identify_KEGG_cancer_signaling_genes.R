# clean up cancer signalling data
require(dplyr)
require(stringr)
require(tidyr)

# read in raw results for KEGG pathway
kegg_cancer <- read.table("~/Projects/fusion_ML/data/pathway_data/raw_KEGG_cancer_pathway_gene_data.txt", sep="\n", header=TRUE)

# clean data into a df
odd_indices <- seq(from = 1, to = nrow(kegg_cancer)-1, by = 2) 
even_indices <- seq(from = 2, to = nrow(kegg_cancer), by = 2) 
kegg_cancer_df <- data.frame(kegg_cancer[odd_indices,], kegg_cancer[even_indices,])
names(kegg_cancer_df) <- c("entrez", "gene_name")

# column 2 isn't actually needed, just convert to ensg now
# inner_join == left_join here, no missing ENSGs
kegg_cancer_df$entrez <- kegg_cancer_df$entrez %>% gsub(pattern=" ", replacement="")
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", sep="\t", header=TRUE)
kegg_cancer_df <- inner_join(kegg_cancer_df, entrez_ensg, by="entrez"); head(kegg_cancer_df)

# clean up and write out
clean_kegg_cancer <- kegg_cancer_df %>% mutate(is_KEGG_cancer_pathway=1) %>% select(ensg, is_KEGG_cancer_pathway) %>% distinct
write.csv(clean_kegg_cancer, "~/Projects/fusion_ML/features/KEGG_cancer_pathway_genes.csv", quote = FALSE, row.names = FALSE)


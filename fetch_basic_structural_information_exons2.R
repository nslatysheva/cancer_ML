#Fetch basic gene structural information
#Average gene length, number of transcripts, number of proteins (should be the same)

require(RMySQL)
require(dplyr)
require(ggplot2)

source("~/Projects/fusion_ML/scripts/database_credentials.R")

gene_lengths <- query("
                        SELECT 
                        ensembl_gene_id,
                        AVG(length) AS avg_length,
                        COUNT(DISTINCT ensembl_protein_id) AS count_ENSPs
                        FROM
                        natashal.ensp_ensg_gene_length
                        GROUP BY ensembl_gene_id
                        ")

gene_lengths <- read.table("~/Projects/fusion_ML/data/structural_data/avg_gene_lengths_num_ENSPs.txt", header=TRUE, sep='\t')
head(gene_lengths); nrow(gene_lengths)

# Get number of transcripts + avg number of exons per gene
isoform_counts <- read.table("~/Projects/fusion_ML/data/structural_data/exons_protein_coding.txt", header=TRUE, sep='\t')
head(isoform_counts); nrow(isoform_counts)

# How many genes are actually here?
n_distinct(isoform_counts$Ensembl.Gene.ID)

num_ensp_avg_exons <-
  isoform_counts %>%
  filter(Ensembl.Protein.ID != "") %>%
  group_by(Ensembl.Gene.ID, Ensembl.Protein.ID) %>%
  summarise(
    num_exons = n_distinct(Ensembl.Exon.ID)) %>%
  group_by(Ensembl.Gene.ID) %>%
  summarise(
    num_proteins = n_distinct(Ensembl.Protein.ID),
    avg_num_exons = mean(num_exons)
  ) %>%
  rename(
    ensg=Ensembl.Gene.ID
  )

nrow(num_ensp_avg_exons)
head(gene_lengths)

# combine gene lengths, count ensps, count transcripts, avg exons
basic_structural_info <- merge(gene_lengths, num_ensp_avg_exons, by="ensg", all.x=TRUE)
head(basic_structural_info, 1000);nrow(basic_structural_info)

nrow(na.omit(basic_structural_info))

basic_structural_info[is.na(basic_structural_info),]
write.table(basic_structural_info, "~/Projects/fusion_ML/features/basic_structural_info.csv", sep=",", quote = FALSE, row.names = FALSE)


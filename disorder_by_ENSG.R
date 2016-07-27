# Molecular feature figure - intrinsic disorder in parents

require(RMySQL)
require(plyr)
require(reshape2)
require(ggplot2)

source("~/Projects/fusion_ML/scripts/database_credentials.R")

# configure queries
query <- function(query) { a = dbGetQuery(con, statement=query) }

# get number by ensp, group by ensg
# note that zero densities are necessarily excluded here for speed
disorder_by_gene <- query("
                            SELECT 
                                ensembl_gene_id,
                                MAX(length) AS max_length,
                                AVG(average_disorder) AS disorder
                            FROM
                                (SELECT 
                                    ensembl_gene_id, protein, length, average_disorder
                                FROM
                                    natashal.disorder_scores_human_proteome_avg_by_protein d
                                JOIN ensp_ensg_gene_length l ON d.protein = l.ensembl_protein_id) x
                            GROUP BY ensembl_gene_id
                             ")

# April 2016 Update
# I think the disorder calc suffers from the max length problem
# correct version:
disorder_by_gene <- query("
                          SELECT 
                          all_lengths.*, max_lengths.max_length
                          FROM
                          (SELECT 
                          ensembl_gene_id,
                          protein,
                          length,
                          AVG(average_disorder) AS disorder
                          FROM
                          (SELECT 
                          ensembl_gene_id, protein, length, average_disorder
                          FROM
                          natashal.disorder_scores_human_proteome_avg_by_protein d
                          JOIN ensp_ensg_gene_length l ON d.protein = l.ensembl_protein_id) x
                          GROUP BY ensembl_gene_id , protein) all_lengths
                          INNER JOIN
                          (SELECT 
                          ensembl_gene_id, MAX(length) AS max_length
                          FROM
                          natashal.ensp_ensg_gene_length
                          GROUP BY ensembl_gene_id) max_lengths ON all_lengths.ensembl_gene_id = max_lengths.ensembl_gene_id
                          AND all_lengths.length = max_lengths.max_length
                          GROUP BY ensembl_gene_id
                          ")

# offline
disorder_by_gene <- read.table("~/Projects/fusion_ML/data/structural_data/disorder/avg_disorder_by_ensg.txt", header=TRUE, sep='\t')

# write out clean feature
disorder_by_gene_clean <- 
  disorder_by_gene %>% select(ensembl_gene_id, disorder) %>% 
  rename(ensg = ensembl_gene_id, avg_disorder = disorder)

write.table(disorder_by_gene_clean, "~/Projects/fusion_ML/features/avg_disorder_by_ensg.csv", sep=",", quote = FALSE, row.names = FALSE)



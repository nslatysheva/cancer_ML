# Pfam domains in genes

require(dplyr)

# Get pfam data
pfam <- read.table("~/Projects/fusion_ML/data/structural_data/pfam/Pfam_domains_with_ENSP.txt", header=TRUE, sep='\t')
head(pfam); nrow(pfam); n_distinct(pfam$Ensembl.Gene.ID)

# Remove NA rows
pfam <- na.omit(pfam)

# Fetch number of pfam domains per protein and coverage of protein covered by aa
num_pfam_domains_by_ensg_ensp <-
  pfam %>%
  mutate(
    domain_accession = paste(Ensembl.Protein.ID, "_", Pfam.ID, "_", Pfam.start, "_", Pfam.end, sep = ''),
    domain_length = Pfam.end-Pfam.start
  ) %>%
  group_by(Ensembl.Gene.ID, Ensembl.Protein.ID) %>% 
  summarise(
    num_domains = n_distinct(domain_accession),
    aa_covered_domains = sum(domain_length),
    avg_domain_length = mean(domain_length)
  ) %>%
  rename(ensg = Ensembl.Gene.ID, ensp=Ensembl.Protein.ID)
  
# group by ensg
num_pfam_domains_by_ensp <-
  num_pfam_domains_by_ensg_ensp %>%
  group_by(ensg) %>%
  summarise(
    avg_num_domains = mean(num_domains),
    avg_domain_length = mean(avg_domain_length)
  )

head(num_pfam_domains_by_ensp)

# write out pfam features
write.table(num_pfam_domains_by_ensp, "~/Projects/fusion_ML/features/pfam.csv", sep=",", quote = FALSE, row.names = FALSE)


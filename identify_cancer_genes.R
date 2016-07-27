# Identify cancer genes
require(dplyr)
require(stringr)
require(tidyr)

# Get data
cancer <- read.table("~/Projects/fusion_ML/data/gene_identity_data/cancer_genes/network_cancer_genes.csv", header=TRUE, sep=',')
nrow(cancer); head(cancer)

# Convert to ensg
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)
cancer_genes <- merge(cancer, entrez_ensg, on="entrez", all.x = TRUE); head(cancer_genes, 500); nrow(cancer_genes)

# Anything we can do to assign ENSGs to the remaining genes? Down from 47 to 26 NAs
unmapped_genes <- filter(cancer_genes, is.na(ensg)) %>% select(entrez); unmapped_genes
write.table(unmapped_genes, "~/Projects/fusion_ML/data/gene_identity_data/cancer_genes/unmapped_genes.txt",  sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Any original entries now missing? All ok.
n_distinct(cancer$symbol); n_distinct(cancer_genes$symbol)
cancer$symbol[!cancer$symbol %in% unique(cancer_genes$symbol)]

# Any original entries duplicated? This is fine.
cancer_genes %>% group_by(symbol) %>%
  summarise(
    n_ensg = n_distinct(ensg)
  ) %>% arrange(-n_ensg)

filter(cancer_genes, symbol == "HLA-A")

# reformat and write out gene list
cancer_list <- 
  cancer_genes %>%
  mutate(
    is_cancer = 1
  ) %>%
  na.omit() %>%
  distinct()

head(cancer_list,10); nrow(cancer_list)

write.table(select(cancer_list, ensg, is_cancer) %>% distinct(), 
            "~/Projects/fusion_ML/features/cancer_genes.csv", sep=",", quote = FALSE, row.names = FALSE)

#### More stringent list of cancer genes from cancer gene census (COSMIC)
# Get data
cancer_stringent <- read.table("~/Projects/fusion_ML/data/gene_identity_data/cancer_genes/cancer_gene_census.csv", header=TRUE, sep=',')
nrow(cancer_stringent); head(cancer_stringent); cancer_stringent$entrez <- as.factor(cancer_stringent$entrez)

# simplify data frame
cancer_stringent <- cancer_stringent %>% select(entrez, somatic, germline, tissue_type)

# Convert to ensg
entrez_ensg <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, sep="\t"); head(entrez_ensg)
cancer_stringent <- left_join(cancer_stringent, entrez_ensg, by="entrez"); head(cancer_stringent, 500); nrow(cancer_stringent)

# Anything we can do to assign ENSGs to the remaining genes? Best I can do already
unmapped_genes <- filter(cancer_stringent, is.na(ensg)) %>% select(entrez); unmapped_genes
write.table(unmapped_genes, "~/Projects/fusion_ML/data/gene_identity_data/cancer_genes/unmapped_genes_cgc.txt",  sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# deal with tissue type field
cgc_summary <- cancer_stringent %>% na.omit() %>% select(ensg, somatic, germline, tissue_type)

# get rid of trailing spaces, separete tissue type column into separate columns for spreading
cancer_types <- 
  cgc_summary %>% 
  mutate(tissue_type = gsub(" ", "", tissue_type)) %>%
  separate(tissue_type, into = c("tissue1", "tissue2", "tissue3", "tissue4"), sep=",",) %>%
  filter(tissue1 != "") %>%
  distinct()

cancer_types[is.na(cancer_types)] <- 0; head(cancer_types,20)

# spread cancer types and process resulting mess
spread_cancer_types <- 
  cancer_types %>% 
  mutate(present1=1) %>% spread(tissue1, present1, fill=0) %>%
  rename(E1=E, L1=L, M1=M, O1=O) %>%
  mutate(present2 = ifelse(tissue2 !=0, 1, 0)) %>% spread(tissue2, present2, fill=0) %>%  
  rename(E2=E, L2=L, M2=M, O2=O) %>% select(-10) %>%
  mutate(present3 = ifelse(tissue3 !=0, 1, 0)) %>% spread(tissue3, present3, fill=0) %>%
  rename(E3=E, L3=L, M3=M, O3=O) %>% select(-13) %>%
  mutate(present4 = ifelse(tissue4 !=0, 1, 0)) %>% spread(tissue4, present4, fill=0) %>% 
  rename(E4=E, L4=L, O4=O) %>% select(-16)

head(spread_cancer_types)

# clean up spread cancer type information into features
spread_cancer_types_clean <-
  spread_cancer_types %>%
  mutate(
    is_cancer_stringent = 1,
    cancer_s_epithelial = ifelse(E1==1 | E2==1 | E3==1 | E4==1, 1, 0),
    cancer_s_leuk_lymph = ifelse(L1==1 | L2==1 | L3==1 | L4==1, 1, 0),
    cancer_s_mesenchym = ifelse(M1==1 | M2==1 | M3==1, 1, 0),
    cancer_s_other_tiss = ifelse(O1==1 | O2==1 | O3==1 | O4==1, 1, 0)
  ) %>%
  select(-(E1:O4)) %>% 
  mutate(
    count_tissue_types = cancer_s_epithelial + cancer_s_leuk_lymph + 
      cancer_s_mesenchym + cancer_s_other_tiss
  ) %>%
  select(ensg, is_cancer_stringent, somatic, germline, cancer_s_epithelial:count_tissue_types)

head(spread_cancer_types_clean)

# final renaming of table
spread_cancer_types_clean <- spread_cancer_types_clean %>%
  rename(somatic_cancer_mutation = somatic,
         germline_cancer_mutation = germline, 
         epithelial_cancer_mutation = cancer_s_epithelial,
         leuk_lymph_cancer_mutation = cancer_s_leuk_lymph, 
         mesenchymal_cancer_mutation = cancer_s_mesenchym, 
         other_tissue_cancer_mutation= cancer_s_other_tiss, 
         count_tissues_cancer_mutation = count_tissue_types)

# write out
write.table(spread_cancer_types_clean, 
            "~/Projects/fusion_ML/features/cancer_genes_stringent.csv", sep=",", quote = FALSE, row.names = FALSE)



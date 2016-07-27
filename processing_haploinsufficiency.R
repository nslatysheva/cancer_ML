# Haploinsufficient genes

require(dplyr)
require(tidyr)

haplo <- read.table("~/Projects/fusion_ML/data/gene_identity_data/happloinsufficiency/happloinsuf.txt", sep="\t", header = TRUE)
haplo <- haplo %>% rename(HGNC = Gene.ID, HGNC_symbol = Gene.Symbol); haplo$HGNC <- as.character(haplo$HGNC)
head(haplo)

# convert to ensg
ensg_hgnc <- read.table("~/Projects/fusion_ML/data/base_data/ensg_hgnc_conversion_global.txt", sep="\t", header = TRUE)
tail(ensg_hgnc)

# merge and patch in missing values if possible
haplo_ensg <- left_join(haplo, ensg_hgnc, by="HGNC_symbol")
unmapped_haplo_symbols <- haplo_ensg %>% filter(is.na(ensg)) %>% select(HGNC_symbol) 
write.table(unmapped_haplo_symbols, "~/Projects/fusion_ML/data/gene_identity_data/happloinsufficiency/unmunmapped_haplo_symbols.txt",
            row.names = FALSE, quote = FALSE)

# simplify dataset
haplo_ensg <- 
  haplo_ensg %>% 
  mutate(loss_phenotype = ifelse(Loss.phenotype.OMIM.ID=='', 0, 1),
         is_happloinsufficient = ifelse(Haploinsufficiency.Score >=2 & Haploinsufficiency.Score !=40, 1, 0)) %>%
  select(ensg, HGNC_symbol, HGNC.x, Haploinsufficiency.Score, Haploinsufficiency.Description, is_happloinsufficient, loss_phenotype)

head(haplo_ensg)

# write out feature set
write.csv(haplo_ensg %>% select(ensg, is_happloinsufficient, loss_phenotype) %>% na.omit(), 
            "~/Projects/fusion_ML/features/happloinsufficiency_loss_phenotype.csv", quote = FALSE, row.names = FALSE)



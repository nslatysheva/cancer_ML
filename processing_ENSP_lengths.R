# get protein sequence lengths
require(seqinr)
require(dplyr)

# load sequences
ensp_sequences <- read.fasta("~/Projects/fusion_ML/data/base_data/all_protein_sequences_biomart_fasta.txt", seqtype = "AA")
head(ensp_sequences,10)
length(ensp_sequences)

# small loop to output ensp with aa length
ensp_lengths <- matrix(data = NA, nrow=length(ensp_sequences), ncol=2)

for (i in 1:length(ensp_sequences)){
  ensp = getName(ensp_sequences[i])
  ensp_length = getLength(ensp_sequences[i])
  cat(ensp, " ", ensp_length, "\n")
  ensp_lengths[i,1] <- ensp; ensp_lengths[i,2] <- ensp_length
}

# process out length information
ensp_lengths_clean <- as.data.frame(ensp_lengths) %>% rename(ensp=V1, ensp_length=V2) 

write.table(ensp_lengths_clean, 
            "~/Projects/fusion_ML/data/base_data/all_ensp_lengths.csv", sep=",", quote = FALSE, row.names = FALSE)



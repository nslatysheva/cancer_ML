# examining pairing patterns
require(dplyr)
require(tidyr)
require(ggplot2)

# read in feature space, delete some irrelevant labels
feature_space <- read.csv("~/Projects/fusion_ML/features/feature_spaces/feature_space_imputation_completed.csv", header=TRUE, stringsAsFactors = FALSE)
head(feature_space); feature_space <- feature_space %>% select(-c(is_parent, is_5_prime_parent, is_3_prime_parent, is_recurrent_parent, num_fusions))

# read in fusion list and ensg conversion
fusions <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/klijn_gene_fusions.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
head(fusions); fusions <- fusions %>% select(id, entrez_id_3, entrez_id_5) %>% distinct

# convert fuson entrez to ensg
ensg_entrez <- read.table("~/Projects/fusion_ML/data/base_data/accession_conversion_tables/ensg_entrez_conversion_global.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
head(ensg_entrez)

# label conversion
fusions_ensg <- merge(merge(fusions, ensg_entrez, by.x = "entrez_id_5", by.y = "entrez", all.x = TRUE), 
                      ensg_entrez, by.x = "entrez_id_3", by.y = "entrez", all.x = TRUE) %>% rename(ensg_3 = ensg.y, ensg_5 = ensg.x) %>% distinct()   
head(fusions_ensg)

# Anything we can do to assign ENSGs to the remaining entrez?
missing_5 <- filter(fusions_ensg, is.na(ensg_5)) %>% select(entrez_id_5)
missing_3 <- filter(fusions_ensg, is.na(ensg_3)) %>% select(entrez_id_3)
missing <- sort(unique(c(missing_5$entrez_id_5, missing_3$entrez_id_3))); missing
write.csv(missing, "~/Projects/fusion_ML/analyses/missing_fusion_parent_entrez_ensg.txt", quote = FALSE, row.names = FALSE)

# prepare clean fusion sets
fusions_clean <- fusions_ensg %>% na.omit %>% select(ensg_5, ensg_3) %>% distinct %>% arrange(ensg_5)

# data dimensions
head(fusions_clean); nrow(fusions_clean); str(fusions_clean)
fusions_clean %>% summarize_each(funs(n_distinct))

# construct matrix of all possible fusions between parents
fusion_matrix <- fusions_clean %>% mutate(present=1) %>% spread(ensg_3, present, fill=0)
head(fusion_matrix)[,1:10]; dim(fusion_matrix)

fusion_matrix %>% filter(ensg_5 == "ENSG00000000971")

# get percentage sparsity
# this will get rounded up to 100 sometimes
percentage_sparsity <- 100*length(fusion_matrix[fusion_matrix==0])/(length(fusion_matrix[fusion_matrix==0]) + length(fusion_matrix[fusion_matrix==1]))
format(percentage_sparsity, digits = 5)

##### Generate random parent sets
n_random_pairs = nrow(fusions_clean)*1.20; n_random_pairs

generate_random_fusions <- function(fusion_product_matrix = fusion_matrix, n_random_pairs=n_random_pairs){
  # make matrix containing randomly generated parent sets
  random_parent_sets <- matrix(NA, nrow=n_random_pairs, ncol=2)
  
  for(i in 1:n_random_pairs){
    # sample random 5' and 3' parent
    # correct for offset by 1 on 3' list
    five_prime_parent_index <- sample(dim(fusion_matrix)[1], 1); print(five_prime_parent_index)
    three_prime_parent_index <- sample(2:dim(fusion_matrix)[2], 1); print(three_prime_parent_index)
  
    # check these aren't actually a real fusion pair
    if(fusion_matrix[five_prime_parent_index,three_prime_parent_index] == 1) { print("warning: actual fusion pair; skipping"); next }
    
    # if all okay, that's a new fusion pair
    five_prime_parent <- as.character(fusion_matrix[five_prime_parent_index,1])
    three_prime_parent <- colnames(fusion_matrix)[three_prime_parent_index]
    random_parent_sets[i,] <- c(five_prime_parent, three_prime_parent)
  }
  # format random parent sets and output
  random_parent_sets <- as.data.frame(random_parent_sets)
  colnames(random_parent_sets) <- c("ensg_5", "ensg_3")
  return(random_parent_sets)
}

random_parent_sets <- generate_random_fusions(fusion_product_matrix = fusion_matrix, n_random_pairs = n_random_pairs)
head(random_parent_sets); nrow(random_parent_sets); nrow(fusions_clean)

# given a list of two parents, generate paired feature space
# skips over rows where either ENSG is not in feature space
head(feature_space)
fill_paired_feature_space <- function(feature_set = feature_space, gene_list){
  
  # create output matrix 
  output <- matrix(NA, nrow=nrow(gene_list), ncol=ncol(feature_space)*2)

  # loop over gene list to pair feature sets
  cat("Merging ...")
  merged <- merge(merge(gene_list, feature_space, by.x = "ensg_5", by.y = "ensg"), feature_space, by.x = "ensg_3", by.y = "ensg")
  cat("Outputting ...")
  output <- merged %>% select(ensg_5, is_oncogene.x:density_phospho_sites.x, ensg_3, is_oncogene.y:density_phospho_sites.y)
  
#   for (i in 1:nrow(gene_list)) {
#     # progress bar
#     print(paste(i, "/", nrow(gene_list), "...", sep=""))
#     
#     # get feature vectors if they exist; otherwise, go to next gene pair
#     # if(gene_list$ensg_5[i] %in% feature_space$ensg & gene_list$ensg_3[i] %in% feature_space$ensg){
#       five_feature_vector <- feature_space %>% filter(ensg==gene_list$ensg_5[i])
#       three_feature_vector <- feature_space %>% filter(ensg==gene_list$ensg_3[i])
#     # }
#     # else {print(paste("skipping gene pair", gene_list$ensg_5[i], "and", gene_list$ensg_3[i], sep=" ")); next}
#     
#     # combine feature vectors and print out
#     row <- cbind(five_feature_vector, three_feature_vector)
#     names(row) <- NULL
#     output[i,] <- as.matrix(row)
#   }
#   
#   # format and output paired features
#   print(paste("nrow before na.omit: ", nrow(output)))
#   output <- output %>% na.omit
#   print(paste("nrow after na.omit: ", nrow(output)))
  # data framing and colnames
  names_5 <- colnames(feature_space) %>% paste("_5_prime", sep="")
  names_3 <- colnames(feature_space) %>% paste("_3_prime", sep="")
  output <- as.data.frame(output); colnames(output) <- c(names_5, names_3)
  # return
  return(output) 
}

# fill feature spaces, clean results and convert to data frame
# fill clean fusion set
paired_features_fusions <- fill_paired_feature_space(gene_list = fusions_clean)
dim(fusions_clean); dim(paired_features_fusions)

# fill randomised fusion set
paired_features_random <- fill_paired_feature_space(gene_list = random_parent_sets)
dim(random_parent_sets); dim(paired_features_random) 
paired_features_random <- head(paired_features_random, nrow(paired_features_fusions))

# double check dimensions
dim(paired_features_fusions); dim(paired_features_random) 

# Combine fusion and random feature pairings
paired_features <- rbind(paired_features_fusions %>% mutate(is_observed_fusion = 1),
                         paired_features_random %>% mutate(is_observed_fusion = 0))

write.csv(paired_features, "~/Projects/fusion_ML/features/feature_spaces/observed_vs_random_paired_features_5_3_ORDER_CONSERVED.csv", quote = FALSE, row.names = FALSE)

##### Generate sets where 5' and 3' order is NOT conserved
# quickly - how often do we get self fusions? Never.
head(fusions_clean); fusions_clean %>% filter(ensg_5==ensg_3)

# matrix of all parents
all_parents <- data.frame(c(fusions_clean$ensg_5, fusions_clean$ensg_3)) %>% distinct; names(all_parents) <- "ensg"
head(all_parents)

# function
generate_random_fusions_mixed <- function(parent_list = all_parents, n_random_pairs=n_random_pairs){
  
  # make matrix containing randomly generated parent sets
  random_parent_sets <- matrix(NA, nrow=n_random_pairs, ncol=2)
  
  for(i in 1:n_random_pairs){
    # sample random 5' and 3' parent
    five_prime_parent_index <- sample(dim(parent_list)[1], 1);  #print(five_prime_parent_index)
    three_prime_parent_index <- sample(dim(parent_list)[1], 1); #print(three_prime_parent_index)
    
    # check these aren't the same one
    if(five_prime_parent_index == three_prime_parent_index) { print("warning: same gene chosen for fusion pair; skipping"); next }
    
    # if all okay, that's a new fusion pair
    five_prime_parent <- as.character(parent_list[five_prime_parent_index,1]); #print(five_prime_parent)
    three_prime_parent <- as.character(parent_list[three_prime_parent_index,1]); #print(three_prime_parent)
    random_parent_sets[i,] <- c(five_prime_parent, three_prime_parent)
  }
  # format random parent sets and output
  random_parent_sets <- as.data.frame(random_parent_sets)
  colnames(random_parent_sets) <- c("ensg_5", "ensg_3")
  return(random_parent_sets)
}

random_parent_sets_mixed <- generate_random_fusions_mixed(parent_list = all_parents, n_random_pairs = n_random_pairs)
head(random_parent_sets_mixed)

# fill feature spaces, clean results and convert to data frame
# take just nrow from random set
paired_features_random_mixed <- fill_paired_feature_space(gene_list = random_parent_sets_mixed)
dim(random_parent_sets_mixed); dim(paired_features_random_mixed)

paired_features_random_mixed <- head(paired_features_random_mixed, nrow(paired_features_fusions))


# Combine fusion and random feature pairings
paired_features_mixed <- rbind(paired_features_fusions %>% mutate(is_observed_fusion = 1),
                               paired_features_random_mixed %>% mutate(is_observed_fusion = 0))

write.csv(paired_features_mixed, "~/Projects/fusion_ML/features/feature_spaces/observed_vs_random_paired_features_ORDER_NOT_CONSERVED.csv", quote = FALSE, row.names = FALSE)




#####
#####  Generate *all possible* gene fusion events, to be tested on with trained models
all_ensg <- read.table("~/Projects/fusion_ML/features/ensg_basic_structural_info.csv", sep=",", header=TRUE) %>% select(ensg) %>% distinct
head(all_ensg); dim(all_ensg)

# generate matrix of all possible fusions
# command takes a while
all_possible_fusions <- all_ensg %>% mutate(ensg_2 = ensg, fusion = 1) %>% spread(key = ensg_2, value = fusion, fill=1)
head(all_possible_fusions)[1:10]; dim(all_possible_fusions)
n_possible_fusions = dim(all_possible_fusions)[1] * (dim(all_possible_fusions)[2]-1); n_possible_fusions
  
  # generate all possible pairs
  all_pairs <- all_possible_fusions %>% 
    gather(partner, del, -ensg) %>% 
    select(-del) %>% rename(ensg_5 = ensg, ensg_3 = partner); head(all_pairs) 

# break up calculation into batches
batch_size = n_possible_fusions/200; batch_size

batches <- c(seq(from = 1, to = n_possible_fusions+1, by = batch_size)); batches

for(i in 1:(length(batches)-1)){
  cat("\nProcessing batch", i, "...")
  # get features
  current_batch <-  all_pairs[batches[i]:batches[i+1], ]; dim(current_batch)
  current_batch_with_features <- fill_paired_feature_space(gene_list = current_batch)

  # write out
  i=1
  filename = paste("~/Projects/fusion_ML/analyses/all_possible_fusions/all_possible_fusion_features", "_batch_", i, "_n_", batches[i], "_to_", batches[i+1], ".csv", sep=""); filename
  write.csv(current_batch_with_features, file = filename, quote = FALSE, row.names = FALSE)
}

# to combine the 200 batches, in shell do:
# cat current_batch_* > all_possible_fusion_features.csv


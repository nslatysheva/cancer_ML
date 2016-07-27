# Feature engineering and aggregation - parent gene level
# ptms shouldn't have avg_length count_ENSPs'
require(dplyr)
require(tidyr)
require(mice) # for imputation (mice=multivariate imputation by chained equations)
require(unbalanced) # for synthetic majority oversampling technique
theme_set(theme_bw(base_size = 20))

#### 1) Load datasets
# Fusion status + recurrency
parent_genes <- read.table("~/Projects/fusion_ML/features/parent_genes.csv", header=TRUE, sep=','); head(parent_genes); 
five_three_prime_parents <- read.table("~/Projects/fusion_ML/features/five_three_prime_parents.csv", header=TRUE, sep=','); head(five_three_prime_parents)
recurrent_parent_genes <- read.table("~/Projects/fusion_ML/features/recurrent_parent_genes.csv", header=TRUE, sep=','); head(recurrent_parent_genes)

## Gene identity
oncogene <- read.table("~/Projects/fusion_ML/features/oncogene_list.csv", header=TRUE, sep=','); head(oncogene)
tsg <- read.table("~/Projects/fusion_ML/features/tsg_genes.csv", header=TRUE, sep=','); names(tsg)[3] <- "num_cancers_downreg_TSG";  head(tsg);
cancer_genes <- read.table("~/Projects/fusion_ML/features/cancer_genes.csv", header=TRUE, sep=','); head(cancer_genes)
cancer_genes_stringent <- read.table("~/Projects/fusion_ML/features/cancer_genes_stringent.csv", header=TRUE, sep=','); head(cancer_genes_stringent)

multifunctionality <- read.table("~/Projects/fusion_ML/features/GO_counts.csv", header=TRUE, sep=','); head(multifunctionality)
essential_genes <- read.table("~/Projects/fusion_ML/features/essential_genes.csv", header=TRUE, sep=','); head(essential_genes)
haplo_loss_pheno <- read.table("~/Projects/fusion_ML/features/happloinsufficiency_loss_phenotype.csv", header=TRUE, sep=','); head(haplo_loss_pheno)

kinase_genes <- read.table("~/Projects/fusion_ML/features/kinase_genes.csv", header=TRUE, sep=','); head(kinase_genes)
kinase_classes <- read.table("~/Projects/fusion_ML/features/kinase_classes.csv", header=TRUE, sep=','); head(kinase_classes)
tf_genes <- read.table("~/Projects/fusion_ML/features/tf_genes.csv", header=TRUE, sep=','); head(tf_genes)
tf_classes <- read.table("~/Projects/fusion_ML/features/tf_classes.csv", header=TRUE, sep=','); head(tf_classes)

epigenetic_genes <- read.table("~/Projects/fusion_ML/features/epigenetic_gene.csv", header=TRUE, sep=',');  colnames(epigenetic_genes)[2:4] <- paste("is_", colnames(epigenetic_genes)[2:4], sep=""); head(epigenetic_genes)
epigenetic_complexes <- read.table("~/Projects/fusion_ML/features/epigenetic_complexes.csv", header=TRUE, sep=','); head(epigenetic_complexes)

## Sub-cellular localization
localization <- read.table("~/Projects/fusion_ML/features/subcellular_localization.csv", header=TRUE, sep=','); head(localization)

## Basic structural information
basic_info <- read.table("~/Projects/fusion_ML/features/ensg_basic_structural_info.csv", header=TRUE, sep=','); head(basic_info)
pfam <- read.table("~/Projects/fusion_ML/features/pfam.csv", header=TRUE, sep=','); head(pfam)
disorder <- read.table("~/Projects/fusion_ML/features/avg_disorder_by_ensg.csv", header=TRUE, sep=','); head(disorder)

## Structural features + PTMs
interacting_segments <- read.table("~/Projects/fusion_ML/features/interacting_tracts.csv", header=TRUE, sep=','); colnames(interacting_segments)[2:3] <- c("longest_interacting_segment", "num_interacting_segments"); head(interacting_segments)
instruct_domains <- read.table("~/Projects/fusion_ML/features/instruct_domains.csv", header=TRUE, sep=','); colnames(instruct_domains)[2] <- "num_INstruct_domains"; head(instruct_domains)
pisa_res <- read.table("~/Projects/fusion_ML/features/PISA_residues.csv", header=TRUE, sep=','); colnames(pisa_res)[2] <- "num_PISA_res"; head(pisa_res)
ELM_lms <- read.table("~/Projects/fusion_ML/features/lm_features_ELM.csv", header=TRUE, sep=',');  colnames(ELM_lms)[2] <- "num_ELM_LMs"; head(ELM_lms)
ANCHOR_lms <- read.table("~/Projects/fusion_ML/features/ANCHOR_lms.csv", header=TRUE, sep=','); colnames(ANCHOR_lms)[2] <- "num_ANCHOR_LMs"; head(ANCHOR_lms)
ptms <- read.table("~/Projects/fusion_ML/features/ptms.csv", header=TRUE, sep=','); names(ptms) <- c("ensg", "num_PTMs", "num_UB_sites", "num_phospho_sites"); head(ptms)
ptm_code <- read.table("~/Projects/fusion_ML/features/ptmcode.csv", header=TRUE, sep=','); names(ptm_code) <- c("ensg", "num_PTMcode_sites"); head(ptm_code)

# Molecular interactions
cancer_pathway <- read.table("~/Projects/fusion_ML/features/KEGG_cancer_pathway_genes.csv", header=TRUE, sep=','); head(cancer_pathway)
interacts_with_cancer_gene <- read.table("~/Projects/fusion_ML/features/interacts_with_cancer.csv", header=TRUE, sep=','); interacts_with_cancer_gene <- interacts_with_cancer_gene %>% select(ensg, interacts_with_cancer, num_cancer_interactions); head(interacts_with_cancer_gene)
interacts_with_OG <- read.table("~/Projects/fusion_ML/features/interacts_with_OG.csv", header=TRUE, sep=','); interacts_with_OG <- interacts_with_OG %>% select(ensg, interacts_with_OG, num_OG_interactions);  head(interacts_with_OG)
interacts_with_TSG <- read.table("~/Projects/fusion_ML/features/interacts_with_TSG.csv", header=TRUE, sep=','); interacts_with_TSG <- interacts_with_TSG %>% select(ensg, interacts_with_TSG, num_TSG_interactions); head(interacts_with_TSG)

# Network centrality
centrality <- read.table("~/Projects/fusion_ML/features/centrality_biogrid.csv", header=TRUE, sep=','); head(centrality)

# Expression features
expression <- read.table("~/Projects/fusion_ML/features/expression_and_breadth.csv", header=TRUE, sep=','); head(expression)

# Merge datasets
feature_space <-
  list(basic_info %>% select(ensg),
       parent_genes,
       five_three_prime_parents,
       recurrent_parent_genes,
       oncogene,
       tsg,
       cancer_genes,
       cancer_genes_stringent,
       multifunctionality,
       essential_genes,
       haplo_loss_pheno,
       kinase_genes, 
       kinase_classes,
       tf_genes,
       tf_classes, 
       epigenetic_genes,
       epigenetic_complexes,
       localization,
       basic_info, 
       pfam,
       disorder,
       interacting_segments,
       instruct_domains,
       pisa_res,
       ELM_lms,
       ANCHOR_lms,
       ptms,
       ptm_code,
       cancer_pathway,
       interacts_with_cancer_gene,
       interacts_with_OG,
       interacts_with_TSG, 
       centrality,        
       expression) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="ensg"), .)

head(feature_space,5); ncol(feature_space); nrow(feature_space); str(feature_space)

### 3) NA value handling
# Some missing values, like the gene identify labels, are implied 0 if missing
variables_to_0 <- c("is_parent", "is_5_prime_parent", "is_3_prime_parent", "is_recurrent_parent", "num_fusions", 
                    "is_oncogene",	"is_TSG",	"num_cancers_downreg_TSG", "is_cancer",	"is_cancer_stringent",	
                    "somatic_cancer_mutation",	"germline_cancer_mutation",	"epithelial_cancer_mutation",	"leuk_lymph_cancer_mutation",
                    "mesenchymal_cancer_mutation", "other_tissue_cancer_mutation","count_tissues_cancer_mutation", 
                    "is_essential",	"is_happloinsufficient", "loss_phenotype",	
                    "is_kinase",	"kinase_AGC",	"kinase_Atypical",	"kinase_CAMK",	"kinase_CK1",
                    "kinase_CMGC",	"kinase_NEK",	"kinase_Other",	"kinase_RGC",	"kinase_STE",	"kinase_TKL",	"kinase_Tyr",
                    "is_TF",	"num_TF_interactions",	"num_TF_activation",	"num_TF_repression", 
                    "is_chromatin_remodelling", "is_histone_modification", "is_other_epigenetic", "is_epigenetic_complex", 
                    "loc_cell_junction", "loc_cytoplasm", "loc_cytoskeleton", "loc_endomembrane_system", "loc_extracellular_region",
                    "loc_nucleolus", "loc_nucleus", "loc_organelle_lumen", "loc_plasma_membrane", "loc_vesicle",
                    "avg_num_domains", "avg_domain_length",
                    "longest_interacting_segment", "num_interacting_segments",
                    "num_INstruct_domains", "num_PISA_res",	"num_ELM_LMs",	"num_ANCHOR_LMs",	
                    "num_PTMs",	"num_UB_sites",	"num_phospho_sites",	"num_PTMcode_sites",
                    "is_KEGG_cancer_pathway", "interacts_with_cancer", "num_cancer_interactions",	
                    "interacts_with_OG", "num_OG_interactions",	"interacts_with_TSG",	"num_TSG_interactions")

# replace NA values in these columns with 0
feature_space[,variables_to_0][is.na(feature_space[,variables_to_0])] <- 0
head(feature_space)

### 4) Calculating densities
# calculate densities of molecular features
feature_space <- feature_space %>% mutate(
  density_domains = avg_num_domains/avg_ensp_length,
  density_INstruct_domains = num_INstruct_domains/avg_ensp_length,
  density_PISA_res = num_PISA_res/avg_ensp_length,
  density_ELM_LMs = num_ELM_LMs/avg_ensp_length,
  density_ANCHOR_LMs = num_ANCHOR_LMs/avg_ensp_length,
  density_PTMs = num_PTMs/avg_ensp_length,
  density_PTMcode_sites = num_PTMcode_sites/avg_ensp_length,
  density_UB_sites = num_UB_sites/avg_ensp_length,
  density_phospho_sites = num_phospho_sites/avg_ensp_length
)

# nrow with/without distinct
nrow(basic_info); nrow(feature_space); nrow(feature_space %>% distinct)

# are there any parent genes not labeled as either 5’ or 3’?
feature_space %>% filter(is_parent == 1 & is_5_prime_parent == 0 & is_3_prime_parent == 0) 
# any other wacky things?
feature_space %>% filter(is_parent == 0 & (is_5_prime_parent == 1 | is_3_prime_parent == 1))
feature_space %>% filter(is_parent == 0 & is_recurrent_parent == 1)

# any ENSGs in multiple times?
feature_space %>% group_by(ensg) %>% summarise(count=n()) %>% arrange(-count)
# if distinct() first?
feature_space %>% distinct %>% group_by(ensg) %>% summarise(count=n()) %>% arrange(-count)

# deduplicate to get rid of some duplicated entries
feature_space <- feature_space %>% distinct

# some statistics on the integrated dataset
nrow(feature_space %>% filter(is_parent == 1) %>% select(ensg) %>% distinct)
nrow(feature_space %>% filter(is_recurrent_parent == 1) %>% select(ensg) %>% distinct)

# Table S2, write out raw dataset if people want it
write.csv(feature_space, "~/Projects/fusion_ML/features/feature_spaces/raw_feature_space_before_imputation_balancing_simplifying.csv", row.names = FALSE, quote = FALSE)

# check columns/feature counts
ncol(feature_space); head(feature_space)
names(feature_space)

################## Further cleaning
# set up data set, get rid of extra expression information
head(feature_space)
feature_space_simplified <- feature_space %>% 
  select(-c(Averaged.RPKM.colon:Averaged.RPKM.sgland, Gini:Pem)) 

# identify factor variables
head(feature_space_simplified)
factor_variables <- c("is_parent", "is_5_prime_parent", "is_3_prime_parent", "is_recurrent_parent", 
                      "is_oncogene",	"is_TSG",	"is_KEGG_cancer_pathway", "is_cancer",	"is_cancer_stringent",	
                      "somatic_cancer_mutation",	"germline_cancer_mutation",	"epithelial_cancer_mutation",	"leuk_lymph_cancer_mutation",
                      "mesenchymal_cancer_mutation", "other_tissue_cancer_mutation",
                      "is_essential",	"is_happloinsufficient","loss_phenotype",	
                      "is_kinase",	"kinase_AGC",	"kinase_Atypical",	"kinase_CAMK",	"kinase_CK1",
                      "kinase_CMGC",	"kinase_NEK",	"kinase_Other",	"kinase_RGC",	"kinase_STE",	"kinase_TKL",	"kinase_Tyr",
                      "is_TF", "is_chromatin_remodelling", "is_histone_modification", "is_other_epigenetic", "is_epigenetic_complex", 
                      "loc_cell_junction", "loc_cytoplasm", "loc_cytoskeleton", "loc_endomembrane_system", "loc_extracellular_region",
                      "loc_nucleolus", "loc_nucleus", "loc_organelle_lumen", "loc_plasma_membrane", "loc_vesicle",
                      "interacts_with_cancer",	"interacts_with_OG", "interacts_with_TSG")

integer_variables <- c("num_fusions", "num_isoforms", "num_cancers_downreg_TSG", "count_tissues_cancer_mutation", 
                       "num_GO_terms", "num_GOSlim_terms", "num_TF_interactions", "num_TF_activation", "num_TF_repression",
                       "longest_interacting_segment", "num_interacting_segments", "num_INstruct_domains", "num_PISA_res", "num_ELM_LMs", "num_ANCHOR_LMs", "num_PTMs", "num_UB_sites", "num_phospho_sites", "num_PTMcode_sites",
                       "degree_centrality", "num_cancer_interactions", "num_OG_interactions", "num_TSG_interactions")

feature_space_simplified_typed <- 
  feature_space_simplified %>%
  mutate_each_(funs(as.factor), factor_variables) %>%
  mutate_each_(funs(as.integer), integer_variables)

# get some simple counts
# parents
nrow(feature_space_simplified_typed %>% filter(is_parent==1)); nrow(feature_space_simplified_typed %>% filter(is_parent==1) %>% select(ensg) %>% distinct)
# non parents
nrow(feature_space_simplified_typed %>% filter(is_parent==0)); nrow(feature_space_simplified_typed %>% filter(is_parent==0) %>% select(ensg) %>% distinct)
# recurrent
nrow(feature_space_simplified_typed %>% filter(is_recurrent_parent==1)); nrow(feature_space_simplified_typed %>% filter(is_parent==1) %>% filter(is_recurrent_parent==1) %>% select(ensg) %>% distinct)
nrow(feature_space_simplified_typed %>% filter(is_recurrent_parent==1))/nrow(feature_space_simplified_typed %>% filter(is_parent==1))
# 5' parents
nrow(feature_space_simplified_typed %>% filter(is_5_prime_parent==1)); nrow(feature_space_simplified_typed %>% filter(is_parent==1) %>% filter(is_5_prime_parent==1) %>% select(ensg) %>% distinct)
nrow(feature_space_simplified_typed %>% filter(is_5_prime_parent==1))/nrow(feature_space_simplified_typed %>% filter(is_parent==1))
# 3' parents
nrow(feature_space_simplified_typed %>% filter(is_3_prime_parent==1)); nrow(feature_space_simplified_typed %>% filter(is_parent==1) %>% filter(is_3_prime_parent==1) %>% select(ensg) %>% distinct)
nrow(feature_space_simplified_typed %>% filter(is_3_prime_parent==1))/nrow(feature_space_simplified_typed %>% filter(is_parent==1))


####### Imputation using mice package ####### 
# mice = multivariate imputation by chained equations
# m = number of multiple imputations
# maxit = scalar giving number of iterations
# method, pmm = predictive mean matching
count_na <- function(data) sum(is.na(data)) 
feature_space_simplified_typed %>% summarize_each(funs(count_na))

# Previously, would not impute betweenness or closeness centrality
# scaling betweenness centrality allows it to work,  scale(test$betweenness_centrality, center=TRUE)
# testing if normalizing these scores in the biogrid centrality script will help - it does :)
feature_space_imputed <- mice(feature_space_simplified_typed,
                              m=5,maxit=5,meth='pmm',seed=1)

# examine predicted values
feature_space_imputed$imp
feature_space_imputed$imp$degree
feature_space_imputed$imp$num_GO_terms

# Check imputed data
head(feature_space_imputed$imp$degree)

# Complete original dataset with imputed data
feature_space_with_imputation <- complete(feature_space_imputed,1)

# check NA values again
feature_space_with_imputation %>% summarize_each(funs(count_na))

###### write out this feature space (note: NOT BALANCED) ###### 
write.csv(feature_space_with_imputation, "~/Projects/fusion_ML/features/feature_spaces/feature_space_imputation_completed.csv", row.names = FALSE, quote = FALSE)

###### write out this dataset for 5' vs 3' analysis, don't need to balance this ###### 
head(feature_space_with_imputation %>% filter(is_parent==1) %>% select(-c(is_parent, is_recurrent_parent, num_fusions)))

# how many genes for 5' and 3' parents?
nrow(feature_space_with_imputation %>% filter(is_5_prime_parent==1))
nrow(feature_space_with_imputation %>% filter(is_3_prime_parent==1))
nrow(feature_space_with_imputation %>% filter(is_5_prime_parent==1 & is_3_prime_parent==1))

# export dataset
# rename categories as integers for scikit's case
# 0 = 5' parent
# 1 = 3' parent
# 2 = 5' and 3' parent
feature_space_5_3_prime_parents <- feature_space_with_imputation %>% filter(is_parent==1) %>%
  mutate(parent_category = ifelse(is_5_prime_parent == 1 & is_3_prime_parent == 1, "2",
                                  ifelse(is_5_prime_parent == 0, "0", "1"))) %>% 
  select(-c(is_parent, is_recurrent_parent, num_fusions, is_5_prime_parent, is_3_prime_parent))

write.csv(feature_space_5_3_prime_parents, "~/Projects/fusion_ML/features/feature_spaces/five_and_three_prime_parents.csv", row.names = FALSE, quote = FALSE)
nrow(feature_space_5_3_prime_parents); nrow(feature_space_with_imputation %>% filter(is_parent==1))

# And if we're only interested in genes that are either 5' or 3', not both
write.csv(feature_space_5_3_prime_parents %>% filter(parent_category != 2), "~/Projects/fusion_ML/features/feature_spaces/five_and_three_prime_parents_not_both.csv", row.names = FALSE, quote = FALSE)
nrow(feature_space_5_3_prime_parents %>% filter(parent_category != 2))

################### Write out downsampled datasets 
# need to downsample majority cases
# parents vs. non parents
summary(feature_space_with_imputation$is_parent)
num_minority <- summary(feature_space_with_imputation$is_parent)[[2]]; num_minority

# downsample majority case
majority_cases <- feature_space_with_imputation %>% filter(is_parent == 0)
downsampled_majority_indices <- sample(1:nrow(majority_cases), size = num_minority, replace = FALSE)

# output
parent_non_parent_downsampled_majority <- rbind(feature_space_with_imputation %>% filter(is_parent == 1),
                                                majority_cases[downsampled_majority_indices, ])

write.csv(parent_non_parent_downsampled_majority, "~/Projects/fusion_ML/features/feature_spaces/downsampled_majority_parent_non_parent_feature_set.csv", row.names = FALSE, quote = FALSE)

# baseline: predicting fusion parents as cancer genes
baseline <- parent_non_parent_downsampled_majority %>% select(is_parent, is_cancer, is_cancer_stringent)
head(baseline)
nrow(baseline %>% filter(is_parent == is_cancer))/nrow(baseline)
nrow(baseline %>% filter(is_parent == is_cancer_stringent))/nrow(baseline)

# random chance
nrow(parent_non_parent_downsampled_majority %>% filter(is_parent==1))/nrow(parent_non_parent_downsampled_majority)

#### recurrent vs. non recurrent
parent_only_data <- feature_space_with_imputation %>% filter(is_parent==1)
summary(parent_only_data$is_recurrent_parent)
summary(parent_only_data$is_recurrent_parent)[[2]]/nrow(parent_only_data)

num_minority <- summary(parent_only_data$is_recurrent_parent)[[2]]; num_minority

# downsample majority case
majority_cases <- parent_only_data %>% filter(is_recurrent_parent == 0)
downsampled_majority_indices <- sample(1:nrow(majority_cases), size = num_minority, replace = FALSE)

# output
recurrent_non_recurrent_downsampled_majority <- rbind(parent_only_data %>% filter(is_recurrent_parent == 1),
                                                majority_cases[downsampled_majority_indices, ])

write.csv(recurrent_non_recurrent_downsampled_majority, "~/Projects/fusion_ML/features/feature_spaces/downsampled_majority_recurrent_non_recurrent_feature_set.csv", row.names = FALSE, quote = FALSE)

########## Generating dataset for fusions, not just parents
fusion_events_ensg_cancer_type_is_metastatic <- read.table("~/Projects/fusion_ML/features/fusion_events_with_ensg_cancer_type_is_metastatic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(fusion_events_ensg_cancer_type_is_metastatic); str(fusion_events_ensg_cancer_type_is_metastatic)

### First, clean dataset
fusion_events <- fusion_events_ensg_cancer_type_is_metastatic %>% select(fusion_id, ensg_5, ensg_3, cancer_type, is_metastatic) %>% distinct %>% arrange(fusion_id)
head(fusion_events); nrow(fusion_events)

# overlap feature set, feature_space_with_imputation, to get features of fusion events
str(feature_space_with_imputation)
fusion_events_with_features <- merge(merge(fusion_events, feature_space_with_imputation, by.x="ensg_5", by.y="ensg", all.x = FALSE),
                                     feature_space_with_imputation, by.x="ensg_3", by.y="ensg", all.x=FALSE)
head(fusion_events_with_features); nrow(fusion_events_with_features)

# fix column names
fusion_events_with_features <- fusion_events_with_features %>% select(-c(is_parent.x, is_parent.y, 
                                       is_5_prime_parent.x, is_5_prime_parent.y, 
                                       is_3_prime_parent.x, is_3_prime_parent.y,
                                       is_recurrent_parent.x, is_recurrent_parent.y,
                                       num_fusions.x, num_fusions.y))

names(fusion_events_with_features) <- names(fusion_events_with_features) %>% gsub(pattern = "\\.x", replacement = "_5_prime") %>% gsub(pattern = "\\.y", replacement = "_3_prime")

# print out this dataset, random pairing 3' 5' stuff
write.csv(fusion_events_with_features, "~/Projects/fusion_ML/features/feature_spaces/fusion_events_with_features_original.csv", row.names = FALSE, quote = FALSE)

#### summary statistics of fusion events
head(fusion_events_with_features)[1:10]
str(fusion_events_with_features)

# compare recurrent + metastatic parent genes overlap
# get metastatic parent set
metastatic_fusions <- fusion_events_with_features %>% filter(is_metastatic == 1)
metastatic_parent_genes <- unique(c(as.character(metastatic_fusions$ensg_5), as.character(metastatic_fusions$ensg_3)))
head(metastatic_parent_genes); length(metastatic_parent_genes)

# get recurrent parent set
recurrent_parent_genes <- (feature_space_with_imputation %>% filter(is_recurrent_parent == 1))$ensg; length(recurrent_parent_genes)
nrow(feature_space_with_imputation %>% filter(is_recurrent_parent == 1) %>% select(ensg) %>% distinct)

# look at feature set for this
feature_space_with_imputation <- feature_space_with_imputation %>% mutate(is_metastatic_parent = ifelse(ensg %in% metastatic_parent_genes, 1, 0))
summary(as.factor(feature_space_with_imputation$is_metastatic_parent))

# looks like recurrent genes are more likely to be metastatic
m <- matrix(c(feature_space_with_imputation %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 1) %>% nrow(), 
              feature_space_with_imputation %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 0) %>% nrow(), 
              feature_space_with_imputation %>% filter(is_recurrent_parent == 0 & is_metastatic_parent== 1) %>% nrow(), 
              feature_space_with_imputation %>% filter(is_parent == 1 & is_recurrent_parent == 0 & is_metastatic_parent== 0) %>% nrow()), nrow=2, byrow = TRUE); m

chisq.test(m)
fisher.test(m)

# check total counts make sense
feature_space_with_imputation %>% filter(is_recurrent_parent == 1) %>% nrow()
feature_space_with_imputation %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 1) %>% nrow() + feature_space_with_imputation %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 0) %>% nrow()
#
feature_space_with_imputation %>% filter(is_metastatic_parent == 1) %>% nrow()
feature_space_with_imputation %>% filter(is_recurrent_parent == 0 & is_metastatic_parent== 1) %>% nrow() + feature_space_with_imputation %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 1) %>% nrow()
# 
feature_space_with_imputation %>% filter(is_parent == 1) %>% nrow()


####### balanced metastasis/cancer type sets ####### 
head(fusion_events_with_features)

# balance just with downsampling
num_minority <- summary(as.factor(fusion_events_with_features$is_metastatic))[[2]]; num_minority

# downsample majority case
majority_cases <- fusion_events_with_features %>% filter(is_metastatic == 0)
downsampled_majority_indices <- sample(1:nrow(majority_cases), size = num_minority, replace = FALSE)

# output
fusion_events_downsampled_majority_metastasis <- rbind(fusion_events_with_features %>% filter(is_metastatic == 1),
                                                       majority_cases[downsampled_majority_indices, ])
write.csv(fusion_events_downsampled_majority_metastasis, "~/Projects/fusion_ML/features/feature_spaces/fusion_events_downsampled_majority_for_metastasis.csv", row.names = FALSE, quote = FALSE)


####### balanced cancer_type set
fusion_events_with_features %>% group_by(cancer_type) %>% summarise(count=n()) %>% arrange(-count)
nrow(fusion_events_with_features)

small_sample_size_cancers <- fusion_events_with_features %>% 
  group_by(cancer_type) %>% summarise(count=n()) %>% arrange(-count) %>%
  filter(count < 100) %>% select(cancer_type)

# Lower number of categories
fusion_events_with_features$cancer_type <- as.character(fusion_events_with_features$cancer_type)
fusion_events_cancer_broad <- fusion_events_with_features %>%
  mutate(cancer_type_broad = ifelse(cancer_type %in% small_sample_size_cancers$cancer_type, "Other", cancer_type)) %>%
  select(-cancer_type)

fusion_events_cancer_broad %>% group_by(cancer_type_broad) %>% summarise(count=n()) %>% arrange(-count)

# exclude "other"
fusion_events_cancer_broad <- fusion_events_cancer_broad %>% filter(cancer_type_broad != "Other")
fusion_events_cancer_broad %>% group_by(cancer_type_broad) %>% summarise(count=n()) %>% arrange(-count)

# come up with a value to downsample the highest ones to
n_to_downsample_to = 166

# Sample down lung
all_lung_fusions <- fusion_events_cancer_broad %>% filter(cancer_type_broad == "Lung")
subset_of_lung_fusions <- all_lung_fusions[sample(1:nrow(all_lung_fusions), n_to_downsample_to, replace = FALSE),]

# Sample down breast
all_breast_fusions <- fusion_events_cancer_broad %>% filter(cancer_type_broad == "Breast")
subset_of_breast_fusions <- all_breast_fusions[sample(1:nrow(all_breast_fusions), n_to_downsample_to, replace = FALSE),]

# Sample down ovary
all_ovary_fusions <- fusion_events_cancer_broad %>% filter(cancer_type_broad == "Ovary")
subset_of_ovary_fusions <- all_ovary_fusions[sample(1:nrow(all_ovary_fusions), n_to_downsample_to, replace = FALSE),]

# Recombine with other dataset
fusion_events_cancer_broad_downsampled <- rbind(fusion_events_cancer_broad %>% filter(cancer_type_broad != "Lung" & cancer_type_broad != "Breast" & cancer_type_broad != "Ovary"),
                                                subset_of_lung_fusions,
                                                subset_of_breast_fusions,
                                                subset_of_ovary_fusions)

# Compare frequencies
fusion_events_cancer_broad_downsampled %>% group_by(cancer_type_broad) %>% summarise(count=n()) %>% arrange(-count)

# awkward blocks follow
before <- cbind(rownames(data.frame(summary(as.factor(fusion_events_with_features$cancer_type)))),
      data.frame(summary(as.factor(fusion_events_with_features$cancer_type))),
      rep("1) original", length(unique(fusion_events_with_features$cancer_type)))); colnames(before) <- c("cancer_type", "count", "dataset")
   
after <- cbind(rownames(data.frame(summary(as.factor(fusion_events_cancer_broad_downsampled$cancer_type_broad)))),
                data.frame(summary(as.factor(fusion_events_cancer_broad_downsampled$cancer_type_broad))),
               rep("2) balanced", length(unique(fusion_events_cancer_broad_downsampled$cancer_type_broad)))); colnames(after) <- c("cancer_type", "count", "dataset")

         
rbind(before, after)

# make plot and write out results
p <- ggplot(rbind(before, after), aes(x=as.factor(cancer_type), y=count)) + 
  geom_bar(stat="identity", aes(fill=as.factor(cancer_type)), color="darkgrey") + 
  facet_grid(dataset ~ ., scales = "free") +  guides(fill=FALSE) + xlab("") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)); p

ggsave(p, file="~/Projects/fusion_ML/figures/Figure S2 - Class rebalancing example/balanced_cancer_types.pdf", width = 7, height = 8)

# create integer labels for cancer types
fusion_events_cancer_broad_downsampled %>% group_by(cancer_type_broad) %>% summarise(count=n()) %>% arrange(-count)

fusion_events_cancer_broad_downsampled_clean <- fusion_events_cancer_broad_downsampled %>% 
  select(-c(fusion_id, is_metastatic, ensg_3, ensg_5))  %>%
  mutate(cancer_type_broad = plyr::revalue(cancer_type_broad, c("Breast"=1, 
                                                                "Lung"=2,
                                                                "Lymphoid"=3,
                                                                "Ovary"=4,
                                                                "Colo-Rectal"=5,
                                                                "Brain"=6,
                                                                "Head-Neck"=7,
                                                                "Skin"=8)))

write.csv(fusion_events_cancer_broad_downsampled_clean, "~/Projects/fusion_ML/features/feature_spaces/fusion_events_cancer_type_broadened.csv", row.names = FALSE, quote = FALSE)


# regularized logistic regression
# require(penalized)
require(penalized)
require(dplyr)
require(tidyr)
require(caret)

### 1) parent vs. non-parent
feature_space <- read.csv("~/Projects/fusion_ML/features/feature_spaces/downsampled_majority_parent_non_parent_feature_set.csv", header=TRUE); head(feature_space)
str(feature_space); nrow(feature_space)

# get rid of pointless variables
feature_space <- feature_space %>% select(-c(ensg, is_5_prime_parent, is_3_prime_parent, is_recurrent_parent, num_fusions))

# process and scale data
# penalized package doesn't actually like factor variables though
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

feature_space_typed <- feature_space # %>%
#   mutate_each_(funs(as.factor), factor_variables) 
#   # mutate_each_(funs(as.integer), integer_variables); str(feature_space_typed)
#   # mutate_each_(funs(as.integer), c(integer_variables, factor_variables))
# str(feature_space_typed)

# scaling
variables_to_scale <- c("num_cancers_downreg_TSG", "count_tissues_cancer_mutation", "num_GO_terms",  "num_GOSlim_terms", 
                        "num_TF_interactions", "num_TF_activation", "num_TF_repression", 
                        "num_isoforms", "avg_ensp_length", "avg_num_exons", "avg_num_domains", "avg_domain_length", "avg_disorder",
                        "longest_interacting_segment", "num_interacting_segments",
                        "num_INstruct_domains", "num_PISA_res", "num_ELM_LMs", "num_ANCHOR_LMs", "num_PTMs", "num_UB_sites", "num_phospho_sites", "num_PTMcode_sites",
                        "num_cancer_interactions", "num_OG_interactions", "num_TSG_interactions", 
                        "degree_centrality", "betweenness_centrality", "closeness_centrality", "eigenvector_centrality",
                        "tau", "mean_expression", "max_expression", 
                        "density_domains", "density_INstruct_domains", "density_PISA_res",  "density_ELM_LMs", 
                        "density_ANCHOR_LMs", "density_PTMs", "density_PTMcode_sites", "density_UB_sites", "density_phospho_sites")

# scale relevant variables
feature_space_typed_scaled <- cbind(feature_space_typed %>% select(-one_of(variables_to_scale)),
      scale(feature_space_typed %>% select(one_of(variables_to_scale)), center = TRUE, scale = TRUE))
head(feature_space_typed_scaled)

# # TESTING EFFECT OF SCALING
# feature_space_typed_scaled <- cbind(feature_space_typed %>% select(is_parent),
#                                     scale(feature_space_typed %>% select(-is_parent)))
# feature_space_typed_scaled$is_parent <- as.factor(feature_space_typed_scaled$is_parent)

# train test split
inTrain <- createDataPartition(y=feature_space_typed_scaled$is_parent, p=0.70, list=FALSE)
training <- feature_space_typed_scaled[inTrain,]
testing <- feature_space_typed_scaled[-inTrain,]


###### 1) fit model, get feature rankings for parent
str(training)[1:10]
pen <- penalized(response = is_parent, 
                 penalized = select(training, -is_parent),
                 data=training, model = "logistic",
                 lambda1 = 5, steps = 40)

# plot path of coefficients
plotpath(pen, standardize = TRUE, labelsize = 0.4)
plotpath(pen)

# show all models
show(pen)
model_with_10 <- pen[[32]]; model_with_10

# Examine coefficient values
coeffs <- as.data.frame(coefficients(model_with_10)); coeffs
feature_rankings <- as.data.frame(cbind(rownames(coeffs), coeffs[,1])) %>% 
  rename(feature=V1, coefficient=V2) %>%
  arrange(-as.numeric(coefficient)); feature_rankings

write.csv(feature_rankings, "~/Projects/fusion_ML/modelling/logit_feature_ranking_parent_imputed_no_balancing.csv", 
          row.names = FALSE, quote = FALSE)

# prediction
predictions <- predict(model_with_10, penalized=select(testing, -is_parent), data = testing); predictions
prediction_df <- (data.frame(predictions, testing$is_parent)) %>% mutate(prediction_binary = ifelse(predictions>=0.5, 1, 0)); head(prediction_df)
nrow(prediction_df %>% filter(testing.is_parent==prediction_binary))/nrow(prediction_df)

# CV to optimized lambda, and then prediction
cv_opt <- optL1(response = is_parent, 
                penalized = select(training, -is_parent),
                data=training, model = "logistic", fold=10, minlambda1 = 0.5, maxlambda1 = 500)

# predict with best param
pen_cv <- penalized(response = is_parent, 
                 penalized = select(training, -is_parent),
                 data=training, model = "logistic",
                 lambda1 = cv_opt$lambda)
predictions_cv <- predict(pen_cv, penalized=select(testing, -is_parent), data = testing); predictions_cv
prediction_cv_df <- (data.frame(predictions_cv, testing$is_parent)) %>% mutate(prediction_binary = ifelse(predictions_cv>=0.5, 1, 0)); head(prediction_cv_df)
nrow(prediction_cv_df %>% filter(testing.is_parent==prediction_binary))/nrow(prediction_cv_df)
head(data.frame(predictions_cv, testing$is_parent))


########## get predictions out for table
# need ENSGs for test predictions
data.frame(feature_space_typed_scaled$is_parent, feature_space$ensg)
nrow(feature_space)

nrow(feature_space_typed_scaled)
predictions_cv
predictions_cv_training
ensg_list <- feature_space_typed_scaled %>% select()

# oddly, need predictions on training
predictions_cv_training <- predict(pen_cv, penalized=select(training, -is_parent), data = training); predictions_cv_training
prediction_cv_training_df <- (data.frame(predictions_cv_training, training$is_parent)) %>% mutate(prediction_binary = ifelse(predictions_cv_training>=0.5, 1, 0)); head(prediction_cv_training_df)
nrow(prediction_cv_training_df %>% filter(training.is_parent==prediction_binary))/nrow(prediction_cv_training_df)


head(prediction_cv_df)
head(prediction_cv_training_df)
rbind(prediction_cv_training_df, prediction_cv_df)

length(predictions_cv_training)
nrow(training)
nrow(testing)


##### 2) recurrent parents
# set up feature space
feature_space_recurrent <- read.csv("~/Projects/fusion_ML/features/feature_spaces/downsampled_majority_recurrent_non_recurrent_feature_set.csv", header=TRUE)
feature_space_recurrent_typed <- feature_space_recurrent %>% select(-c(ensg, is_parent, is_5_prime_parent, is_3_prime_parent, num_fusions)) 
head(feature_space_recurrent_typed)[1:10]
nrow(feature_space_recurrent_typed)

# scale relevant variables
feature_space_recurrent_typed_scaled <- cbind(feature_space_recurrent_typed %>% select(-one_of(variables_to_scale)),
                                    scale(feature_space_recurrent_typed %>% select(one_of(variables_to_scale)), center = TRUE, scale = TRUE))

# train test split
inTrain <- createDataPartition(y=feature_space_recurrent_typed_scaled$is_recurrent_parent, p=0.70, list=FALSE)
training <- feature_space_recurrent_typed_scaled[inTrain,]
testing <- feature_space_recurrent_typed_scaled[-inTrain,]

pen_recurrent <- penalized(response = is_recurrent_parent, 
                           penalized = select(training, -is_recurrent_parent),
                           data=training, model = "logistic",
                           lambda1 = 5, steps = 40)

# plot path of coefficients
plotpath(pen_recurrent, standardize = TRUE, labelsize = 0.5)
plotpath(pen_recurrent)

# show all models
show(pen_recurrent)
model_with_10 <- pen_recurrent[[32]]; model_with_10

# Examine coefficient values
coeffs <- as.data.frame(coefficients(model_with_10)); coeffs
feature_rankings_recurrent <- as.data.frame(cbind(rownames(coeffs), coeffs[,1])) %>% 
  rename(feature=V1, coefficient=V2) %>%
  arrange(-as.numeric(coefficient)); feature_rankings_recurrent

write.csv(feature_rankings_recurrent, "~/Projects/fusion_ML/modelling/logit_feature_ranking_recurrent_NOT_BALANCED.csv", 
          row.names = FALSE, quote = FALSE)

# prediction
predictions <- predict(model_with_10, penalized=select(testing, -is_recurrent_parent), data = testing); predictions
prediction_df <- (data.frame(predictions, testing$is_recurrent_parent)) %>% mutate(prediction_binary = ifelse(predictions>=0.5, 1, 0)); head(prediction_df)
nrow(prediction_df %>% filter(testing.is_recurrent_parent==prediction_binary))/nrow(prediction_df)

# CV to optimized lambda, and then prediction
cv_opt <- optL1(response = is_recurrent_parent, 
                penalized = select(training, -is_recurrent_parent),
                data=training, model = "logistic", 
                fold=10, minlambda1 = 0.5, maxlambda1 = 150)

# predict with best param
pen_cv <- penalized(response = is_recurrent_parent, 
                    penalized = select(training, -is_recurrent_parent),
                    data=training, model = "logistic",
                    lambda1 = cv_opt$lambda)

predictions_cv <- predict(pen_cv, penalized=select(testing, -is_recurrent_parent), data = testing); predictions_cv
prediction_cv_df <- (data.frame(predictions_cv, testing$is_recurrent_parent)) %>% mutate(prediction_binary = ifelse(predictions_cv>=0.5, 1, 0)); head(prediction_cv_df)
nrow(prediction_cv_df %>% filter(testing.is_recurrent_parent==prediction_binary))/nrow(prediction_cv_df)
head(data.frame(predictions_cv, testing$is_recurrent_parent))






###### 3) Get feature rankings for 3' versus 5'
# set up feature space
feature_space_5_3 <- read.csv("~/Projects/fusion_ML/features/feature_spaces/five_and_three_prime_parents_not_both.csv", header=TRUE)
feature_space_5_3 <- feature_space_5_3 %>% select(-ensg)

# process
feature_space_5_3_typed <- feature_space_5_3 

# scale relevant variables
feature_space_5_3_typed_scaled <- cbind(feature_space_5_3_typed %>% select(-one_of(variables_to_scale)),
                                              scale(feature_space_5_3_typed %>% select(one_of(variables_to_scale)), center = TRUE, scale = TRUE))

# train test split
inTrain <- createDataPartition(y=feature_space_5_3_typed_scaled$parent_category, p=0.70, list=FALSE)
training <- feature_space_5_3_typed_scaled[inTrain,]
testing <- feature_space_5_3_typed_scaled[-inTrain,]

# penalized regression
pen_5_3 <- penalized(response = parent_category, 
                           penalized = select(training, -parent_category),
                           data=training, model="logistic",
                           lambda1 = 5, steps = 55, standardize = TRUE)

# plot path of coefficients
plotpath(pen_5_3, standardize = TRUE, labelsize = 0.5)
plotpath(pen_5_3)

# show all models
show(pen_5_3)
model_with_10 <- pen_5_3[[47]]; model_with_10

# Examine coefficient values
coeffs_5_3 <- as.data.frame(coefficients(model_with_10)); coeffs_5_3
feature_rankings_5_3 <- as.data.frame(cbind(rownames(coeffs_5_3), coeffs_5_3[,1])) %>% 
  rename(feature=V1, coefficient=V2) %>%
  arrange(-as.numeric(coefficient)); feature_rankings_5_3

write.csv(feature_rankings_5_3, "~/Projects/fusion_ML/modelling/logit_feature_ranking_5_3_prime.csv", 
          row.names = FALSE, quote = FALSE)

# predict
# best parameter via cv1
predictions <- predict(model_with_10, penalized=select(testing, -parent_category), data = testing); predictions
prediction_df <- (data.frame(predictions, testing$parent_category)) %>% mutate(prediction_binary = ifelse(predictions>=0.5, 1, 0)); head(prediction_df)
nrow(prediction_df %>% filter(testing.parent_category==prediction_binary))/nrow(prediction_df)

# CV to optimized lambda, and then prediction
cv_opt <- optL1(response = parent_category, 
                penalized = select(training, -parent_category),
                data=training, model = "logistic", fold=10, minlambda1 = 25, maxlambda1 = 100)

# predict with best param
pen_cv <- penalized(response = parent_category, 
                    penalized = select(training, -parent_category),
                    data=training, model = "logistic",
                    lambda1 = cv_opt$lambda)

predictions_cv <- predict(pen_cv, penalized=select(testing, -parent_category), data = testing); predictions_cv
prediction_cv_df <- (data.frame(predictions_cv, testing$parent_category)) %>% mutate(prediction_binary = ifelse(predictions_cv>=0.5, 1, 0)); head(prediction_cv_df)
nrow(prediction_cv_df %>% filter(testing.parent_category==prediction_binary))/nrow(prediction_cv_df)
# head(data.frame(predictions_cv, testing$parent_category))











###### 4) Get feature rankings for metastasis
# set up feature space
feature_space_metastasis <- read.csv("~/Projects/fusion_ML/features/feature_spaces/fusion_events_downsampled_majority_for_metastasis.csv", header=TRUE)
head(feature_space_metastasis); factor_variables[1] <- "is_metastatic"

# new variables to scale
variables_to_scale_metastasis <- c(paste(variables_to_scale, "_3_prime", sep=""),
                                   paste(variables_to_scale, "_5_prime", sep=""))
# scale relevant variables
feature_space_5_3_scaled <- cbind(feature_space_metastasis %>% select(-one_of(variables_to_scale_metastasis)),
                                  scale(feature_space_metastasis %>% select(one_of(variables_to_scale_metastasis)), center = TRUE, scale = TRUE))

feature_space_5_3_scaled <- feature_space_5_3_scaled %>% select(-c(ensg_3, ensg_5, fusion_id, cancer_type))
head(feature_space_5_3_scaled)

# train test split
inTrain <- createDataPartition(y=feature_space_5_3_scaled$is_metastatic, p=0.70, list=FALSE)
training <- feature_space_5_3_scaled[inTrain,]
testing <- feature_space_5_3_scaled[-inTrain,]

head(training)
# regression
pen_metastasis <- penalized(response = is_metastatic, 
                           penalized = select(training, -is_metastatic),
                           data=training, model="logistic",
                           lambda1 = 10, steps = 70, standardize = TRUE)

plotpath(pen_metastasis, standardize = TRUE, labelsize = 0.5)
plotpath(pen_metastasis)

# show all models
show(pen_metastasis)
model_with_10 <- pen_metastasis[[45]]; model_with_10

# Examine coefficient values
coeffs_metastasis <- as.data.frame(coefficients(model_with_10)); coeffs_metastasis
feature_rankings_metastasis <- as.data.frame(cbind(rownames(coeffs_metastasis), coeffs_metastasis[,1])) %>% 
  rename(feature=V1, coefficient=V2) %>%
  arrange(-as.numeric(coefficient)); feature_rankings_metastasis

write.csv(feature_rankings_metastasis, "~/Projects/fusion_ML/modelling/logit_feature_ranking_metastasis_imputed_downsampling_balanced.csv", 
          row.names = FALSE, quote = FALSE)

# predict
predictions <- predict(model_with_10, penalized=select(testing, -is_metastatic), data = testing); predictions
prediction_df <- (data.frame(predictions, testing$is_metastatic)) %>% mutate(prediction_binary = ifelse(predictions>=0.5, 1, 0)); head(prediction_df)
nrow(prediction_df %>% filter(testing.is_metastatic==prediction_binary))/nrow(prediction_df)

# CV to optimized lambda, and then prediction
cv_opt <- optL1(response = is_metastatic, 
                penalized = select(training, -is_metastatic),
                data=training, model = "logistic", fold=10, minlambda1 = 200, maxlambda1 = 500)

# predict with best param
pen_cv <- penalized(response = is_metastatic, 
                    penalized = select(training, -is_metastatic),
                    data=training, model = "logistic",
                    lambda1 = cv_opt$lambda)

predictions_cv <- predict(pen_cv, penalized=select(testing, -is_metastatic), data = testing); predictions_cv
prediction_cv_df <- (data.frame(predictions_cv, testing$is_metastatic)) %>% mutate(prediction_binary = ifelse(predictions_cv>=0.5, 1, 0)); head(prediction_cv_df)
nrow(prediction_cv_df %>% filter(testing.is_metastatic==prediction_binary))/nrow(prediction_cv_df)



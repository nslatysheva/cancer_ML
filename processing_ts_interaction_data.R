# Processing iid data
require(dplyr)
require(tidyr)
require(igraph)

iid <- read.csv("~/Projects/fusion_ML/data/network_data/iid.human.2016-03.txt", sep="\t")
head(iid); nrow(iid)

# get just uniprot and tissue presence
iid <- iid %>% select(uniprot1, uniprot2, adipose.tissue:uterus) 
# replace empty cells with explicit NA
iid[iid==""] <- NA

nrow(iid)
test <- iid

head(test)

# define network metrics
metrics <- c(
  "n_vertices", "n_edges", 
  "diameter", "avg_path_length", 
  "num_articulation_points", "cohesion", "clique_number", 
  "edge_connect", "vertex_connect", 
  "degree_median", "degree_mean", "degree_sd", 
  "betweenness_median", 
  "closeness_median",
  "network_computation_time")

# set up graph object

convert_to_clean_graph <- function(data, feature_number){
  graph <- 
    data %>% select(uniprot1, uniprot2, feature_number) %>% na.omit() %>%
    select(uniprot1, uniprot2) %>% graph.data.frame(directed = FALSE)
  return(graph)
}

num_tissues <- ncol(test[3:ncol(test)])



ts_network_report <- matrix(nrow = 30, ncol = 15)

for (i in 1:num_tissues){
  
  t_start <- proc.time()
  
  graph = convert_to_clean_graph(data = test, feature_number = i+2)
  
  # number of vertices, edges
  n_vertices <- length(V(graph)); cat("Finished ", metrics[1], "\n")
  n_edges <- length(E(graph)); cat("Finished ", metrics[2], "\n")
  
  # get diameter
  diameter <- diameter(graph); cat("Finished ", metrics[3], "\n")
  avg_path_length <- average.path.length(graph); cat("Finished ", metrics[4], "\n")
  # articulation points
  num_articulation_points <- length(articulation.points(graph)); cat("Finished ", metrics[5], "\n")
 
  # cohesion, cliques
  cohesion <- cohesion(graph); cat("Finished ", metrics[6], "\n")
  clique_number <- clique.number(graph); cat("Finished ", metrics[7], "\n")
  # group adhesion (edge connectivity)
  edge_connect <- edge.connectivity(graph); cat("Finished ", metrics[8], "\n")
  # group cohesions (vertex connectivity)
  vertex_connect <-  vertex.connectivity(graph); cat("Finished ", metrics[9], "\n")
  
  #### average centrality measures
  # degree centrality 
  degree_median <- median(degree(graph)); cat("Finished ", metrics[10], "\n")
  degree_mean <- mean(degree(graph)); cat("Finished ", metrics[11], "\n")
  degree_sd <- sd(degree(graph)); cat("Finished ", metrics[12], "\n")
  # betweenness centrality 
  betweenness_median <- median(betweenness(graph)); cat("Finished ", metrics[13], "\n")
  #betweenness_mean <- mean(betweenness(graph))
  #betweenness_sd <- sd(betweenness(graph))
  # closeness centrality
  closeness_median <- median(closeness(graph)); cat("Finished ", metrics[14], "\n")
  #closeness_mean <- mean(closeness(graph))
  #closeness_sd <- sd(closeness(graph))
  
  # build network summary report
  
  ts_network_report[i,1] <- n_vertices
  ts_network_report[i,2] <- n_edges
  ts_network_report[i,3] <- diameter
  ts_network_report[i,4] <- avg_path_length
  ts_network_report[i,5] <- num_articulation_points
  ts_network_report[i,6] <- cohesion
  ts_network_report[i,7] <- clique_number
  ts_network_report[i,8] <- edge_connect
  ts_network_report[i,9] <- vertex_connect
  ts_network_report[i,10] <- degree_median
  ts_network_report[i,11] <- degree_mean
  ts_network_report[i,12] <- degree_sd
  ts_network_report[i,13] <- betweenness_median
  ts_network_report[i,14] <- closeness_median
  ts_network_report[i,15] <- (proc.time() - t_start)[3]
  
  cat("Finished network number ", i)
  
}

# output final result
ts_network_report_df <- as.data.frame(ts_network_report, row.names = colnames(test)[3:32])
colnames(ts_network_report_df) <- c(
  "n_vertices", "n_edges", 
  "diameter", "avg_path_length", 
  "num_articulation_points", "cohesion", "clique_number", 
  "edge_connect", "vertex_connect", 
  "degree_median", "degree_mean", "degree_sd", 
  "betweenness_median", 
  "closeness_median",
  "network_computation_time"
)

ts_network_report_df

write.table(ts_network_report_df, "~/Projects/fusion_ML/data/network_data/ts_ts_network_report_df.csv", sep=",", quote = FALSE, row.names = TRUE)

# trying to modularize everything
initialize_network_report_matrix <- function(number_of_tissues, number_of_metrics){
  ts_network_report <- matrix(nrow = number_of_tissues, ncol = number_of_metrics)
} 



build_network_report_matrix <- function(network_results, number_of_metrics) {
  
}

output_network_report <- function(network_report_matrix, row_names, col_names){
  ts_network_report_df <- as.data.frame(network_report_matrix, row.names = row_names)
  colnames(ts_network_report_df) <- col_names
  return(ts_network_report_df)
}


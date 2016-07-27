# Summarizing klijn data and venn diagram of fusion classes

require(dplyr)
require(tidyr)
require(ggplot2)

# Plot number of cell lines per tissue supergroup
klijn <- read.csv("~/Projects/fusion_ML/data/fusion_identity_data/klijn_cell_lines_simple.csv", header=TRUE); head(klijn)

summary_klijn <- klijn %>% group_by(tissue_supergroup) %>% 
  summarise(tissue_count = n()) %>% 
  arrange(-tissue_count) %>%
  mutate(y_position = 1) %>%
  mutate(new_tissue_label = paste(tissue_count, "_", tissue_supergroup))

theme_set(theme_bw(base_size = 20))

# Let's just do bar plots
p <- ggplot(summary_klijn, 
            aes(x=reorder(tissue_supergroup, -tissue_count), y=tissue_count, fill=reorder(new_tissue_label, -tissue_count))) + 
  geom_bar(width = 0.5, stat = "identity", show.legend = FALSE) + xlab("Tissue Supergroup") + ylab("Number of Cell Lines") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p 

ggsave(plot = p, filename = "tissue_type_counts.pdf", path = "~/Projects/fusion_ML/figures/", width = 6, height = 4)



############### Summarize tissue migrations with circos plot 
# http://mkweb.bcgsc.ca/tableviewer/
klijn %>% group_by(primary_tissue) %>% summarise(n = n())
klijn %>% group_by(metastatic_tissue) %>% summarise(n = n())

combo <- klijn %>% 
  group_by(primary_tissue, metastatic_tissue) %>%
  summarise(n = n()) %>%
  na.omit() %>% filter(metastatic_tissue != "NULL"); combo

circos_data <- combo %>% spread(metastatic_tissue, n, fill=0)
write.table(circos_data, "~/Projects/fusion_ML/figures/circos_data.txt", sep="\t", quote=FALSE, row.names=FALSE)



############ Summarize number of gene fusions in primary and metastatic 
fusions <- read.table("~/Projects/fusion_ML/data/fusion_identity_data/gene_fusions_joined_cell_line_data.txt", header=TRUE, sep="\t", na.strings = "<NA>"); head(fusions)
head(fusions)

fusion_summary <- fusions %>% mutate(is_metastatic = NA) %>%
  mutate(is_metastatic = ifelse(metastatic_tissue=="NA", "primary", "metastatic")) %>% 
  group_by(is_metastatic) %>% summarise(n_fusions=n()) 

fusion_summary$is_metastatic <- as.factor(fusion_summary$is_metastatic)

p <- ggplot(fusion_summary, aes(x=reorder(is_metastatic, -n_fusions), y=n_fusions, fill=is_metastatic)) + 
  geom_bar(width = 0.5, stat = "identity", show.legend = FALSE, color="grey") + xlab("Metastatic") + ylab("Gene fusion count") +
  coord_cartesian(ylim=c(0,1600)) + scale_fill_manual(values=c("#FF6600", "#C2F23D")); p 

ggsave(plot = p, filename = "fusion_counts.pdf", path = "~/Projects/fusion_ML/figures/", width = 6, height = 6)


############ Check gene set counts + contingency table
feature_space <- read.csv("~/Projects/fusion_ML/features/feature_spaces/feature_space.csv", header=TRUE)

# how many genes of each category?
# 3067 parent genes
feature_space %>% filter(is_parent == 1) %>% distinct(ensg) %>% nrow()
feature_space %>% filter(is_parent == 1) %>% distinct() %>% nrow()
feature_space %>% filter(is_parent == 1) %>% nrow()
# 868 recurrent parents
feature_space %>% filter(is_parent == 1 & is_recurrent_parent == 1) %>% nrow()
feature_space %>% filter(is_parent == 1 & is_recurrent_parent == 0) %>% nrow()
# recurrent and not metastatic
feature_space %>% filter(is_recurrent_parent == 1 & is_metastatic_parent==0) %>% nrow()
# 392 metastatic parents
feature_space %>% filter(is_parent == 1 & is_metastatic_parent == 1) %>% nrow()
feature_space %>% filter(is_metastatic_parent == 1) %>% nrow()
# 177 recurrent and metastatic parents
feature_space %>% filter(is_recurrent_parent == 1 & is_metastatic_parent == 1) %>% nrow()
# 215 non recurrent and metastatic
feature_space %>% filter(is_recurrent_parent == 0 & is_metastatic_parent == 1) %>% nrow()
# other parents
feature_space %>% filter(is_recurrent_parent == 0 & is_metastatic_parent == 0 & is_parent==1) %>% nrow()

# total number of genes 
feature_space %>% nrow()

# looks like recurrent genes are more likely to be metastatic
m <- matrix(c(feature_space %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 1) %>% nrow(), 
              feature_space %>% filter(is_recurrent_parent == 1 & is_metastatic_parent== 0) %>% nrow(), 
              feature_space %>% filter(is_recurrent_parent == 0 & is_metastatic_parent== 1) %>% nrow(), 
              feature_space %>% filter(is_parent == 1, is_recurrent_parent == 0 & is_metastatic_parent== 0) %>% nrow()), nrow=2); m

chisq.test(m)
fisher.test(m)





# Failed attempt at something like a dot plot
# str(summary_klijn)
# 
# p <- ggplot(summary_klijn, aes(new_tissue_label, y_position), fill=tissue_supergroup) + 
#   geom_point(aes(size = tissue_count, color=tissue_supergroup)); p
# p
# bad_grid <- p + theme(axis.text.y=element_blank(), 
#                       axis.title.y=element_blank(),legend.position="none",
#                       panel.background=element_blank(),
#                       panel.border=element_blank(),
#                       panel.grid.major=element_blank(),
#                       panel.grid.minor=element_blank(),
#                       plot.background=element_blank(),
#                       axis.text.x = element_text(angle = 90, hjust = 1),
#                       axis.ticks.y=element_blank())
# 
# bad_grid
# ggsave(plot = bad_grid, filename = "tissue_type_counts.pdf", path = "~/Projects/fusion_ML/figures/", width = 6, height = 4)

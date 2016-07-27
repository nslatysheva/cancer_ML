# Molecular feature figure - PTMcode residues in parents

require(RMySQL)
require(dplyr)
require(reshape2)
require(ggplot2)

source("~/Projects/fusion_ML/scripts/database_credentials.R")

# configure queries
query <- function(query) { a = dbGetQuery(con, statement=query) }

# get number by ensp, group by ensg
# note that zero densities are necessarily excluded here for speed
ptm_density_by_gene <- query("
                             SELECT 
                             all_lengths.*, max_lengths.max_length
                             FROM
                             (SELECT 
                             ensembl_gene_id,
                             ensembl_protein_id,
                             ensp_length,
                             COUNT(DISTINCT accession) AS count_ptm,
                             COUNT(DISTINCT accession) / ensp_length AS ptm_density
                             FROM
                             (SELECT 
                             *
                             FROM
                             natashal.ensg_ensp_uniprot_mapping_97p_genes_89p_prot base
                             JOIN PTMcode_inter_protein ptms ON base.gene = ptm.protein1) x
                             GROUP BY ensembl_gene_id , ensembl_protein_id) all_lengths
                             INNER JOIN
                             (SELECT 
                             ensembl_gene_id, MAX(ensp_length) AS max_length
                             FROM
                             natashal.ensg_ensp_uniprot_mapping_97p_genes_89p_prot base
                             GROUP BY ensembl_gene_id) max_lengths ON all_lengths.ensembl_gene_id = max_lengths.ensembl_gene_id
                             AND all_lengths.ensp_length = max_lengths.max_length
                             GROUP BY ensembl_gene_id
                             ")



# let's get the raw text file
ptmcode <- read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/PTMcode2_associations_between_proteins_HUMAN_final-1.txt", header=TRUE, sep='\t')

# get number of pfam domains by ensp, group by ensg
ensp <- query(" SELECT * FROM natashal.ensp_ensg_gene_length")
ensp <- rename(ensp, c("associated_gene_name"="gene"))

# offline mode
ensp <-  read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/remote_ensp.txt", header=TRUE, sep='\t')
ensg <-  read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/remote_ensg.txt", header=TRUE, sep='\t')
ensg_parent <- read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/remote_ensg_parent.txt", header=TRUE, sep='\t')
ensg_OG <- read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/remote_ensg_OGs.txt", header=TRUE, sep='\t')
ensg_TSG <- read.table("~/Projects/chimera_project/analyses/July_molecular_figure_updates/remote_ensg_TSGs.txt", header=TRUE, sep='\t')


# rename and merge
# limit settable in case not working with full 10 million row table
limit <- nrow(ptmcode)
x <- head(ptmcode, limit)
x <- x[,c("protein1", "ptm_location1", "ptm_type1")]
x <- rename(x, c("protein1"="gene"))
head(x)

# merges and summaries by ensg
ensp <- rename(ensp, c("associated_gene_name" = "gene"))
y <- merge(ensp, x, by="gene"); head(y)

# doing length things properly
z <- ddply(y, c("ensembl_gene_id", "ensembl_protein_id"), function(x) c(num_ptms = length(unique(x$ptm_location1))))

max_lengths <- ddply(ensp, .(ensembl_gene_id), summarise, max_length=max(length))
max_ensp_only <- merge(max_lengths, ensp, by="ensembl_gene_id", all.y = FALSE)
ensg_with_max_length <- unique(max_ensp_only[,c("ensembl_gene_id", "max_length")])

ptm_count_max_lengths <- unique(merge(z[,c("ensembl_gene_id", "num_ptms")], ensg_with_max_length, by="ensembl_gene_id"))

ptm_density_by_gene <- ptm_count_max_lengths
ptm_density_by_gene$density <- ptm_density_by_gene$num_ptms/ptm_density_by_gene$max_length; head(ptm_density_by_gene)

# fill in genes with 0 PTMs
summary(ptm_density_by_gene)

# merge, replace NA
ptm_density_by_gene <- merge(ensg, ptm_density_by_gene, by='ensembl_gene_id', all.x=TRUE); head(ptm_density_by_gene)
ptm_density_by_gene[is.na(ptm_density_by_gene$density),]$density <- 0
summary(ptm_density_by_gene)

# get gene sets
ensg_parent <- query("
                     SELECT DISTINCT
                     (ensembl_gene_id)
                     FROM
                     natashal.chitars_with_two_mapped_proteins c
                     JOIN
                     ensp_ensg_gene_length p ON c.ensp = p.ensembl_protein_id
                     ")

ensg_OG <- query("SELECT DISTINCT(ensg) FROM natashal.oncogenes")
ensg_TSG <- query("SELECT DISTINCT(ensembl_gene_id) FROM natashal.tumour_suppresor_genes;")

# overlap gene sets
ptm_density_by_gene$category <- NA
# label parents
ptm_density_by_gene[which(ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id),]$category <- 'parent'
# label non-parents
ptm_density_by_gene[which(!ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id),]$category <- 'non-parent'
# label parent OGs
ptm_density_by_gene[which((ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_OG$ensg)),]$category <- "OG_parent"
# label parent TSGs
ptm_density_by_gene[which((ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_TSG$ensembl_gene_id)),]$category <- "TSG_parent"

ptm_density_by_gene$category <- as.factor(ptm_density_by_gene$category)
summary(ptm_density_by_gene$category)
ddply(ptm_density_by_gene, "category", function(x) summary(x$density))

# values for plotting
# medians <- ddply(ptm_density_by_gene, "category", function(x) c(median=median(x$ANCHOR_ptm_density)))
third_quartile <- ddply(ptm_density_by_gene, "category", function(x) c(third_quartile=quantile(x$density)[4]))
colnames(third_quartile) <- c("category", "third_quartile")

# plot
theme_set(theme_gray(base_size = 32))

p <- ggplot(ptm_density_by_gene, aes(x = reorder(category, -density), y = density)) + 
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = category)) +
  scale_fill_manual(values=c('#f0f0f0', '#67000d', '#cb181d', '#a50f15')) +
  xlab("") + ylab("") + guides(fill=FALSE) +
  coord_cartesian(ylim=c(-0.004,0.125)) +
  annotate("text",x=1,y=0.8*third_quartile$third_quartile[third_quartile$category=="OG_parent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="OG_parent"))) + 
  annotate("text",x=2,y=0.8*third_quartile$third_quartile[third_quartile$category=="TSG_parent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="TSG_parent"))) + 
  annotate("text",x=3,y=0.7*third_quartile$third_quartile[third_quartile$category=="parent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="parent"))) + 
#   annotate("text",x=4,y=0.8*third_quartile$third_quartile[third_quartile$category=="non-parent"], size=8, colour = "black", label=nrow(subset(ptm_density_by_gene, category=="non-parent"))) 
  annotate("text",x=4,y=0.006, size=8, colour = "black", label=nrow(subset(ptm_density_by_gene, category=="non-parent"))) 


p 

ggsave(p, height = 8, width = 6, file = "Figure 4H_interacting February.pdf", path="~/Projects/chimera_project/analyses/July_molecular_figure_updates/")

pairwise.wilcox.test(ptm_density_by_gene$density, ptm_density_by_gene$category)

# summary stats
# statistics
# the statistics
ptm_summary <- ddply(ptm_density_by_gene, "category", function(x) summary(x$density))

# fold increases of non-parents
ptm_summary[which(ptm_summary$category=="OG_parent"),]$Mean/ptm_summary[which(ptm_summary$category=="non-parent"),]$Mean
ptm_summary[which(ptm_summary$category=="TSG_parent"),]$Mean/ptm_summary[which(ptm_summary$category=="non-parent"),]$Mean
ptm_summary[which(ptm_summary$category=="parent"),]$Mean/ptm_summary[which(ptm_summary$category=="non-parent"),]$Mean


##### coordinating ptms in fusion and non-fusion OGs and TSGs
# overlap gene sets
ptm_density_by_gene$category <- NA
# label parent OGs
ptm_density_by_gene[which((ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_OG$ensg)),]$category <- "1_OG_parent"
# label parent TSGs
ptm_density_by_gene[which((ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_TSG$ensembl_gene_id)),]$category <- "3_TSG_parent"
# label nonparent OGs
ptm_density_by_gene[which((!ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_OG$ensg)),]$category <- "2_OG_nonparent"
# label nonparent TSGs
ptm_density_by_gene[which((!ptm_density_by_gene$ensembl_gene_id %in% ensg_parent$ensembl_gene_id) & (ptm_density_by_gene$ensembl_gene_id %in% ensg_TSG$ensembl_gene_id)),]$category <- "4_TSG_nonparent"

# clean up
ptm_density_by_gene$category <- as.factor(ptm_density_by_gene$category)
ptm_density_by_gene <- ptm_density_by_gene[,c("ensembl_gene_id", "density", "category")]
ptm_density_by_gene<-na.omit(ptm_density_by_gene)
summary(ptm_density_by_gene$category)
ddply(ptm_density_by_gene, "category", function(x) summary(x$density))

# values for plotting
third_quartile <- ddply(ptm_density_by_gene, "category", function(x) c(third_quartile=quantile(x$density)[4]))
colnames(third_quartile) <- c("category", "third_quartile")

# plot
theme_set(theme_gray(base_size = 32))

p <- ggplot(ptm_density_by_gene, aes(x = category, y = density)) + 
  geom_boxplot(outlier.shape=NA, notch = TRUE, aes(fill = category)) +
  scale_fill_manual(values=c('#084081', '#0868ac', '#004529', '#238443')) +
  xlab("") + ylab("") + guides(fill=FALSE) +
  coord_cartesian(ylim=c(-0.004,0.125)) +
  annotate("text",x=1,y=0.8*third_quartile$third_quartile[third_quartile$category=="1_OG_parent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="1_OG_parent"))) + 
  annotate("text",x=2,y=0.8*third_quartile$third_quartile[third_quartile$category=="2_OG_nonparent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="2_OG_nonparent"))) + 
  annotate("text",x=3,y=0.8*third_quartile$third_quartile[third_quartile$category=="3_TSG_parent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="3_TSG_parent"))) + 
  annotate("text",x=4,y=0.8*third_quartile$third_quartile[third_quartile$category=="4_TSG_nonparent"], size=8, colour = "white", label=nrow(subset(ptm_density_by_gene, category=="4_TSG_nonparent"))) 


p 

wilcox.test(subset(ptm_density_by_gene, category=="1_OG_parent")$density, subset(ptm_density_by_gene, category=="2_OG_nonparent")$density)
wilcox.test(subset(ptm_density_by_gene, category=="3_TSG_parent")$density, subset(ptm_density_by_gene, category=="4_TSG_nonparent")$density)


ggsave(p, height = 8, width = 6, file = "coord_ptms_onc_tsg February.pdf", path="~/Projects/chimera_project/analyses/July_molecular_figure_updates/")

# stats

ptm_summary <- ddply(ptm_density_by_gene, "category", function(x) summary(x$density))

# fold increases of non-parents
ptm_summary[which(ptm_summary$category=="1_OG_parent"),]$Mean/ptm_summary[which(ptm_summary$category=="2_OG_nonparent"),]$Mean
ptm_summary[which(ptm_summary$category=="3_TSG_parent"),]$Mean/ptm_summary[which(ptm_summary$category=="4_TSG_nonparent"),]$Mean




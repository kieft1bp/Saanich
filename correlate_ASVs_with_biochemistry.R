library(ggplot2)
library(phyloseq)
library(plyr)
library(yarrr)
library(reshape2)
library(tidyverse)
library(dplyr)

load("./dada2run2_track_ps.RData")
ps_working = ps

# Update metadata and tax
new_metadata = read.table("dada2_run2_metadata_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
sample_data(ps_working) = new_metadata
new_tax = read.table("dada2_run2_taxonomy_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
tax_table(ps_working) = as.matrix(new_tax)

# Subset the metadata by only predictive values
flat_meta_table_subset = new_metadata[,c(6,34,35,28:32,10:15,23,25)]
flat_meta_table_subset = sapply( flat_meta_table_subset, as.numeric )
rownames(flat_meta_table_subset) = rownames(new_metadata)
flat_meta_table_subset_complete = flat_meta_table_subset[complete.cases(flat_meta_table_subset),]
calculate_correlation_meta = as.data.frame(flat_meta_table_subset_complete)
calculate_correlation_meta$Sample = rownames(calculate_correlation_meta)

# Subset counts table and get only topmost abundant taxa to plot
physeq_TAXON_matrix = as(otu_table(ps_working), "matrix")
physeq_TAXON_matrix = physeq_TAXON_matrix[rownames(physeq_TAXON_matrix) %in% colnames(calculate_bioenv_meta),]
ASV_means = as.data.frame(as.matrix(colMeans(physeq_TAXON_matrix)))
ASV_means = ASV_means[order(-ASV_means$V1),,drop=F]
top_ASVs = rownames(ASV_means)[1:50]
calculate_correlation_taxa = as.data.frame(physeq_TAXON_matrix)
calculate_correlation_taxa = calculate_correlation_taxa[, colnames(calculate_correlation_taxa) %in% top_ASVs]
calculate_correlation_taxa$Sample = rownames(calculate_correlation_taxa)

# Combine the metadata and top taxa data
calculate_correlation_combined = merge(calculate_correlation_meta, calculate_correlation_taxa, by="Sample")
rownames(calculate_correlation_combined) = calculate_correlation_combined$Sample

# Run correlation tests on the variable matrix
find_correlations = calculate_correlation_combined[,2:ncol(calculate_correlation_combined)]
find_correlations_calc = cor(find_correlations, method = "pearson")

# Clean the correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(find_correlations_calc)
upper_tri = cormat

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap object
plot_correlation_matrix <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  xlab("") +
  ylab("") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 1),
        axis.text.y = element_text(size=6)) +
  coord_fixed(ratio = 1:1)

# Plot
plot_correlation_matrix

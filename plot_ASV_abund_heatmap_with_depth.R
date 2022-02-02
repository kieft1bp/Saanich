library(ggplot2)
library(phyloseq)
library(plyr)
library(vegan)
library(reshape2)
library(tidyverse)
library(dplyr)

setwd("~/Hallam_lab/Saanich")
load("./dada2run2_track_ps.RData")
ps_working = ps
#ps_working = subset_samples(ps_working, Depth != "85m") # This is not a consistent depth across samples
ps_working = subset_samples(ps_working, Depth != "80m") # This is not a consistent depth across samples

family = "Thioglobaceae"  # Clade_I, Clade_II, Clade_III, Clade_IV, Unassigned_Marinimicrobia, SAR116_clade, Thioglobaceae, Arcobacteraceae, Nitrosopumilaceae, Actinomarinaceae, Ectothiorhodospiraceae, Rhodobacteraceae, Unassigned_SAR324 clade, Methylomonadaceae, Flavobacteriaceae, Microtrichaceae, Unassigned_Marine Group II, Scalinduaceae
number_of_ASVs = 50    # This will take the top X ASVs of the family. This will be unique to a family, but it's not realistic to plot all the rare ones. If you get a "subscript out of bounds" error, try lowering this value
order_ASVs_by = "abundance"   # 'abundance' or 'clustering'

# Update metadata and tax
new_metadata = read.table("dada2_run2_metadata_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
sample_data(ps_working) = new_metadata
new_tax = read.table("dada2_run2_taxonomy_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
tax_table(ps_working) = as.matrix(new_tax)

# Make a relative ps object
ps_relative  = transform_sample_counts(ps_working, function(x) x / sum(x) )
flat_taxonomy_table = as.data.frame(as(tax_table(ps_relative), "matrix"))
flat_taxonomy_table = flat_taxonomy_table %>% mutate_if(is.character, str_replace_all, ' ', '_')
#write.table(flat_taxonomy_table, "taxonomy_table.tsv", row.names = T, sep="\t", quote = F)
flat_taxonomy_table$ASV = rownames(flat_taxonomy_table)
flat_taxonomy_table$ASV_num = gsub("ASV","",flat_taxonomy_table$ASV)
flat_taxonomy_table$ASV_num = as.numeric(flat_taxonomy_table$ASV_num)

# Subset ASV table and aggregate by phylum
physeq_TAXON = subset_taxa(ps_relative, Family==family)
#taxa_names(physeq_TAXON) <- paste(tax_table(physeq_TAXON)[,6]," (",rownames(tax_table(physeq_TAXON)),")", sep="")

# Export flat files from phyloseq relative abundance object

flat_counts_table = as.data.frame(as(otu_table(physeq_TAXON), "matrix"))
flat_meta_table = as.data.frame(as(sample_data(physeq_TAXON), "matrix"))
flat_meta_table$sampleID = rownames(flat_meta_table)

# Combine ASV table and metadata
physeq_TAXON_matrix = as.matrix(otu_table(physeq_TAXON))
ASV_means = as.data.frame(as.matrix(colMeans(physeq_TAXON_matrix)))
ASV_means = ASV_means[order(-ASV_means$V1),,drop=F]
top_X_ASVs = rownames(ASV_means)
#top_X_ASVs = top_X_ASVs[1:50]
#top_X_ASVs = top_X_ASVs[26:50]
physeq_TAXON_df = as.data.frame(physeq_TAXON_matrix)
physeq_TAXON_df = physeq_TAXON_df[, colnames(physeq_TAXON_df) %in% top_X_ASVs]
physeq_TAXON_melted = melt(as.matrix(physeq_TAXON_df))
colnames(physeq_TAXON_melted) = c("sampleID","ASV","relabund")
physeq_TAXON_melted$relabund = physeq_TAXON_melted$relabund*100
TAXON_plot_data = join(physeq_TAXON_melted,flat_meta_table,by="sampleID")
is.na(TAXON_plot_data) = "NAN"

# Aggregate by mean for a given parameter
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV2", "ASV6", "ASV10", "ASV11", "ASV16","ASV453"),"ASV1",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV5", "ASV13", "ASV20", "ASV23", "ASV15"),"ASV3",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV80", "ASV147", "ASV167", "ASV192", "ASV193"),"ASV58",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV2265", "ASV2407", "ASV2435", "ASV4261", "ASV4566"),"ASV1239",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV676", "ASV985", "ASV1092", "ASV1545", "ASV1736"),"ASV559",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV2721", "ASV3357", "ASV3810", "ASV5797", "ASV6401"),"ASV2690",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV2080", "ASV4384", "ASV4650", "ASV4797", "ASV4895"),"ASV1738",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV580", "ASV912", "ASV1049", "ASV1059", "ASV1188"),"ASV380",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV245", "ASV397", "ASV442", "ASV495", "ASV604"),"ASV210",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV12520", "ASV12306", "ASV9101"),"ASV7836",as.character(TAXON_plot_data$ASV))
#TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV4980", "ASV5270", "ASV5541"),"ASV3623",as.character(TAXON_plot_data$ASV))
#TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV1898"),"ASV1842",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV1432","ASV1842"),"ASV1285",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV7657"),"ASV7657",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV5053"),"ASV4504",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV2020","ASV1921"),"ASV1334",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV103"),"ASV103",as.character(TAXON_plot_data$ASV))
TAXON_plot_data$ASV = ifelse(TAXON_plot_data$ASV %in% c("ASV6432","ASV5255"),"ASV4631",as.character(TAXON_plot_data$ASV))
top_X_ASVs = c("ASV1","ASV3","ASV58","ASV1239","ASV559","ASV2690","ASV1738","ASV380","ASV210","ASV1285","ASV1334","ASV103")
TAXON_plot_data = TAXON_plot_data[TAXON_plot_data$ASV %in% top_X_ASVs,]

TAXON_plot_data_agg = aggregate(relabund ~ depth + ASV + cruise, TAXON_plot_data, sum)
#TAXON_plot_data_agg = aggregate(relabund ~ depth + ASV, TAXON_plot_data_agg, mean)

# If you use two variables, it's possible to unmelt (dcast) the df and create hierarchical clustering of ASVs
# Make the ASV distance matrix
TAXON_plot_data_agg_dcast = dcast(TAXON_plot_data_agg, ASV ~ depth, value.var = "relabund")
rownames(TAXON_plot_data_agg_dcast) = TAXON_plot_data_agg_dcast$ASV
TAXON_plot_data_agg_dcast$ASV = NULL
TAXON_plot_data_agg_dcast_dist = vegdist(TAXON_plot_data_agg_dcast)
# Run the hierarchical clustering algorithm and get the new order of ASVs
heatmap_cluster = hclust(TAXON_plot_data_agg_dcast_dist)
heatmap_cluster_labs = as.data.frame(heatmap_cluster$labels)
heatmap_cluster_labs$id <- as.integer(row.names(heatmap_cluster_labs)) 
heatmap_cluster_labs = heatmap_cluster_labs[order(match(heatmap_cluster_labs$id, heatmap_cluster$order)), , drop = FALSE]
heatmap_labs_ordered = as.character(heatmap_cluster_labs[,1])

# Factorize some of the variables
# if(order_ASVs_by == "abundance"){
#   TAXON_plot_data_agg$ASV = factor(TAXON_plot_data_agg$ASV, levels = top_X_ASVs)
# }
# if(order_ASVs_by == "clustering"){
#   TAXON_plot_data_agg$ASV = factor(TAXON_plot_data_agg$ASV, levels = heatmap_labs_ordered)
# }
#TAXON_plot_data_agg = TAXON_plot_data_agg[order(TAXON_plot_data_agg$taxon),]
#TAXON_plot_data_agg$ASV = factor(TAXON_plot_data_agg$ASV, levels = unique(TAXON_plot_data_agg$ASV))
TAXON_plot_data_agg$ASV = factor(TAXON_plot_data_agg$ASV, levels = c("ASV1738","ASV2690","ASV58","ASV1334","ASV103","ASV3","ASV1239","ASV559","ASV1285","ASV380","ASV1","ASV210"))
TAXON_plot_data_agg$depth = as.numeric(TAXON_plot_data_agg$depth)
TAXON_plot_data_agg$depth = factor(TAXON_plot_data_agg$depth, levels = c(10,20,40,60,75,85,90,97,100,110,120,135,150,165,185,200))

#TAXON_plot_data_agg$relabund = ifelse(TAXON_plot_data_agg$relabund == 0, NA, as.numeric(TAXON_plot_data_agg$relabund))

# Plot
ggplot(TAXON_plot_data_agg, aes(x=ASV, y=depth, fill = sqrt(relabund))) +
  geom_tile(color = "black") +
  facet_wrap(~ASV, scales = "free_x", nrow=1) +
  scale_fill_gradient(low = "gray90", high = "black", na.value="white", name="Sqrt. % Abund.") +
  xlab("") +
  ylab("Depth (m)") +
  theme_minimal() +
  ggtitle(paste(family," ASV abundances",sep="")) +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(color="black",angle = 60, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size=10,color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank())


ggplot(data = TAXON_plot_data_agg, aes(x=depth, y=relabund, group = ASV)) +
  geom_smooth(se = T, method = 'loess', span=0.3, size=0.5, color="black") +
  geom_point(size=0.5, color="black") +
  #geom_boxplot(aes(fill=depth)) +
  facet_grid(~ASV, scales = "free_x") +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  xlab("Depth (m)") +
  ylab("Relative Community Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color="black",angle = 60, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(color="black",size=10),
        panel.background = element_rect(fill=NA,color="black"))
  

#ggsave(paste(family,"_ASVs_vertical_distribution.svg",sep = ""), ASV_lineplot_with_O2, device = "svg",scale = 1, width = 8, height = 4, units = "in", dpi = 600)



flat_meta_table = as.data.frame(as(sample_data(ps_relative), "matrix"))
flat_meta_table_subset = flat_meta_table[,c(6,34,35,28:32,10:15,23,25)]
flat_meta_table_subset = sapply( flat_meta_table_subset, as.numeric )
rownames(flat_meta_table_subset) = rownames(new_metadata)
flat_meta_table_subset_complete = flat_meta_table_subset[complete.cases(flat_meta_table_subset),]
calculate_correlation_meta = as.data.frame(flat_meta_table_subset_complete)
calculate_correlation_meta$Beam_Transmission_Cstar = NULL
calculate_correlation_meta$PAR_Irradiance_Licor = NULL

depth_averages_for_meta_correlation = aggregate(. ~ depth, calculate_correlation_meta, mean)
depth_averages_for_meta_correlation = depth_averages_for_meta_correlation[depth_averages_for_meta_correlation$depth != 80,]

depth_averages_for_taxa_correlation = aggregate(relabund ~ depth + ASV, TAXON_plot_data_agg, mean)
depth_averages_for_taxa_correlation = dcast(depth_averages_for_taxa_correlation, depth ~ ASV, value.var = "relabund")

taxa_plus_meta_correlation = merge(depth_averages_for_taxa_correlation, depth_averages_for_meta_correlation, by = "depth")
rownames(taxa_plus_meta_correlation) = taxa_plus_meta_correlation$depth

# Run correlation tests on the variable matrix
find_correlations = taxa_plus_meta_correlation[,2:ncol(taxa_plus_meta_correlation)]
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
newnames = colnames(upper_tri)
newnames = gsub("mean_ch4","Methane",newnames)
newnames = gsub("si$","Silicon",newnames)
newnames = gsub("h2s","Hydrogen Sulfide",newnames)
newnames = gsub("po4","Phosphate",newnames)
newnames = gsub("nh4","Ammonium",newnames)
newnames = gsub("_Wetstar","",newnames)
newnames = gsub("no2","Nitrite",newnames)
newnames = gsub("no3","Nitrate",newnames)
newnames = gsub("mean_n2o","Nitrous Oxide",newnames)
newnames = gsub("_uM","",newnames)
colnames(upper_tri) = newnames
rownames(upper_tri) = newnames

melted_cormat <- melt(upper_tri, na.rm = TRUE)


# Create a ggheatmap object
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Rho\n") +
  xlab("") +
  ylab("") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(color = "black", angle = 90, size = 10, hjust = 1, vjust =0.5),
        axis.text.y = element_text(size=10, color="black"),
        legend.text = element_text(color="black", size=10)) +
  coord_fixed(ratio = 1:1)




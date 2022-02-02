library(vegan)
library(phyloseq)
library(plyr)
library(caret)
library(reshape2)
library(tidyverse)
library(dplyr)
library(kableExtra)

load("./dada2run2_track_ps.RData")
ps_working = ps

# Update metadata and tax
new_metadata = read.table("dada2_run2_metadata_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
is.na(new_metadata) = "NAN"
sample_data(ps_working) = new_metadata
new_tax = read.table("dada2_run2_taxonomy_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
tax_table(ps_working) = as.matrix(new_tax)

# Make a relative ps object
ps_relative  = transform_sample_counts(ps_working, function(x) x / sum(x) )

# Export flat files from phyloseq relative abundance object
flat_counts_table = as.data.frame(as(otu_table(ps_relative), "matrix"))
flat_meta_table = as.data.frame(as(sample_data(ps_relative), "matrix"))
flat_meta_table$sampleID = rownames(flat_meta_table)

# Subset the metadata just by biochemical variables
flat_meta_table_subset = flat_meta_table[,c(6,34,35,28:32,10:15,23,25)]
flat_meta_table_subset = sapply( flat_meta_table_subset, as.numeric )
rownames(flat_meta_table_subset) = rownames(flat_meta_table)
flat_meta_table_subset_complete = flat_meta_table_subset[complete.cases(flat_meta_table_subset),]
calculate_bioenv_meta = as.data.frame(t(flat_meta_table_subset_complete))

# Extract a list of sample and biochem var names
sample_names = colnames(calculate_bioenv_meta)
metadata_names = rownames(calculate_bioenv_meta)

# Find correlations between variables in order to remove biochem variables that highly correlate with other biochem parameters
find_correlations = cor(t(calculate_bioenv_meta), method = "pearson")
remove_correlated_vars = findCorrelation(find_correlations, cutoff = 0.75, verbose = F,names = T,exact = F)
calculate_bioenv_meta$metadata = rownames(calculate_bioenv_meta)
calculate_bioenv_meta = calculate_bioenv_meta[!(calculate_bioenv_meta$metadata %in% remove_correlated_vars),]
calculate_bioenv_meta$metadata = NULL

# Create an ASV table that is subsetted by only the new sample_names list (which is smaller than the original number of samples)
calculate_bioenv_ASV = flat_counts_table[rownames(flat_counts_table) %in% sample_names,]
calculate_bioenv_ASV = calculate_bioenv_ASV[rowSums(calculate_bioenv_ASV) > 0,]

# Put metadata and ASV dataframe into the correct orientation
calculate_bioenv_meta = as.data.frame(t(calculate_bioenv_meta))
calculate_bioenv_ASV = as.data.frame(calculate_bioenv_ASV)
calculate_bioenv_meta = calculate_bioenv_meta[rownames(calculate_bioenv_meta) %in% rownames(calculate_bioenv_ASV),]

# Run the bioenv model
bioenv_correlation_model = bioenv(comm = wisconsin(calculate_bioenv_ASV), env = calculate_bioenv_meta)

bioenv_correlation_model_results = cbind(data.frame(c(unlist(summary(bioenv_correlation_model)[3]))), data.frame(c(as.numeric(unlist(summary(bioenv_correlation_model)[2])))))
colnames(bioenv_correlation_model_results) = c("vars","corr")
bioenv_correlation_model_results = bioenv_correlation_model_results[order(-bioenv_correlation_model_results$corr),]
bioenv_correlation_model_results$vars = gsub(" "," + ",bioenv_correlation_model_results$vars)
rownames(bioenv_correlation_model_results) = NULL
colnames(bioenv_correlation_model_results) = c("Variable Combinations","Explained Variance")

table_2_print = kable(bioenv_correlation_model_results, "html", caption = "Table of BIOENV models") %>% kable_styling("striped") %>% scroll_box(width = "100%")
table_2_print

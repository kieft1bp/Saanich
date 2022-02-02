library(ggplot2)
library(phyloseq)
library(plyr)
library(vegan)
library(reshape2)
library(tidyverse)
library(dplyr)
library(psadd)

load("./dada2run2_track_ps.RData")
ps_working = ps
ps_working = subset_samples(ps_working, Cruise == "SI067")    # Best just to check out one cruise at first
ps_working = filter_taxa(ps_working, function(x) var(x) > 0, TRUE)

# Update metadata and tax
new_metadata = read.table("dada2_run2_metadata_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
sample_data(ps_working) = new_metadata
new_tax = read.table("dada2_run2_taxonomy_clean.tsv", header=T, row.names = 1, sep="\t", stringsAsFactors = F)
tax_table(ps_working) = as.matrix(new_tax)

# Plot Krona plot into an HTML that should appear in the browser
plot_krona(ps_working, "krona_plot.html", "Depth", trim = F)




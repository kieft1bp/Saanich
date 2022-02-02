library(ggplot2)
library(phyloseq)

load("./dada2run2_track_ps.RData")
ps_working = ps
ps_working = subset_samples(ps_working, Cruise == "SI022")
ps_working = filter_taxa(ps_working, function(x) var(x) > 0, TRUE)

ps_nmds = ordinate(ps, "NMDS", "bray")
p = plot_ordination(ps, ps_nmds, type = "samples", color = "Depth") + geom_point()
p

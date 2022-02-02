library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)

# Read in MAG to ASV mapping summary table
MAG_to_ASV_summary = read.table("MAG_to_ASV_summary.tsv", header=T, quote="")
MAG_to_ASV_summary$MAG = gsub("_k.*$","",MAG_to_ASV_summary$MAG)

# Create a code for determining RNA content
MAG_to_ASV_summary$RNA = ifelse((MAG_to_ASV_summary$X5S == 0 & MAG_to_ASV_summary$X23S == 0), 124, 4)
MAG_to_ASV_summary$RNA = ifelse((MAG_to_ASV_summary$X23S != 0 & MAG_to_ASV_summary$X5S != 0), 8, as.integer(MAG_to_ASV_summary$RNA) )
MAG_to_ASV_summary$RNA = as.integer(MAG_to_ASV_summary$RNA)
MAG_to_ASV_summary = unique(MAG_to_ASV_summary[,c(1:3,7:11)])

# Filter by quality
MAG_to_ASV_summary = MAG_to_ASV_summary[MAG_to_ASV_summary$ASV_LEN > 200,]
MAG_to_ASV_summary = MAG_to_ASV_summary[!grepl("bin.41",MAG_to_ASV_summary$MAG),]

# Create ordered table
MAG_to_ASV_summary = MAG_to_ASV_summary[order(-MAG_to_ASV_summary$ASV_PID),]
MAG_to_ASV_summary = MAG_to_ASV_summary[order(MAG_to_ASV_summary$MAG),]
MAG_to_ASV_summary = MAG_to_ASV_summary[match(unique(MAG_to_ASV_summary$MAG), MAG_to_ASV_summary$MAG),]

# Create plotting aesthetics
ASV = as.character(unique(MAG_to_ASV_summary$ASV))
colors = as.character(brewer.pal(length(ASV),"Set1"))
unique_ASV_to_color = cbind(ASV,colors)
MAG_to_ASV_summary = merge(MAG_to_ASV_summary,unique_ASV_to_color, by = "ASV")
MAG_to_ASV_summary$ASV = factor(MAG_to_ASV_summary$ASV, levels = unique(MAG_to_ASV_summary$ASV))
MAG_to_ASV_summary$MAG = factor(MAG_to_ASV_summary$MAG, levels = unique(MAG_to_ASV_summary$MAG))

# Plot Comp by Cont
ggplot(MAG_to_ASV_summary, aes(x=Completeness, y = Contamination, shape = RNA, color = colors, label = MAG)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2, max.overlaps = 17) +
  geom_point(shape = 1, size = 4.2) +
  scale_shape_identity() +
  theme_minimal() +
  scale_color_identity(guide = "legend", labels = c("ASV3 (Arctic96)","ASV210 (SUP05)","ASV1 (SUP05)")) +
  scale_x_continuous(limits = c(70,100)) +
  scale_y_continuous(limits = c(0,10)) +
  xlab("Completeness (%)") +
  ylab("Contamination (%)") +
  theme(axis.text = element_text(color="black", size= 10),
        panel.border = element_rect(fill = NA, color = "black"))

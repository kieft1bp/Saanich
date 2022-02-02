# R Script for dada2 pipeline to convert 16S RNA data to ASVs from Saanich Inlet
# Load packages
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Indicate path to samples on the shamwow2 server
path <- "/mnt/nfs/sharknado/LimsData/OMZs/Saanich/2008_to_2012_16S_time_series/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and 
# SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = 6))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = 6))
# Extract sample names, assuming filenames have format: CRUISE_DEPTH_PLATE_DATE_SAMPLE_RX.fastq, 
# where "CRUISE_DEPTH_PLATE_DATE_SAMPLE" = to the SAMPLENAME in SAMPLENAME_XXX.fastq
#Therefore, we have to apply the last "_" instead of the first (as shown in the tutorial)
sample.fields <- lapply(strsplit(basename(fnFs), "_"), `[`, c(1,2,3,4,5))
sample.names <- sapply(sample.fields, function(fields) paste(fields[1], fields[2], fields[3], fields[4], fields[5], sep="_"))

# Filter & Trimming Steps: 
# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Set standard filtering parameters; maxEE = max number of expected errors allowed 
# Filtering and trimming processes will occur with these parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=6,
                     compress=TRUE, multithread=6, trimLeft = c(15,15)) 

# Note: parameters were adjusted to include trimLeft = c(15,15) to avoid chimeras caused by primer seqs

# Inspect number of samples after filter/trim
head(out)

# Learn the Error Rates: whereby the DADA2 algorithm uses a model (err) to adjust 
# for the different set of error rates for different datasets, using function errF()
errF <- learnErrors(filtFs, multithread=6)
errR <- learnErrors(filtRs, multithread=6)
# Plotting errors in PDF file
pdf("plotErrors.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

# Applying CORE SAMPLE INFERENCE algorithm to filtered/trimmed sequences to ID
# real variants
dadaFs <- dada(filtFs, err=errF, multithread=6)
dadaRs <- dada(filtRs, err=errR, multithread=6)

# Merging Steps:
# Merging paired reads (F and R together) to obtain full denoised sequences 
# Aligning denoised forward reads with reverse-complement of  the corresponding denoised 
# reverse reads to construct ontigs 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Constructing a sequence table (higher resolution than OTU table; ASV table)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Chimera Removal Step: 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=6, verbose=TRUE)

# Inspect dimensions of the chimera-removed df
dim(seqtab.nochim)

# Calculate proportion of non-chimeric merged sequence variants/total 
# merged sequence variants (proportion). This should be close to 1
sum(seqtab.nochim)/sum(seqtab)

#Track Reads through Pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Read counts are tracked and stored in the object "track" at each step
head(track)

# Taxonomoy Assignment steps:
taxa <- assignTaxonomy(seqtab.nochim, "/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/Dada2/silva_database/silva_nr99_v138.1_train_set.fa.gz", multithread=6)
# Assigning species: 
taxa <- addSpecies(taxa, "/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/Dada2/silva_species_assignment_v138.1fa.gz")

# Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# HANDOFF to PHYLOSEQ
samples.out <- rownames(seqtab.nochim)
cruise <- sapply(strsplit(samples.out, "_"), `[`, 1)
depth <- sapply(strsplit(samples.out, "_"), `[`, 2)
plate <- sapply(strsplit(samples.out, "_"), `[`, 3)
date <- sapply(strsplit(samples.out, "_"), `[`, 4)
code <-sapply(strsplit(samples.out, "_"), `[`, 5)

samdf <- data.frame(Cruise=cruise, Depth=depth, Plate=plate, Date=date, Code=code) 
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps = filter_taxa(ps, function(x) var(x) > 0, TRUE)

# Keep absolute count PS object
ps_absolute = ps

# Create relative abundance PS object
ps = transform_sample_counts(ps, function(x) x / sum(x) )

save.image("Saanich_ASV_dada2_object.rdata")

# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# THIS IS THE SECOND STEP AFTER DADA2
## rownames-match.r and replicates-contaminated.r must be run before this!
# ------------------------------------------------------------------


## Sets up the Environment and Libraries

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(tibble)

# ------------------------------------------------------------------
# Loads in Data
# ------------------------------------------------------------------
# Loads in output from DADA2 + filtered seqtab.nochim and samdf from rownames-match.r
load("replicates-contaminated.RData")

# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim_filt, taxa_are_rows=FALSE), 
                   sample_data(samdf_filt), 
                   tax_table(taxam))

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.16s))
names(dna) <- taxa_names(ps.16s)
ps.raw <- merge_phyloseq(ps.16s, dna)
taxa_names(ps.16s) <- paste0("ASV", seq(ntaxa(ps.16s)))

nsamples(ps.16s)

# Filters out any Mammalia and NA
# CHANGE TO FILTER ANYTHING THAT ISN'T ACTINOPTERI  
ps.16s <- subset_taxa(ps.16s, Class!="Mammalia")
ps.16s <- subset_taxa(ps.16s, Kingdom!="Bacteria")

# Remove samples with total abundance == 0
ps.16s <- prune_samples(sample_sums(ps.16s) > 0, ps.16s)

## MERGE TO SPECIES HERE (TAX GLOM)
ps.16s = tax_glom(ps.16s, "Species", NArm = FALSE)

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)
lowcount.minor <- genefilter_sample(ps.16s, f1, A=1)
ps.16s.minor <- prune_taxa(lowcount.minor, ps.16s)

# Filtering to remove taxa with less than 1% of diet in 10 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.16s, f1, A=10)
ps.16s.major <- prune_taxa(lowcount.filt, ps.16s)

# Plots stacked bar plot of abundance - to confirm presence of NA's
# plot_bar(ps.16s, fill="Species")

# Plots stacked bar plot of abundance - to confirm presence of NA's
abs <- plot_bar(ps.16s, fill="Species")
abs

# Saves relative abundance plot
ggsave("Deliverables/ALL/ABSOLUTE.png", plot = abs, width = 30, height = 8, units = "in", dpi = 300)


# Calculates proportional abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.16s.minor.rel <- transform_sample_counts(ps.16s.minor, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.16s.major.rel <- transform_sample_counts(ps.16s.major, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.16s.minor.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.16s.major.rel))), arr.ind = TRUE)

#SAVES

save(ps.16s, ps16s.rel, samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "srkw-ps.16s.RData")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps16s.rel, ps.16s, file = "SRKW-diet-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.minor, ps.16s.minor.rel, file = "SRKW-diet-FILT-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.major, ps.16s.major.rel, file = "SRKW-diet-MAJOR-16SALL.Rdata")

# ------------------------------------------------------------------
# PLOTS RELATIVE ABUNDANCE
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps16s.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

gen.rel.plot <- plot_bar(ps16s.rel, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

fam.rel.plot <- plot_bar(ps16s.rel, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot


pod.faucet <- plot_bar(ps16s.rel, x = "Sample_name", fill = "Species") +
  facet_wrap(~Pod, ncol = 4, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),  # Remove duplicate below
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines")
  )
pod.faucet


#saves plots 
ggsave("Deliverables/ALL/srkw-species.png", plot = sp.rel.plot, width = 40, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/ALL/srkw-genus.png", plot = gen.rel.plot, width = 30, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/ALL/srkw-family.png", plot = fam.rel.plot, width = 30, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# PREY ANALYSIS BY POD
# ------------------------------------------------------------------

# plots fauceted plot by pod
by.pod <- plot_bar(ps16s.rel, fill = "Species") +
  facet_wrap(~ Pod, ncol = 4, scales = "free_x", strip.position = "top") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Removes x-axis sample names
    strip.text = element_text(size = 24),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10), size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 24)
  )

by.pod



ggsave("Deliverables/ALL/srkw-by.pod.png", plot = by.pod, width = 30, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.16s.major))
colnames(otu.abs) <- as.data.frame(tax_table(ps.16s.major))$Species

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf_filt$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps.16s.major.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps.16s.major.rel))$Species

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf_filt$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
write.csv(otu.abs, "./Deliverables/ALL/SRKW_absolute_speciesxsamples-MAJOR.csv", row.names = TRUE)
write.csv(otu.prop, "./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv", row.names = TRUE)


# Calculates the total percent abundance for each prey species (across all samples)
options(scipen=999)#turns off scientific notation

# Remove the first column if it's sample IDs
df_no_id <- otu.abs[ , -1]

# Remove any rows/columns with all zeros or NAs
df_no_id <- df_no_id[rowSums(df_no_id) > 0, colSums(df_no_id) > 0]

# Calculate total abundance (sum across all samples per species)
total_abundance <- colSums(df_no_id)
total_abundance_df <- as.data.frame(total_abundance)
sum <- sum(total_abundance_df)

23348669/29075704

# Convert to percent abundance
percent_abundance <- total_abundance / sum(total_abundance_df$total_abundance)

# View result
percent_abundance_df <- as.data.frame(percent_abundance)


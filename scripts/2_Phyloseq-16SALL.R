# ------------------------------------------------------------------
# CURRENTLY WRITTEN TO JUST INCLUDE SDZWA SAMPLES
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# THIS IS THE SECOND STEP AFTER DADA2
## rownames-match.r must be run before this!
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
load("rownames-match.RData")

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

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)
lowcount.filt <- genefilter_sample(ps.16s, f1, A=1)
ps.16s.filt <- prune_taxa(lowcount.filt, ps.16s)

# must be at least 1% of diet in 4 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.16s, f1, A=4)
ps.16s.major <- prune_taxa(lowcount.filt, ps.16s)

# Saves phyloseq obj
saveRDS(ps.16s, "srkw-ps.16s.ALL")


# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.16s, fill="Species")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.16s = tax_glom(ps.16s, "Species", NArm = FALSE)

# Plots stacked bar plot of abundance - to confirm presence of NA's
abs <- plot_bar(ps.16s, fill="Species")
abs

# Saves absolute abundance plot
ggsave("Deliverables/ALL/ABSOLUTE.png", plot = abs, width = 30, height = 8, units = "in", dpi = 300)


# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.16s.filt.rel <- transform_sample_counts(ps.16s.filt, function(x) {
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

#SAVES
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps16s.rel, ps.16s, file = "SRKW-diet-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.filt, ps.16s.filt.rel, file = "SRKW-diet-FILT-16SALL.Rdata")
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
  facet_wrap(~pod, ncol = 4, scales = "free_x", strip.position = "right") +
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

# # ------------------------------------------------------------------
# # FORMATS DATA AND PLOTS RELATIVE ABUNDANCE BY PRE OR POST HORMONE ANALYSIS
# # ------------------------------------------------------------------
# 
# # Ensure sample_data is a data frame and ID columns as character
# hormone <- as(sample_data(ps.16s), "data.frame")
# hormone$ID <- as.character(hormone$ID)
# hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
# hormone$Sample_name <- as.character(hormone$Sample_name)
# 
# # Specify your chosen sample names explicitly
# samples_to_keep <- c(
#   "F24MAY27-02A-S11",
#   "F24MAY27-02A-Post-S12",
#   "F24JUN01-01C-S14",
#   "F24JUN01-01C-Post-S15",
#   "F24OCT29-06B-S20",
#   "F24OCT29-06B-Post-S21",
#   "06Oct2024-3-S25",
#   "06Oct2024-3-Post-S26",
#   "F25SEP08-03A-S16",
#   "F25SEP08-03A-Post-S17",
#   "F25NOV10-02A-S20",
#   "F25NOV10-02A-Post-S21",
#   "F25Sept16-01MS-S25",
#   "F25Sept16-01MS-Post-S26",
#   "F26JAN13-08-MS-S35",
#   "F26JAN13-08-MS-Post-S36"
# )
# 
# ps16s.abs.filtered <- prune_samples(samples_to_keep, ps.16s)
# 
# # Plot with nicer x labels
# faucet.abs <- plot_bar(ps16s.abs.filtered, fill = "Species") +
#   facet_wrap(~ Pre.Post_hormone, ncol = 1, scales = "free_x", strip.position = "right") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     panel.spacing = unit(0.5, "lines"),
#     axis.title.x = element_text(margin = margin(t = 10))
#   )
# 
# faucet.abs
# 
# # Ensure sample_data is a data frame and ID columns as character
# hormone <- as(sample_data(ps16s.rel), "data.frame")
# hormone$ID <- as.character(hormone$ID)
# hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
# hormone$Sample_name <- as.character(hormone$Sample_name)
# 
# #Checks for NaN's
# which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)
# 
# # Prune phyloseq object to retain only these samples
# ps16s.filtered <- prune_samples(samples_to_keep, ps16s.rel)
# 
# # Plot with nicer x labels
# faucet <- plot_bar(ps16s.filtered, fill = "Species") +
#   facet_wrap(~ Pre.Post_hormone, ncol = 1, scales = "free_x", strip.position = "right") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     panel.spacing = unit(0.5, "lines"),
#     axis.title.x = element_text(margin = margin(t = 10))
#   )
# 
# faucet
# 
# 
# 
# #saves plots
# ggsave("Deliverables/Plate1/srkw-pre-post-hormone.relative-all.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
# ggsave("Deliverables/Plate1/srkw-pre-post-hormone.absolute-all.png", plot = faucet.abs, width = 16, height = 8, units = "in", dpi = 300)
# 
# #SAVES
# save(seqtab.nochim, freq.nochim, track, taxam, ps16s.rel, ps.16s, file = "SRKW-diet-16SALL.Rdata")

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
otu.abs <- as.data.frame(otu_table(ps.16s))
colnames(otu.abs) <- as.data.frame(tax_table(ps.16s))$Species

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf_filt$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps16s.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps16s.rel))$Species

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
write.csv(otu.abs, "./Deliverables/ALL/SRKW_absolute_speciesxsamples.csv", row.names = TRUE)
write.csv(otu.prop, "./Deliverables/ALL/SRKW_relative_speciesxsamples.csv", row.names = TRUE)


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


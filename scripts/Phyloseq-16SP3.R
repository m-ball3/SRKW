# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
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

# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/Plate3/SRKW-diet-16SP3.Rdata")

# Gets sample metadata
samdf <- read.csv("metadata/Plate3/SRKW_Diet_Meta_Plate3.csv")

# Creates pod column from ID column
samdf <- samdf %>%
  mutate(pod = case_when(
    grepl("^J", ID) ~ "J",
    grepl("^K", ID) ~ "K",
    grepl("^L", ID) ~ "L",
    TRUE            ~ NA_character_
  ))
# ------------------------------------------------------------------
# Ensures rownames are the same
# ------------------------------------------------------------------

# Sets row names to LabID
rownames(samdf) <- samdf$Sample_ID

# Removes anythinng after underscopre for consistency
rownames(seqtab.nochim) <- sub("_.*", "", rownames(seqtab.nochim))
                             

# Replaces all dashes with underscores for consistency
rownames(seqtab.nochim) <- gsub("-", "_", rownames(seqtab.nochim))

# Checks for identical sample rownames in both
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)

# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxam))

# Creates a dataframe that maps ADFG IDs to sequences in taxa

## Converts seqtab.nochim to long format (SampleID = WADE sample ID, Sequence = sequence column names)
seqtab_long <- as.data.frame(seqtab.nochim)
seqtab_long$WADE_ID <- rownames(seqtab_long)

seqtab_long <- seqtab_long %>%
  pivot_longer(
    cols = -WADE_ID,
    names_to = "Sequence",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0)   # Keep only entries where the sequence is present in the sample

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
saveRDS(ps.16s, "srkw-ps.16s")


# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.16s, fill="Species")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.16s = tax_glom(ps.16s, "Species", NArm = FALSE)

# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.16s, fill="Species")

# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

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


pod.faucet <- plot_bar(ps16s.rel, x = "Sample_ID", fill = "Species") +
  facet_wrap(~pod, ncol = 1, scales = "free_x", strip.position = "right") +
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
ggsave("Deliverables/Plate3/srkw-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/Plate3/srkw-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/Plate3/srkw-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/Plate3/srkw-sp-bypod.png, plot = pod.faucet, width = 16, height = 8, units = "in", dpi = 300)")
# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.16s))
colnames(otu.abs) <- as.data.frame(tax_table(ps.16s))$Species

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps16s.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps16s.rel))$Species

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
write.csv(otu.abs, "./Deliverables/Plate3/SRKW_absolute_speciesxsamples.csv", row.names = FALSE)
write.csv(otu.prop, "./Deliverables/Plate3/SRKW_relative_speciesxsamples.csv", row.names = FALSE)


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
ps.16s <- merge_phyloseq(ps.raw, dna)
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

# BUILD ASV TABLE FOR MANUAL ASSIGNMENTS----------------------------------------
# # Get ASV IDs for the major taxa
# asv_ids <- taxa_names(ps.16s)
# 
# # Get sequences from ps.raw for those ASVs
# seqs <- as.character(refseq(ps.16s)[asv_ids])
# 
# # Get taxonomy for those ASVs
# tax.df <- as.data.frame(tax_table(ps.16s))
# 
# # Build the table
# asv_table <- data.frame(
#   ASV = asv_ids,
#   sequence = seqs,
#   Species = tax.df$Species,
#   stringsAsFactors = FALSE
# )
# 
# head(asv_table)
# 
# # Writes to CSV 
## PASS THIS TO 3.5_manual-taxonomy-assignment_srkw
# write.csv(asv_table, "./Deliverables/ALL/asv_table.csv", row.names = TRUE)

# USE MANUALLY CURATED TABLE TO MANUALLY ASSIGN NAs
# OUTPUT FROM 3.5_manual-taxonomy-assignment_srkw
asv_table_fixed <- read.csv("./Deliverables/ALL/asv_table-fixed.csv")

# Extract taxonomy as a data frame
tax.df <- as.data.frame(tax_table(ps.16s))

# Make sure the rownames in tax.df are ASV IDs
tax.df$ASV <- rownames(tax.df)

# Join fixed species from asv_table_fixed
tax_fixed_df <- tax.df %>%
  left_join(
    asv_table_fixed %>% select(ASV, Species_fixed),
    by = "ASV"
  )

# Use Species_fixed where available, otherwise keep original Species
tax_fixed_df$Species <- ifelse(
  is.na(tax_fixed_df$Species_fixed),
  tax_fixed_df$Species,
  tax_fixed_df$Species_fixed
)

# Drop helper columns and convert back to matrix for tax_table
tax_fixed_mat <- as.matrix(
  tax_fixed_df %>%
    select(-ASV, -Species_fixed)
)
rownames(tax_fixed_mat) <- tax_fixed_df$ASV

# Write back into ps.16s
tax_table(ps.16s) <- tax_fixed_mat

## RE-MERGE TO SPECIES HERE (TAX GLOM)
ps.16s = tax_glom(ps.16s, "Species", NArm = FALSE)

# Calculates proportional abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x >= 0.01) 

# “Minor” taxa: ≥1% in at least 1 sample
lowcount.minor <- genefilter_sample(ps16s.rel, f1, A = 1)
ps.16s.minor <- prune_taxa(lowcount.minor, ps16s.rel)

# “Major” taxa: ≥1% in at least 10 samples
lowcount.filt <- genefilter_sample(ps16s.rel, f1, A = 10)
ps.16s.major <- prune_taxa(lowcount.filt, ps16s.rel)

# Recalculates proportional abundance of each species for major and minor filtered objs
# Transforms NaN (0/0) to 0
ps.16s.minor <- transform_sample_counts(ps.16s.minor, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.16s.major <- transform_sample_counts(ps.16s.major , function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Plots stacked bar plot of abundance - to confirm presence of NA's
# plot_bar(ps.16s, fill="Species")

# Plots stacked bar plot of abundance - to confirm presence of NA's
abs <- plot_bar(ps.16s, fill="Species")
abs

# Saves relative abundance plot
ggsave("Deliverables/ALL/ABSOLUTE.png", plot = abs, width = 30, height = 8, units = "in", dpi = 300)


#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.16s.minor))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.16s.major))), arr.ind = TRUE)

#SAVES

save(ps.16s, ps16s.rel, samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "srkw-ps.16s.RData")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps16s.rel, ps.16s, file = "SRKW-diet-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.minor, file = "SRKW-diet-MINOR-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.major, file = "SRKW-diet-MAJOR-16SALL.Rdata")
save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.major, ps.16s.major, ps.16s, ps16s.rel, file = "SRKW-diet-MASTER.Rdata")

# ------------------------------------------------------------------
# PLOTS RELATIVE ABUNDANCE
# ------------------------------------------------------------------

# ALL AFTER MAMMALIA AND BACTERIA FILTERING------------------------

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



# INCLUDES FILTERING FOR MINOR SP------------------------

# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot.minor <- plot_bar(ps.16s.minor, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot.minor

gen.rel.plot.minor <- plot_bar(ps.16s.minor, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot.minor

fam.rel.plot.minor <- plot_bar(ps.16s.minor, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot.minor


# INCLUDES FILTERING FOR MAJOR SP------------------------

# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot.major <- plot_bar(ps.16s.major, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot.major

gen.rel.plot.major <- plot_bar(ps.16s.major, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot.major

fam.rel.plot.major <- plot_bar(ps.16s.major, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot.major


#saves plots 
ggsave("Deliverables/ALL/srkw-species.png", plot = sp.rel.plot, width = 40, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-species.minor.png", plot = sp.rel.plot.minor, width = 40, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-species.major.png", plot = sp.rel.plot.major, width = 40, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/ALL/srkw-genus.png", plot = gen.rel.plot, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-genus.minor.png", plot = gen.rel.plot.minor, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-genus.major.png", plot = gen.rel.plot.major, width = 30, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/ALL/srkw-family.png", plot = fam.rel.plot, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-family.minor.png", plot = fam.rel.plot.minor, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-family.major.png", plot = fam.rel.plot.major, width = 30, height = 8, units = "in", dpi = 300)

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

# plots fauceted plot by pod for minor sp
by.pod.minor <- plot_bar(ps.16s.minor, fill = "Species") +
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

by.pod.minor

# plots fauceted plot by pod for major sp
by.pod.major <- plot_bar(ps.16s.major, fill = "Species") +
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

by.pod.major

ggsave("Deliverables/ALL/srkw-by.pod.png", plot = by.pod, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-by.pod.minor.png", plot = by.pod.minor, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-by.pod.major.png", plot = by.pod.major, width = 30, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------
# CREATES ABSOLUTE SAMPLES X SPECIES TABLE
otu.abs <- as.data.frame(otu_table(ps.16s))
# colnames(otu.abs) <- as.data.frame(tax_table(ps.16s))$Species

# relabels the column names to the species name or keeps ASV name when NA
tax.df <- as.data.frame(tax_table(ps.16s))
species <- tax.df$Species
asv_ids <- taxa_names(ps.16s)

# Identify missing species names (NA or "")
missing <- is.na(species) | species == ""

# Build final labels: Species when present, ASV ID when missing
final_labels <- species
final_labels[missing] <- asv_ids[missing]

colnames(otu.abs) <- final_labels

# Writes to CSV
write.csv(otu.abs, "./Deliverables/ALL/SRKW_absolute_speciesxsamples.csv", row.names = TRUE)

#MAJOR PREY ITEMS---------------------------------------------------

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop.major <- as.data.frame(otu_table(ps.16s.major))
# colnames(otu.prop.major) <- as.data.frame(tax_table(ps.16s.major))$Species

# relabels the column names to the species name or keeps ASV name when NA
tax.df <- as.data.frame(tax_table(ps.16s.major))
species <- tax.df$Species
asv_ids <- taxa_names(ps.16s.major)

# Identify missing species names (NA or "")
missing <- is.na(species) | species == ""

# Build final labels: Species when present, ASV ID when missing
final_labels <- species
final_labels[missing] <- asv_ids[missing]

colnames(otu.prop.major) <- final_labels

# Rounds to three decimal places
is.num <- sapply(otu.prop.major, is.numeric)
otu.prop.major[is.num] <- lapply(otu.prop.major[is.num], round, 3)

# Writes to CSV
write.csv(otu.prop.major, "./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv", row.names = TRUE)


# Calculates the total percent abundance for each prey species (across all samples)
options(scipen=999)#turns off scientific notation

# # Remove the first column if it's sample IDs
# df_no_id <- otu.abs[ , -1]
# 
# # Remove any rows/columns with all zeros or NAs
# df_no_id <- df_no_id[rowSums(df_no_id) > 0, colSums(df_no_id) > 0]
# 
# # Calculate total abundance (sum across all samples per species)
# total_abundance <- colSums(df_no_id)
# total_abundance_df <- as.data.frame(total_abundance)
# sum <- sum(total_abundance_df)
# 
# 23348669/29075704

# # Convert to percent abundance
# percent_abundance <- total_abundance / sum(total_abundance_df$total_abundance)
# 
# # View result
# percent_abundance_df <- as.data.frame(percent_abundance)



#MINOR PREY ITEMS-----------------------------------------------
# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop.minor <- as.data.frame(otu_table(ps.16s.minor))

# relabels the column names to the species name or keeps ASV name when NA
tax.df <- as.data.frame(tax_table(ps.16s.minor))
species <- tax.df$Species
asv_ids <- taxa_names(ps.16s.minor)

# Identify missing species names (NA or "")
missing <- is.na(species) | species == ""

# Build final labels: Species when present, ASV ID when missing
final_labels <- species
final_labels[missing] <- asv_ids[missing]

colnames(otu.prop.minor) <- final_labels
# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop.minor, is.numeric)
otu.prop.minor[is.num] <- lapply(otu.prop.minor[is.num], round, 3)

# Writes to CSV
write.csv(otu.prop.minor, "./Deliverables/ALL/SRKW_relative_speciesxsamples-MINOR.csv", row.names = TRUE)


#ALL PREY ITEMS-----------------------------------------------
# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop.all <- as.data.frame(otu_table(ps16s.rel))

# relabels the column names to the species name or keeps ASV name when NA
tax.df <- as.data.frame(tax_table(ps16s.rel))
species <- tax.df$Species
asv_ids <- taxa_names(ps16s.rel)

# Identify missing species names (NA or "")
missing <- is.na(species) | species == ""

# Build final labels: Species when present, ASV ID when missing
final_labels <- species
final_labels[missing] <- asv_ids[missing]

colnames(otu.prop.all) <- final_labels
# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop.all, is.numeric)
otu.prop.all[is.num] <- lapply(otu.prop.all[is.num], round, 3)

# Writes to CSV
write.csv(otu.prop.all, "./Deliverables/ALL/SRKW_relative_speciesxsamples-ALL.csv", row.names = TRUE)


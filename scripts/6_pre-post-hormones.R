# ------------------------------------------------------------------
# PRE AND POST HORMONE ANALYSIS
# THIS IS THE OPTIONAL FIFTH STEP AFTER DADA2
## Phyloseq-16SALL.R must be run before this!
# ------------------------------------------------------------------

getRversion()

# Sets up environment
library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)

# loads
load("rarefaction.RData")
packageVersion("vegan")

# ------------------------------------------------------------------
# FORMATS DATA AND PLOTS RELATIVE ABUNDANCE BY PRE OR POST HORMONE ANALYSIS
# ------------------------------------------------------------------

# Ensure sample_data is a data frame and ID columns as character
hormone <- as(sample_data(ps_rarefied), "data.frame")
hormone$ID <- as.character(hormone$ID)
hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
hormone$Sample_name <- as.character(hormone$Sample_name)

# Specify your chosen sample names explicitly
samples_to_keep <- c(
  "F24MAY27-02A-S11",
  "F24MAY27-02A-Post-S12",
  "F24JUN01-01C-S14",
  "F24JUN01-01C-Post-S15",
  "F24OCT29-06B-S20",
  "F24OCT29-06B-Post-S21",
  "06Oct2024-3-S25",
  "06Oct2024-3-Post-S26",
  "F25SEP08-03A-S16",
  "F25SEP08-03A-Post-S17",
  "F25NOV10-02A-S20",
  "F25NOV10-02A-Post-S21",
  "F25Sept16-01MS-S25",
  "F25Sept16-01MS-Post-S26",
  "F26JAN13-08-MS-S35",
  "F26JAN13-08-MS-Post-S36"
)

ps16s.abs.filtered <- prune_samples(samples_to_keep, ps.16s)

# Creates a df that has the sample name, and then the pre or post (for a side by side comparison on the stacked bar plots)
# Extract sample_data as a data.frame
sd <- data.frame(sample_data(ps16s.abs.filtered))

# Remove underscore and everything after it
sd$Sample_name <- sub("_.*", "", sd$Sample_name)
sd$Sample_name <- sub("-.*", "", sd$Sample_name)

# Assign back to the phyloseq object
sample_data(ps16s.abs.filtered) <- sample_data(sd)

# creates a vector of custom colors for the stacked bar plot
# Manual color palette based on the graph/reference style:
# warm pink/magenta-heavy palette + a few orange/green/blue accents
species_cols <- c(
  "Anarrhichthys ocellatus" = "#F28E8C",
  "Anoplopoma fimbria" = "#249FE8",
  "Atheresthes stomias" = "#00E381",
  "Citharichthys sordidus" = "#CDA000",
  "Citharichthys stigmaeus" = "#B8A200",
  "Clupea pallasii" = "#9FAA00",
  "Cymatogaster aggregata" = "#7FAE00",
  "Diaphus theta" = "#57B000",
  "Engraulis mordax" = "#18B600",
  "Gadus chalcogrammus" = "#08BA4A",
  "Gasterosteus aculeatus" = "#10BE73",
  "Glyptocephalus zachirus" = "#16C08E",
  "Hippoglossus stenolepis" = "#19BEA7",
  "Hydrolagus colliei" = "#1ABBB9",
  "Liparis pulchellus" = "#18B2D9",
  "Microstomus pacificus" = "#1CA7E0",
  "Myoxocephalus polyacanthocephalus" = "#249FE8",
  "Oncorhynchus gorbuscha" = "#4C95EA",
  "Oncorhynchus keta" = "#008F51",
  "Oncorhynchus kisutch" = "#091A80",
  "Oncorhynchus mykiss" = "#B36CE2",
  "Oncorhynchus nerka" = "#11A5FF",
  "Oncorhynchus tshawytscha" = "#1233FF",
  "Ophiodon elongatus" = "#005933",
  "Raja binoculata" = "#EC5FB9",
  "Raja rhina" = "#F062A5",
  "Salmo salar" = "#F46A8A"
)

#3B3B3B

# Plot with nicer x labels
faucet.abs <- plot_bar(ps16s.abs.filtered, fill = "Species") +
  facet_wrap(~ Sample_name, ncol = 20, scales = "free_x", strip.position = "top") +
  scale_fill_manual(
    values = species_cols,
    na.value = "grey60",
    drop = FALSE
  )+
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(), 
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none"
  )

faucet.abs

# Ensure sample_data is a data frame and ID columns as character
hormone <- as(sample_data(ps_rarefied.rel), "data.frame")
hormone$ID <- as.character(hormone$ID)
hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
hormone$Sample_name <- as.character(hormone$Sample_name)

#Checks for NaN's
which(is.nan(as.matrix(otu_table(ps_rarefied.rel))), arr.ind = TRUE)

# Prune phyloseq object to retain only these samples
ps16s.rel.filtered <- prune_samples(samples_to_keep, ps_rarefied.rel)

# Creates a df that has the sample name, and then the pre or post (for a side by side comparison on the stacked bar plots)
# Extract sample_data as a data.frame
sd2 <- data.frame(sample_data(ps16s.rel.filtered))

# Remove underscore and everything after it
sd2$Sample_name <- sub("_.*", "", sd2$Sample_name)
sd2$Sample_name <- sub("-.*", "", sd2$Sample_name)

# Assign back to the phyloseq object
sample_data(ps16s.rel.filtered) <- sample_data(sd2)

## COULD ADD A COLUMN TO SAMDF OR WHATEVER THAT STRIPS THE COMMON NAME FROM BOTH, MAKES A PRE AND POST
## AND THEN THE PLOT WILL HAVE THEM SIDE BY SIDE GROUPED TOGETHER LIKE IN THE REPLICATES PLOT IN ADFG
# Plot with nicer x labels
faucet <- plot_bar(ps16s.rel.filtered, fill = "Species") +
  facet_wrap(~ Sample_name, ncol = 20, scales = "free_x", strip.position = "top") +
  scale_fill_manual(
    values = species_cols,
    na.value = "grey60",
    drop = FALSE
  )+
  ylab("Proportion") + 
  theme_bw() +
  theme(
    text = element_text(size = 14), 
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text.x     = element_blank(), 
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )

faucet

# modifies the x axis sample labels to just say pre or post
# Extract the data used by plot_bar
pd <- faucet$data

# Build a named character vector: names = original Sample values, values = "Pre"/"Post"
sample_labels <- c(
  # F24MAY27
  "F24MAY27-02A-S11"        = "Pre",
  "F24MAY27-02A-Post-S12"   = "Post",
  
  # F24JUN01
  "F24JUN01-01C-S14"        = "Pre",
  "F24JUN01-01C-Post-S15"   = "Post",
  
  # F24OCT29
  "F24OCT29-06B-S20"        = "Pre",
  "F24OCT29-06B-Post-S21"   = "Post",
  
  # 06Oct2024
  "06Oct2024-3-S25"         = "Pre",
  "06Oct2024-3-Post-S26"    = "Post",
  
  # F25SEP08.03A
  "F25SEP08-03A-S16"        = "Pre",
  "F25SEP08-03A-Post-S17"   = "Post",
  
  # F25NOV10.02A
  "F25NOV10-02A-S20"        = "Pre",
  "F25NOV10-02A-Post-S21"   = "Post",
  
  # F25Sept16.01MS
  "F25Sept16-01MS-S25"      = "Pre",
  "F25Sept16-01MS-Post-S26" = "Post",
  
  # F26JAN13.08.MS
  "F26JAN13-08-MS-S35"      = "Pre",
  "F26JAN13-08-MS-Post-S36" = "Post"
)

# Apply labels to the x-axis
faucet <- faucet +
  scale_x_discrete(labels = sample_labels)

faucet

combined <- faucet.abs / faucet
combined

#saves plots
ggsave("Deliverables/ALL/srkw-pre-post-hormone.relative-all.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-pre-post-hormone.absolute-all.png", plot = faucet.abs, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-pre-post-hormone.stacked.png", plot = combined, width = 16, height = 8, units = "in", dpi = 300)


# # STATS : BRAY-CURTIS PERMANOVA ON RELATIVE ABUNDANCES
# FROM WIKIPEDIA: The null hypothesis states that the centroids (averages) and/or the dispersion (spread) of the groups are equivalent in the multivariate space.
## H0 - there is no difference in the pre samples vs post samples
## HA - there is a difference in the pre samples vs post samples

# filters samdf_filt
samdf_horm <- samdf_filt[rownames(samdf_filt) %in% samples_to_keep, ]


# Calculate a distance matrix (Bray-Curtis)
bray_matrix <- phyloseq::distance(ps16s.filtered, method = "bray")

# Run the PERMANOVA
permanova_bray. <- adonis2(bray_matrix ~ Pre.Post_hormone, data = samdf_horm, permutations = 999)

print(permanova_bray)


# STATS : JACCARD PERMANOVA ON ABSOLUTE ABUNDANCES
jac_matrix <- phyloseq::distance(ps16s.abs.filtered, method = "jaccard")

# Run the PERMANOVA

# Calculate a distance matrix (Jaccard)
permanova_jac <- adonis2(jac_matrix ~ Pre.Post_hormone, data = samdf_horm, permutations = 999)

print(permanova_jac)

# STATS : JACCARD PERMANOVA ON RELATIVE ABUNDANCES

# Calculate a distance matrix (Jaccard)
jac_matrix.rel <- phyloseq::distance(ps16s.rel.filtered, method = "jaccard")

# Run the PERMANOVA
permanova_jac.rel <- adonis2(jac_matrix.rel ~ Pre.Post_hormone, data = samdf_horm, permutations = 999)

print(permanova_jac.rel)





##### NEW CODE BUT NOT READY YET
library(phyloseq)
library(ggplot2)
library(grid)

# Keep sample order exactly as you want it shown
samples_to_keep.test <- c(
  "F24MAY27-02A-S11",
  "F24MAY27-02A-Post-S12",
  "F24JUN01-01C-S14",
  "F24JUN01-01C-Post-S15",
  "F24OCT29-06B-S20",
  "F24OCT29-06B-Post-S21",
  "06Oct2024-3-S25",
  "06Oct2024-3-Post-S26",
  "F25SEP08-03A-S16",
  "F25SEP08-03A-Post-S17",
  "F25NOV10-02A-S20",
  "F25NOV10-02A-Post-S21",
  "F25Sept16-01MS-S25",
  "F25Sept16-01MS-Post-S26",
  "F26JAN13-08-MS-S35",
  "F26JAN13-08-MS-Post-S36"
)

# Filter objects
ps16s.abs.filtered.test <- prune_samples(samples_to_keep, ps.16s)
ps16s.rel.filtered.test <- prune_samples(samples_to_keep, ps_rarefied.rel)

# Force sample order in metadata
sample_data(ps16s.abs.filtered.test)$Sample_name <- factor(
  sample_data(ps16s.abs.filtered.test)$Sample_name,
  levels = samples_to_keep.test
)

sample_data(ps16s.rel.filtered.test)$Sample_name <- factor(
  sample_data(ps16s.rel.filtered.test)$Sample_name,
  levels = samples_to_keep.test
)



# Common theme
base_bar_theme.test <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Absolute abundance plot
faucet.abs.test <- plot_bar(ps16s.abs.filtered.test, x = "Sample_name", fill = "Species") +
  facet_wrap(~ Pre.Post_hormone, ncol = 1, scales = "free_x", strip.position = "right") +
  scale_fill_manual(
    values = species_cols.test,
    na.value = "grey60",
    drop = FALSE
  ) +
  labs(x = "Sample", y = "Abundance") +
  base_bar_theme.test

faucet.abs.test

# Relative abundance plot
faucet.test <- plot_bar(ps16s.rel.filtered.test, x = "Sample_name", fill = "Species") +
  facet_wrap(~ Pre.Post_hormone, ncol = 2, scales = "free_x", strip.position = "top") +
  scale_fill_manual(
    values = species_cols.test,
    na.value = "grey60",
    drop = FALSE
  ) +
  labs(x = "Sample", y = "Abundance") +
  base_bar_theme.test

faucet.test


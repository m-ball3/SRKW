# ------------------------------------------------------------------
# PRE AND POST HORMONE ANALYSIS
# THIS IS THE OPTIONAL FIFTH STEP AFTER DADA2
## Phyloseq-16SALL.R must be run before this!
# ------------------------------------------------------------------

# Sets up environment
library(phyloseq)
library(tidyverse)
library(vegan)

# loads
load("srkw-ps.16s.Rdata")

# ------------------------------------------------------------------
# FORMATS DATA AND PLOTS RELATIVE ABUNDANCE BY PRE OR POST HORMONE ANALYSIS
# ------------------------------------------------------------------

# Ensure sample_data is a data frame and ID columns as character
hormone <- as(sample_data(ps.16s), "data.frame")
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

# Plot with nicer x labels
faucet.abs <- plot_bar(ps16s.abs.filtered, fill = "Species") +
  facet_wrap(~ Pre.Post_hormone, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  )

faucet.abs

# Ensure sample_data is a data frame and ID columns as character
hormone <- as(sample_data(ps16s.rel), "data.frame")
hormone$ID <- as.character(hormone$ID)
hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
hormone$Sample_name <- as.character(hormone$Sample_name)

#Checks for NaN's
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

# Prune phyloseq object to retain only these samples
ps16s.filtered <- prune_samples(samples_to_keep, ps16s.rel)

# Plot with nicer x labels
faucet <- plot_bar(ps16s.filtered, fill = "Species") +
  facet_wrap(~ Pre.Post_hormone, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  )

faucet

#saves plots
ggsave("Deliverables/Plate1/srkw-pre-post-hormone.relative-all.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/Plate1/srkw-pre-post-hormone.absolute-all.png", plot = faucet.abs, width = 16, height = 8, units = "in", dpi = 300)


# STATS : PERMANOVA?
# FROM WIKIPEDIA: The null hypothesis states that the centroids (averages) and/or the dispersion (spread) of the groups are equivalent in the multivariate space.
## H0 - there is no difference in the pre samples vs post samples
## HA - there is a difference in the pre samples vs post samples

# filters samdf_filt
samdf_horm <- samdf_filt[rownames(samdf_filt) %in% samples_to_keep, ]

# prunes the phyloseq obj
ps.16s_horm <- prune_samples(
  sample_names(ps.16s) %in% samples_to_keep,
  ps.16s
)


# Assuming 'community_data' is your species/microbiome table 
# and 'metadata' contains your grouping variable

# Calculate a distance matrix (e.g., Bray-Curtis)
dist_matrix <- phyloseq::distance(ps.16s_horm, method = "bray")

# Run the PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Pre.Post_hormone, data = samdf_horm, permutations = 999)

print(permanova_result)


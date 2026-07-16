# ------------------------------------------------------------------
# PRE AND POST HORMONE ANALYSIS
# THIS IS THE OPTIONAL FIFTH STEP AFTER DADA2
## 3_phyloseq_srkw must be run before this!
## 5_rarefaction_srkw.R optional to run before this, but should not 
## be used for clr and Aitchison distance!!
# ------------------------------------------------------------------

getRversion()

# Sets up environment
library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)
library(robCompositions)

# loads
load("srkw-ps.16s.Rdata")
packageVersion("vegan")

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

# -----------------------------
# ABSOLUTE ABUNDANCE PHYLOSEQ
# -----------------------------

ps16s.abs.filtered <- prune_samples(samples_to_keep, ps.16s)

# Creates a df that has the sample name, and then the pre or post (for a side by side comparison on the stacked bar plots)
# Extract sample_data as a data.frame
sd_abs <- data.frame(sample_data(ps16s.abs.filtered))

# Remove underscore and everything after it
sd_abs$Sample_name <- sub("_.*", "", sd_abs$Sample_name)
sd_abs$Sample_name <- sub("-.*", "", sd_abs$Sample_name)

# renames samples to 1, 2, 3...
# Order samples in the exact sequence you want
sd_abs <- sd_abs[match(samples_to_keep, rownames(sd_abs)), ]

# Assign an index 1..n per paired pre/post sample.
# Assuming two rows per Sample_name, you want an index per unique Sample_name:
sd_abs$Sample_index <- as.integer(factor(sd_abs$Sample_name))  # 1,2,3,...

# makes it so that Pre comes before Post
sd_abs$PrePost <- factor(sd_abs$Pre.Post_hormone,
                         levels = c("Pre", "Post"))  # enforce Pre then Post

# Assign back to the phyloseq object
sample_data(ps16s.abs.filtered) <- sample_data(sd_abs)

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
faucet.abs <- plot_bar(ps16s.abs.filtered,
                       x = "PrePost", fill = "Species") +
  facet_wrap(~ Sample_index, ncol = 20,
             scales = "free_x", strip.position = "top") +
  scale_fill_manual(
    values = species_cols,
    na.value = "grey60",
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    text          = element_text(size = 14),
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    legend.position  = "none"
  )

faucet.abs

# -----------------------------
# RELATIVE ABUNDANCE PHYLOSEQ
# -----------------------------

# Ensure sample_data is a data frame and ID columns as character
hormone <- as(sample_data(ps16s.rel), "data.frame")
hormone$ID <- as.character(hormone$ID)
hormone$Pre.Post_hormone <- as.character(hormone$Pre.Post_hormone)
hormone$Sample_name <- as.character(hormone$Sample_name)

#Checks for NaN's
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

# Prune phyloseq object to retain only these samples
ps16s.rel.filtered <- prune_samples(samples_to_keep, ps16s.rel)

# Creates a df that has the sample name, and then the pre or post (for a side by side comparison on the stacked bar plots)
# Extract sample_data as a data.frame
sd_rel <- data.frame(sample_data(ps16s.rel.filtered))

# Remove underscore and everything after it
sd_rel$Sample_name <- sub("_.*", "", sd_rel$Sample_name)
sd_rel$Sample_name <- sub("-.*", "", sd_rel$Sample_name)

# Order samples in the exact sequence you want
sd_rel <- sd_rel[match(samples_to_keep, rownames(sd_rel)), ]

# Assign an index 1..n per paired pre/post sample.
sd_rel$Sample_index <- as.integer(factor(sd_rel$Sample_name))  # 1,2,3,...

# makes it so that Pre comes before Post
sd_rel$PrePost <- factor(sd_rel$Pre.Post_hormone,
                         levels = c("Pre", "Post"))  # enforce Pre then Post

# Assign back to the phyloseq object
sample_data(ps16s.rel.filtered) <- sample_data(sd_rel)

## COULD ADD A COLUMN TO SAMDF OR WHATEVER THAT STRIPS THE COMMON NAME FROM BOTH, MAKES A PRE AND POST
## AND THEN THE PLOT WILL HAVE THEM SIDE BY SIDE GROUPED TOGETHER LIKE IN THE REPLICATES PLOT IN ADFG
# Plot with nicer x labels
faucet <- plot_bar(ps16s.rel.filtered,
                   x = "PrePost",
                   fill = "Species") +
  facet_wrap(~ Sample_index, ncol = 20,
             scales = "free_x", strip.position = "top") +
  scale_fill_manual(
    values = species_cols,
    na.value = "grey60",
    drop = FALSE
  ) +
  xlab("Sample") +
  ylab("Proportion") +
  theme_bw() +
  theme(
    text            = element_text(size = 14),
    strip.background = element_blank(),
    strip.text.x     = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10)),
    legend.position  = "none"
  )

faucet

combined <- faucet.abs / faucet
combined

#saves plots
ggsave("Deliverables/ALL/srkw-pre-post-hormone.relative-all-NORARERIFY.png",plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-pre-post-hormone.absolute-all-NORARERIFY.png",plot = faucet.abs, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ALL/srkw-pre-post-hormone.stacked-NORARERIFY.png",plot = combined, width = 16, height = 8, units = "in", dpi = 300)

# # STATS : AITCHISON PERMANOVA ON RELATIVE ABUNDANCES 
# # NOT RAREFIED; USED CLR METHOD PER GLOOR ET AL., 2017
# FROM WIKIPEDIA: The null hypothesis states that the centroids (averages) and/or the dispersion (spread) of the groups are equivalent in the multivariate space.
## H0 - there is no difference in the pre samples vs post samples
## HA - there is a difference in the pre samples vs post samples

# filters samdf_filt
samdf_horm <- samdf_filt[rownames(samdf_filt) %in% samples_to_keep, ]

# creates a new column in samdf_horm that has the Sample_name with the suffix dropped
## this is so we can identify strata in the permanova
### meaning that we can account for the paired nature of the data
samdf_horm$samps <- samdf_horm$Sample_name |>
  sub("_Post.*", "", x = _) |>
  sub("_S[0-9]+", "", x = _) |>
  sub("-Post", "", x = _)

# Extract OTU/ASV table as a plain matrix (for passing to vegdist)
otu_mat <- as(otu_table(ps16s.rel.filtered), "matrix")  # or as.matrix(otu_table(...))

# Decide orientation: vegdist expects rows = samples, columns = variables
# If taxa_are_rows(ps16s.rel.filtered) is TRUE, transpose:
if (taxa_are_rows(ps16s.rel.filtered)) {
  otu_mat <- t(otu_mat)
}

# Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"
robust_aitch_dist <- vegdist(otu_mat, method = "robust.aitchison")

# Run the PERMANOVA
permanova_aitch <- adonis2(robust_aitch_dist ~ Pre.Post_hormone, 
                           data = samdf_horm, 
                           strata = samdf_horm$samps, 
                           permutations = 999)

print(permanova_aitch)


# 
# # OPTIONAL PCA -------------------------------------------------------------------------------------------
# 
# # rCLR transform (manually creates the matrix that vegdist creates)
# otu_rclr <- decostand(otu_mat, method = "rclr")
# 
# # Runs PCA on rCLR-transformed data
# pca_res <- prcomp(otu_rclr, center = TRUE, scale. = FALSE)
# pca_res
# 
# # Scores (sample coordinates in PCA space)
# pca_scores <- pca_res$x
# pca_scores
# 
# # Loadings (taxon contributions)
# pca_loadings <- pca_res$rotation
# pca_loadings
# 
# 
# # PCA scores (samples × PCs)
# pca_scores <- as.data.frame(pca_res$x)
# 
# # Add sample IDs as a column
# pca_scores$Sample <- rownames(pca_scores)
# 
# # Pull metadata from your phyloseq object
# meta <- data.frame(sample_data(ps16s.rel.filtered))
# 
# meta$Sample <- rownames(meta)
# 
# # Join scores with metadata
# pca_plot_df <- pca_scores %>%
#   left_join(meta, by = "Sample")
# 
# var_explained <- summary(pca_res)$importance[2, ]  # proportion per PC [web:244]
# pc1_pct <- round(var_explained[1] * 100, 1)
# pc2_pct <- round(var_explained[2] * 100, 1)
# 
# p_pca <- ggplot(pca_plot_df,
#                 aes(x = PC1, y = PC2,
#                     color = Pre.Post_hormone)) +
#   geom_point(size = 3, alpha = 0.9) +
#   labs(
#     x = paste0("PC1 (", pc1_pct, "% variance)"),
#     y = paste0("PC2 (", pc2_pct, "% variance)"),
#     color = "Hormone status"
#   ) +
#   theme_bw() +
#   theme(
#     text          = element_text(size = 14),
#     panel.grid    = element_blank(),
#     legend.position = "right"
#   )
# 
# p_pca
# 
# 
# p_pca_faceted <- ggplot(pca_plot_df,
#                         aes(x = PC1, y = PC2,
#                             color = Pre.Post_hormone)) +
#   geom_point(size = 3, alpha = 0.9) +
#   facet_wrap(~ Sample_index, ncol = 4) +
#   labs(
#     x = paste0("PC1 (", pc1_pct, "% variance)"),
#     y = paste0("PC2 (", pc2_pct, "% variance)"),
#     color = "Hormone status"
#   ) +
#   theme_bw() +
#   theme(
#     text            = element_text(size = 12),
#     strip.background = element_blank(),
#     strip.placement  = "outside",
#     panel.grid       = element_blank(),
#     legend.position  = "bottom"
#   )
# 
# p_pca_faceted

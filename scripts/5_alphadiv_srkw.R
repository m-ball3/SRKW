# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS ALPHA DIVERSITY METRICS ON SRKW DATA
# THIS IS THE FIRST STATISTICAL INVESTIGATION 
## ALL P VALUES REPORTED ARE THE ADJUSTED P!!
# ------------------------------------------------------------------------------

# install.packages(
#   "rbiom",
#   type = "binary",
#   repos = "https://cran.rstudio.com"
# )

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(tibble)
library(stringr)
library(phyloseq)
library(car)
library(vegan)
library(rstatix)
library(Hmisc)
library(rbiom)
library(mia) # via bioconductor
library(patchwork)

getwd()
load("./SRKW-diet-MASTER.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.16s.major), "matrix")

# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------

# formatting for metadata
metadata <- data.frame(sample_data(ps.16s.major))

class(metadata)

# optionally removes NA in pod_from_id
metadata <- metadata %>%
  dplyr::filter(!is.na(pod_from_id))

# optionally removes NA in season
metadata <- metadata %>%
  dplyr::filter(!is.na(season))

# makes rownames a column
metadata <- metadata %>%
  rownames_to_column("samp")

# makes asv_mat match
asv_mat <- as.data.frame(asv_mat) %>%
  rownames_to_column(var = "samp") %>%
  dplyr::filter(samp %in% metadata$samp) %>%
  column_to_rownames(var = "samp")

metadata %>%
  dplyr::count(pod_from_id, name = "n_samples") %>%
  dplyr::arrange(pod_from_id)

# sets host colors
pod_colors <- c(
  "J" = "#F8766D",  # ggplot2 default red
  "K" = "#00BA38",  # ggplot2 default green
  "L"  = "#619CFF"   # ggplot2 default blue
)

# ------------------------------------------------------------------------------
# Gets richness, dominance, and information metrics for alpha diversity
# ------------------------------------------------------------------------------

alpha_div <- data.frame(
  samp = rownames(asv_mat),
  
  # Observed ASV richness
  Observed = rowSums(asv_mat > 0),
  
  # Berger-Parker dominance: dominant ASV reads / total sample reads
  Berger_Parker = apply(asv_mat, 1, max) / rowSums(asv_mat),
  
  # Shannon diversity
  Shannon = vegan::diversity(asv_mat, index = "shannon"),
  
  # Useful QC covariate
  n_reads = rowSums(asv_mat),
  
  row.names = NULL
)

# adds alpha_div metrics to the metadata tables
metadata <- metadata %>%
  dplyr::left_join(alpha_div, by = "samp")

# ------------------------------------------------------------------------------
# Richness
# ------------------------------------------------------------------------------

# Robbins
## the ratio of singletons to total taxa


# Observed 
## the number of ASVs present

hist(metadata$Observed)

## wondering if I should run a Kruskal-Wallis test to test for actual differences here

# test for normalcy
norm <- shapiro.test(metadata$Observed)
norm # p-value = 2.296e-11 (significant - not normally distributed)

# test homogeneity of variance FOR HOST
#?leveneTest()
var <- leveneTest(Observed ~ season, metadata)
var # p-value = 0.01366 (significant - non-homonegeous variance)
# Does not meet assumptions for parametric tests

# tests for correlations between library depth and observed richness
cor.test(
  metadata$Observed,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
)

ggplot(metadata, aes(x = n_reads, y = Observed, color = season)) +
  geom_point(alpha = 0.75) +
  scale_x_log10() +
  theme_classic() +
  labs(
    x = "Retained reads per sample (log10 scale)",
    y = "Observed ASVs"
  )

# p-value = 0.2469; though there is a slightly positive trend, there is no evidence
# that the read depth is highly correlated with the observed richness
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)



# BY SEASON----------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_season <- kruskal.test(Observed ~ season, data = metadata)
obs_test_season #  p-value = 2.2e-16 (significant)

# Dunn test to determine which groups are significantly different
#?dunn_test()
obs_dunn_results_season <- metadata %>%
  dunn_test(
    Observed ~ season,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_season

# Add x- and y-coordinates for the comparison brackets
obs_dunn_plot_season <- obs_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
obs_dunn_plot_season

# Box-and-whisker plot, with individual samples overlaid
# BY SEASON
obs_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Observed, fill = pod_from_id)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    obs_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = pod_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  labs(
    x = "Season",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by season, colored by pod"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom", 
        # Text sizes in points
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 12)),
        axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 12)),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        
        # Make ticks and plot lines readable on slides
        axis.line = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.25, "cm")
  )

obs_plot_season

ggsave(
  "Deliverables/ALL/alphadiv/OBS-SEASON_major.png",
  plot = obs_plot_season,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)

# BY POD----------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_pod <- kruskal.test(Observed ~ pod_from_id, data = metadata)
obs_test_pod #  p-value = 0.2654 (non-significant)

# Dunn test to determine which groups are significantly different
#?dunn_test()
obs_dunn_results_pod <- metadata %>%
  dunn_test(
    Observed ~ pod_from_id,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_pod

# Add x- and y-coordinates for the comparison brackets
obs_dunn_plot_pod <- obs_dunn_results_pod %>%
  add_xy_position(x = "pod_from_id")

# Inspect the table used for plotting
obs_dunn_plot_pod

# Box-and-whisker plot, with individual samples overlaid
# BY SEASON
obs_plot_pod <- ggplot(
  metadata,
  aes(x = pod_from_id, y = Observed, fill = pod_from_id)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    obs_dunn_plot_pod,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = pod_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  labs(
    x = "Pod",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by pod, colored by pod"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom", 
        # Text sizes in points
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 12)),
        axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 12)),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        
        # Make ticks and plot lines readable on slides
        axis.line = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.25, "cm")
  )

obs_plot_pod

ggsave(
  "Deliverables/ALL/alphadiv/OBS-POD_major.png",
  plot = obs_plot_pod,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# Separates into Seasons
## DO THIS ONLY IF POD IS SIGNIFICANT
# ------------------------------------------------------------------------------

summer_metadata <- metadata %>%
  dplyr::filter(season == "Summer")
summer_asv_mat <- asv_mat %>%
  rownames_to_column(var = "samp") %>%
  dplyr::filter(samp %in% summer_metadata$samp) %>%
  column_to_rownames(var = "samp")

autumn_metadata <- metadata %>%
  dplyr::filter(season == "Autumn")
autumn_asv_mat <- asv_mat %>%
  rownames_to_column(var = "samp") %>%
  dplyr::filter(samp %in% autumn_metadata$samp) %>%
  column_to_rownames(var = "samp")

spring_metadata <- metadata %>%
  dplyr::filter(season == "Spring")
spring_asv_mat <- asv_mat %>%
  rownames_to_column(var = "samp") %>%
  dplyr::filter(samp %in% spring_metadata$samp) %>%
  column_to_rownames(var = "samp")

winter_metadata <- metadata %>%
  dplyr::filter(season == "Winter")
summer_asv_mat <- asv_mat %>%
  rownames_to_column(var = "samp") %>%
  dplyr::filter(samp %in% winter_metadata$samp) %>%
  column_to_rownames(var = "samp")

# ------------------------------------------------------------------------------
# Separates into Pods
## DO THIS ONLY IF POD IS SIGNIFICANT
# ------------------------------------------------------------------------------
# 
# J_metadata <- metadata %>%
#   dplyr::filter(pod_from_id == "J")
# J_asv_mat <- asv_mat %>%
#   rownames_to_column(var = "samp") %>%
#   dplyr::filter(samp %in% J_metadata$samp) %>%
#   column_to_rownames(var = "samp")
# 
# K_metadata <- metadata %>%
#   dplyr::filter(pod_from_id == "K")
# K_asv_mat <- asv_mat %>%
#   rownames_to_column(var = "samp") %>%
#   dplyr::filter(samp %in% K_metadata$samp) %>%
#   column_to_rownames(var = "samp")
# 
# L_metadata <- metadata %>%
#   dplyr::filter(pod_from_id == "L")
# L_asv_mat <- asv_mat %>%
#   rownames_to_column(var = "samp") %>%
#   dplyr::filter(samp %in% L_metadata$samp) %>%
#   column_to_rownames(var = "samp")


# ------------------------------------------------------------------------------
# Dominance/Evenness
# ------------------------------------------------------------------------------

# Berger-Parker is a ratio between the abundance of the most abundant ASV and the 
# total abundance. 
# It is the within-sample relative abundance of the single most abundant ASV.
# Higher values indicate stronger dominance and lower community evenness.

# Compared with other dominance formula metrics, it is easy to calculate and has a 
# clear biological interpretation.

# test for normalcy
bp.norm <- shapiro.test(metadata$Berger_Parker)
bp.norm # p-value =  2.2e-16 (significant; nonnormal)

# test homogeneity of variance
#?leveneTest()
bp.var <- leveneTest(Berger_Parker ~ season, metadata)
bp.var # p-value = 2.2e-16 (non-signficant - homogeneous variance)

# Does meet assumptions for parametric tests

# tests if read depth is correlated with Berger-Parker ratios
cor.test(
  metadata$Berger_Parker,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
) 

# p-value = 0.1459; though there is a slight positive trend, there is no evidence
# that the read depth is highly correlated with the Berger-Parker ratio
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)

# BY HOST---------------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_season <- kruskal.test(Berger_Parker ~ season, data = metadata)
bp_test_season #  p-value = 2.2e-16 significant 

# Dunn test to determine which groups are significantly different
bp_dunn_results_season <- metadata %>%
  dunn_test(
    Berger_Parker ~ season,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_season

### SIGNIFICANTLY DIFFERENT DOMINANCE


# Add x- and y-coordinates for the comparison brackets
# Compute non-overlapping bracket positions
bp_dunn_plot_season <- bp_dunn_results_season %>%
  dplyr::filter(p.adj <= 0.05) %>%
  rstatix::add_x_position(x = "season") %>%
  dplyr::arrange(group1, group2) %>%
  dplyr::mutate(
    y.position = 1 + seq(0.05, by = 0.07, length.out = dplyr::n())
  )

bp_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Berger_Parker, fill = pod_from_id)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    bp_dunn_plot_season,
    label = "p.adj.signif",
    tip.length = 0.01,
    size = 8,
    bracket.size = 0.8,
    vjust = 0.2
  ) +
  scale_fill_manual(
    values = pod_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Season",
    y = "Dominance of ASVs per sample",
    title = "Berger–Parker dominance by season, colored by pod"
  ) +
  scale_y_continuous(
    limits = c(0, 1 + 0.07 * nrow(bp_dunn_plot_season) + 0.08),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 12)),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.line = element_line(linewidth = 1.1),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length = unit(0.25, "cm")
  )

bp_plot_season

ggsave(
  "Deliverables/ALL/alphadiv/BP-SEASON_major.png",
  plot = bp_plot_season,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)

# BY POD----------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_pod <- kruskal.test(Berger_Parker ~ pod_from_id, data = metadata)
bp_test_pod #  p-value = 0.2654 (non-significant)

# Dunn test to determine which groups are significantly different
#?dunn_test()
bp_dunn_results_pod <- metadata %>%
  dunn_test(
    Berger_Parker ~ pod_from_id,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_pod

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_pod <- bp_dunn_results_pod %>%
  add_xy_position(x = "pod_from_id")

# Inspect the table used for plotting
bp_dunn_plot_pod

# Box-and-whisker plot, with individual samples overlaid
# BY SEASON
bp_plot_pod <- ggplot(
  metadata,
  aes(x = pod_from_id, y = Berger_Parker, fill = pod_from_id)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    bp_dunn_plot_pod,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = pod_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  labs(
    x = "Pod",
    y = "Dominance of ASVs per sample",
    title = "Berger-Parker Dominance by pod, colored by pod"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom", 
        # Text sizes in points
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 12)),
        axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 12)),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        
        # Make ticks and plot lines readable on slides
        axis.line = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.25, "cm")
  )

bp_plot_pod

ggsave(
  "Deliverables/ALL/alphadiv/BP-POD_major.png",
  plot = bp_plot_pod,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)

# DOMINANT ASV -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# GET TAXONOMY AND DOMINANT ASV FOR EACH SAMPLE
# ------------------------------------------------------------------------------

# Extract taxonomy table and retain ASV IDs as a column
tax_df <- as.data.frame(phyloseq::tax_table(ps.16s.major)) %>%
  tibble::rownames_to_column("ASV")

# Identify the column index of the most abundant ASV for each sample.
# ties.method = "first" consistently selects the first ASV if two ASVs tie.
dominant_idx <- max.col(asv_mat, ties.method = "first")

# Create a sample-level table containing each sample's dominant ASV.
dominant_asv_per_sample <- tibble::tibble(
  samp = rownames(asv_mat),
  ASV = colnames(asv_mat)[dominant_idx],
  n_reads_asv = asv_mat[
    cbind(seq_len(nrow(asv_mat)), dominant_idx)
  ],
  sample_reads = rowSums(asv_mat)
) %>%
  dplyr::mutate(
    dominance = n_reads_asv / sample_reads
  ) %>%
  dplyr::left_join(
    tax_df,
    by = "ASV"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        samp,
        pod_from_id,
        season,
      ),
    by = "samp"
  ) %>%
  dplyr::select(
    samp,
    pod_from_id,
    season,
    ASV,
    n_reads_asv,
    sample_reads,
    dominance,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::arrange(dplyr::desc(dominance))

dominant_asv_per_sample

# ------------------------------------------------------------------------------
# CREATE READABLE DOMINANT-ASV LABELS
# ------------------------------------------------------------------------------

# ASV sequences or ASV IDs are often difficult to read in a figure legend.
# This label uses genus where available, otherwise Family, Order, Class, or
# Phylum, and appends the ASV ID to distinguish multiple ASVs in the same taxon.

dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::mutate(
    ASV_short = substr(ASV, 1, 12),
    
    taxon_label = dplyr::case_when(
      !is.na(Species) &
        Species != "" &
        !tolower(Species) %in% c("na", "unknown", "unclassified") ~
        paste0("s__ ", Species, " | ", ASV_short),
      
      !is.na(Genus) &
        Genus != "" &
        !tolower(Genus) %in% c("na", "unknown", "unclassified") ~
        paste0("g__ ", Genus, " | ", ASV_short),
      
      !is.na(Family) &
        Family != "" &
        !tolower(Family) %in% c("na", "unknown", "unclassified") ~
        paste0("f__ ", Family, " | ", ASV_short),
      
      !is.na(Order) &
        Order != "" &
        !tolower(Order) %in% c("na", "unknown", "unclassified") ~
        paste0("o__ ", Order, " | ", ASV_short),
      
      !is.na(Class) &
        Class != "" &
        !tolower(Class) %in% c("na", "unknown", "unclassified") ~
        paste0("c__ ", Class, " | ", ASV_short),
      
      !is.na(Phylum) &
        Phylum != "" &
        !tolower(Phylum) %in% c("na", "unknown", "unclassified") ~
        paste0("p__ ", Phylum, " | ", ASV_short),
      
      TRUE ~ paste0("Unclassified | ", ASV_short)
    )
  )

# ------------------------------------------------------------------------------
# OVERALL DOMINANT-ASV SUMMARY
# ------------------------------------------------------------------------------

dominant_asv_summary <- dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

dominant_asv_summary

write.csv(
  dominant_asv_summary,
  "Deliverables/ALL/alphadiv.dominant_asv_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# SPLIT DOMINANT-ASV DATA BY SEASON
# ------------------------------------------------------------------------------

summer_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(season == "Summer")

autumn_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(season == "Autumn")

spring_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(season == "Spring")

winter_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(season == "Winter")

# Confirm that dominant-ASV sample tables match their host metadata objects.
setdiff(summer_dominant_asv_per_sample$samp, summer_metadata$samp)
setdiff(autumn_dominant_asv_per_sample$samp, autumn_metadata$samp)
setdiff(spring_dominant_asv_per_sample$samp, spring_metadata$samp)
setdiff(winter_dominant_asv_per_sample$samp, winter_metadata$samp)

# ------------------------------------------------------------------------------
# SEASON-SPECIFIC DOMINANT-ASV SUMMARY TABLES
# ------------------------------------------------------------------------------

summer_dominant_asv_summary <- summer_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(summer_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

autumn_dominant_asv_summary <- autumn_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(autumn_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

spring_dominant_asv_summary <- spring_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(spring_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

winter_dominant_asv_summary <- winter_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(winter_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

summer_dominant_asv_summary
autumn_dominant_asv_summary
spring_dominant_asv_summary
winter_dominant_asv_summary

library(openxlsx)

season_summaries <- list(
  Summer = summer_dominant_asv_summary,
  Autumn = autumn_dominant_asv_summary,
  Spring = spring_dominant_asv_summary,
  Winter = winter_dominant_asv_summary
)

wb <- createWorkbook()

for (season in names(season_summaries)) {
  addWorksheet(wb, season)
  writeData(wb, sheet = season, x = season_summaries[[season]])
  freezePane(wb, sheet = season, firstRow = TRUE)
  setColWidths(wb, sheet = season, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  "Deliverables/ALL/alphadiv/dominant_asv_summary_by_season.xlsx",
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# IDENTIFY TOP DOMINANT ASVs FOR EACH SEASON
# ------------------------------------------------------------------------------

# To prevent unreadable stacked-bar legends, retain the top 5 dominant ASVs
# within each host. All remaining ASVs are grouped into "Other dominant ASVs".
# Change n = 5 if you want more or fewer individually displayed ASVs.

n_top_asvs <- 5

summer_top_dominant_asvs <- summer_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

autumn_top_dominant_asvs <- autumn_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

spring_top_dominant_asvs <- spring_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

winter_top_dominant_asvs <- winter_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

# ------------------------------------------------------------------------------
# FUNCTION: PLOT DOMINANT-ASV IDENTITY WITHIN A HOST
# ------------------------------------------------------------------------------

# The plot displays:
#   x-axis: a grouping variable such as season, Locale, or Age_Group
#   y-axis: proportion of samples in that group
#   fill: ASV that had the highest within-sample read count
#
# position = "fill" makes every bar sum to 100%.

plot_dominant_asv <- function(
    data,
    grouping_var,
    season,
    top_asvs,
    x_label,
    title
) {
  
  grouping_var <- rlang::ensym(grouping_var)
  grouping_var_name <- rlang::as_string(grouping_var)
  
  plot_df <- data %>%
    dplyr::filter(
      !is.na(.data[[grouping_var_name]]),
      trimws(as.character(.data[[grouping_var_name]])) != ""
    ) %>%
    dplyr::mutate(
      dominant_asv_plot = dplyr::if_else(
        taxon_label %in% top_asvs,
        taxon_label,
        "Other dominant ASVs"
      )
    ) %>%
    dplyr::count(
      !!grouping_var,
      dominant_asv_plot,
      name = "n_samples"
    ) %>%
    dplyr::group_by(!!grouping_var) %>%
    dplyr::mutate(
      group_n = sum(n_samples),
      proportion_samples = n_samples / group_n
    ) %>%
    dplyr::ungroup()
  
  ggplot(
    plot_df,
    aes(
      x = !!grouping_var,
      y = proportion_samples,
      fill = dominant_asv_plot
    )
  ) +
    geom_col(
      position = "fill",
      width = 0.75,
      colour = "black",
      linewidth = 0.15
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_fill_manual(
      values = c(
        setNames(
          scales::hue_pal()(length(top_asvs)),
          top_asvs
        ),
        "Other dominant ASVs" = "grey70"
      ),
      drop = FALSE
    ) +
    labs(
      x = x_label,
      y = "Proportion of samples",
      fill = "Dominant ASV",
      title = title
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      legend.position = "right"
    )
}

# ------------------------------------------------------------------------------
# DOMINANT ASV BY SEASON, 
# ------------------------------------------------------------------------------

# SUMMER BY POD
summer_dominant_asv_plot_pod <- plot_dominant_asv(
  data = summer_dominant_asv_per_sample,
  grouping_var = pod_from_id,
  season = "Summer",
  top_asvs = summer_top_dominant_asvs,
  x_label = "Pod",
  title = "Dominant ASVs in summer by pod"
)

summer_dominant_asv_plot_pod

# SEASON BY POD
autumn_dominant_asv_plot_pod <- plot_dominant_asv(
  data = autumn_dominant_asv_per_sample,
  grouping_var = pod_from_id,
  season = "Autumn",
  top_asvs = autumn_top_dominant_asvs,
  x_label = "Pod",
  title = "Dominant ASVs in autumn by pod"
)

autumn_dominant_asv_plot_pod

# SEASON BY POD
spring_dominant_asv_plot_pod <- plot_dominant_asv(
  data = spring_dominant_asv_per_sample,
  grouping_var = pod_from_id,
  season = "Spring",
  top_asvs = spring_top_dominant_asvs,
  x_label = "Pod",
  title = "Dominant ASVs in spring by pod"
)

spring_dominant_asv_plot_pod

# SEASON BY POD
winter_dominant_asv_plot_pod <- plot_dominant_asv(
  data = winter_dominant_asv_per_sample,
  grouping_var = pod_from_id,
  season = "Winter",
  top_asvs = winter_top_dominant_asvs,
  x_label = "Pod",
  title = "Dominant ASVs in winter by pod"
)

winter_dominant_asv_plot_pod


# Combine seasonal dominant-ASV plots
dominant_asv_plots_season <- (
  spring_dominant_asv_plot_pod /
    summer_dominant_asv_plot_pod /
    autumn_dominant_asv_plot_pod / 
    winter_dominant_asv_plot_pod
  
) &
  theme(
    text = element_text(size = 34),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 32),
    plot.title = element_text(size = 35),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )

dominant_asv_plots_season

ggsave(
  "Deliverables/ALL/alphadiv/dominant_asv_by-season_each-pod.png",
  plot = dominant_asv_plots_season,
  width = 18,
  height = 18,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# END DOMINANT-ASV SECTION
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Phylogenetics
# ------------------------------------------------------------------------------

# Faith

### biological interpretation and conclusion

# ------------------------------------------------------------------------------
# Information
# ------------------------------------------------------------------------------
# BY PREDATOR
# ------------------------------------------------------------------------------

# Kruskal-Wallis test
shan_test_season <- kruskal.test(
  Shannon ~ season,
  data = metadata
)

shan_test_season # p = 2.2e-16 (significant)

# Dunn test to determine which host pairs differ
shan_dunn_results_season <- metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(season)
  ) %>%
  rstatix::dunn_test(
    Shannon ~ season,
    p.adjust.method = "hochberg"
  )

shan_dunn_results_season

# Add x- and y-coordinates for comparison brackets
shan_dunn_plot_season <- shan_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

shan_dunn_plot_season

# Box-and-whisker plot by host, with individual samples overlaid
shan_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Shannon, fill = pod_from_id)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    shan_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = pod_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Host",
    y = "Shannon entropy per sample",
    title = "Shannon entropy by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom", 
        # Text sizes in points
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 12)),
        axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 12)),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        
        # Make ticks and plot lines readable on slides
        axis.line = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.25, "cm")
  )

shan_plot_season


ggsave(
  "Deliverables/ALL/alphadiv/SHAN-SEASON_major.png",
  plot = shan_plot_season,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)


# interpretation of all metrics collectively


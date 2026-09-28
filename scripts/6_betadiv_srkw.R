# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS BETA DIVERSITY METRICS ON UBIOME DATA
# THIS IS THE SECOND STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
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
  column_to_rownames(var = "samp")%>%
  as.matrix()

shared_samps <- intersect(metadata$samp, rownames(asv_mat))

metadata <- metadata %>%
  dplyr::filter(samp %in% shared_samps) %>%
  dplyr::arrange(match(samp, shared_samps))

asv_mat <- asv_mat[shared_samps, , drop = FALSE]

stopifnot(identical(metadata$samp, rownames(asv_mat)))

# sets host colors
pod_colors <- c(
  "J" = "#F8766D",  # ggplot2 default red
  "K" = "#00BA38",  # ggplot2 default green
  "L"  = "#619CFF"   # ggplot2 default blue
)


# ------------------------------------------------------------------------------
# Aitchison Distances
# ------------------------------------------------------------------------------
# # STATS : AITCHISON PERMANOVA ON RELATIVE ABUNDANCES 
# # NOT RAREFIED; USED CLR METHOD PER GLOOR ET AL., 2017
# FROM WIKIPEDIA: The null hypothesis states that the centroids (averages) and/or the dispersion (spread) of the groups are equivalent in the multivariate space.
## H0 - there is no difference in the pre samples vs post samples
## HA - there is a difference in the pre samples vs post samples

# Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"
robust_aitch_dist <- vegdist(asv_mat, method = "robust.aitchison")
head(rownames(robust_aitch_dist))

# Run the PERMANOVA-------------------------------------------------------------

# BY SEASON
permanova_aitch_season <- adonis2(robust_aitch_dist ~ season, 
                                data = metadata, 
                                permutations = 999)

print(permanova_aitch_season) # p-value = 0.001 

metadata <- metadata %>%
  dplyr::filter(samp %in% labels(robust_aitch_dist)) %>%
  dplyr::arrange(match(samp, labels(robust_aitch_dist)))

stopifnot(
  identical(metadata$samp, labels(robust_aitch_dist))
)


pairwise_permanova <- function(dist_obj, metadata, group_var,
                               sample_id = "LabID",
                               permutations = 999,
                               p_adjust = "BH") {
  
  # Ensure metadata and distance labels are identically ordered
  metadata <- metadata %>%
    dplyr::filter(.data[[sample_id]] %in% labels(dist_obj)) %>%
    dplyr::arrange(match(.data[[sample_id]], labels(dist_obj)))
  
  stopifnot(identical(metadata[[sample_id]], labels(dist_obj)))
  
  groups <- sort(unique(metadata[[group_var]]))
  pairs <- combn(groups, 2, simplify = FALSE)
  
  purrr::map_dfr(pairs, function(pair) {
    
    meta_sub <- metadata %>%
      dplyr::filter(.data[[group_var]] %in% pair)
    
    dist_sub <- as.dist(
      as.matrix(dist_obj)[meta_sub[[sample_id]], meta_sub[[sample_id]]]
    )
    
    model <- vegan::adonis2(
      dist_sub ~ group,
      data = data.frame(group = meta_sub[[group_var]]),
      permutations = permutations
    )
    
    tibble::tibble(
      group1 = pair[1],
      group2 = pair[2],
      n1 = sum(meta_sub[[group_var]] == pair[1]),
      n2 = sum(meta_sub[[group_var]] == pair[2]),
      pseudo_F = model$F[1],
      R2 = model$R2[1],
      p = model$`Pr(>F)`[1]
    )
  }) %>%
    dplyr::mutate(
      p_BH = p.adjust(p, method = p_adjust),
      significance = dplyr::case_when(
        p_BH < 0.001 ~ "***",
        p_BH < 0.01 ~ "**",
        p_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(p_BH, p)
}


set.seed(123)

pairwise_aitch_season <- pairwise_permanova(
  dist_obj = robust_aitch_dist,
  metadata = metadata,
  group_var = "season",
  sample_id = "samp",
  permutations = 9999,
  p_adjust = "BH"
)

pairwise_aitch_season

# tests for dispersion
disp <- vegan::betadisper(
  robust_aitch_dist,
  group = metadata$season
)

set.seed(123)
vegan::permutest(disp, permutations = 9999) # p-value = 0.1382
## Multivariate dispersion did not differ significantly among host groups
## (PERMDISP: F = 2.-118, p = 0.0517), although the result was borderline.
## Therefore, the significant PERMANOVA result is not accompanied by formal
## evidence of heterogeneous dispersion at alpha = 0.05.

## However, because this result is so borderline, caution:
## this means that significant results in differences may not be due to 
## multivariate microbial compositional differences,
## but may be due to unequal dispersion


# Run the PERMANOVA BY SEASON----------------------------------------------------------------------

permanova_aitch_pod <- adonis2(robust_aitch_dist ~ pod_from_id, 
                                          data = metadata, 
                                          permutations = 999)

print(permanova_aitch_pod) # p-value = 0.028 significant

metadata <- metadata %>%
  dplyr::filter(samp %in% labels(robust_aitch_dist)) %>%
  dplyr::arrange(match(samp, labels(robust_aitch_dist)))

stopifnot(
  identical(metadata$samp, labels(robust_aitch_dist))
)

set.seed(123)

pairwise_aitch_pod <- pairwise_permanova(
  dist_obj = robust_aitch_dist,
  metadata = metadata,
  group_var = "pod_from_id",
  sample_id = "samp",
  permutations = 9999,
  #p_adjust = "BH"
)

pairwise_aitch_pod

# ------------------------------------------------------------------------------
# PCoA ORDIANTION PLOTS
# ------------------------------------------------------------------------------

# Season palette
season_colors <- c(
  "Winter" = "#4E79A7",
  "Spring" = "#59A14F",
  "Summer" = "#E15759",
  "Autumn" = "#B07AA1"
)

# ------------------------------------------------------------------------------
# FUNCTION: CREATE PCoA SCORES FROM A DISTANCE MATRIX
# ------------------------------------------------------------------------------

make_pcoa_df <- function(dist_obj, metadata_df, grouping_var, sample_id = "samp") {
  
  # Make sample ID a plain character vector
  metadata_df <- metadata_df %>%
    dplyr::mutate(
      "{sample_id}" := as.character(.data[[sample_id]])
    ) %>%
    dplyr::filter(.data[[sample_id]] %in% labels(dist_obj)) %>%
    dplyr::arrange(match(.data[[sample_id]], labels(dist_obj)))
  
  # Ensure exact alignment between metadata and distance labels
  stopifnot(
    identical(
      as.character(metadata_df[[sample_id]]),
      labels(dist_obj)
    )
  )
  
  # PCoA / classical multidimensional scaling
  # add = TRUE applies a correction if the distance produces negative eigenvalues.
  pcoa_out <- cmdscale(
    dist_obj,
    k = 2,
    eig = TRUE,
    add = TRUE
  )
  
  # Calculate explained variation from positive eigenvalues only.
  eig_positive <- pcoa_out$eig[pcoa_out$eig > 0]
  
  axis1_pct <- round(
    100 * eig_positive[1] / sum(eig_positive),
    1
  )
  
  axis2_pct <- round(
    100 * eig_positive[2] / sum(eig_positive),
    1
  )
  
  # Construct an ordination data frame and join relevant metadata.
  pcoa_df <- as.data.frame(pcoa_out$points) %>%
    tibble::rownames_to_column(var = sample_id) %>%
    dplyr::rename(
      PCoA1 = V1,
      PCoA2 = V2
    ) %>%
    dplyr::left_join(
      metadata_df,
      by = sample_id
    )
  
  list(
    data = pcoa_df,
    axis1_pct = axis1_pct,
    axis2_pct = axis2_pct,
    grouping_var = grouping_var
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: CREATE A PCoA PLOT
# ------------------------------------------------------------------------------

plot_pcoa <- function(
    pcoa_list,
    grouping_var,
    colors,
    legend_title,
    title,
    subtitle = NULL,
    ellipse = TRUE
) {
  
  grouping_var <- rlang::ensym(grouping_var)
  
  p <- ggplot(
    pcoa_list$data,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = !!grouping_var
    )
  )
  
  # Ellipses require at least three samples per group. They are descriptive only.
  if (ellipse) {
    p <- p +
      stat_ellipse(
        aes(colour = !!grouping_var),
        type = "t",
        level = 0.95,
        linewidth = 0.7,
        show.legend = FALSE
      )
  }
  
  p +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.35,
      size = 5,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    scale_colour_manual(
      values = colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = paste0("PCoA 1 (", pcoa_list$axis1_pct, "%)"),
      y = paste0("PCoA 2 (", pcoa_list$axis2_pct, "%)"),
      fill = legend_title,
      title = title,
      subtitle = subtitle
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

# ==============================================================================
# PCoA BY SEASON
# ==============================================================================

# metadata is already aligned to robust_aitch_dist after your stopifnot above.
# Re-align defensively in case you re-ran earlier chunks out of order.

metadata_season_pcoa <- metadata %>%
  dplyr::filter(samp %in% labels(robust_aitch_dist)) %>%
  dplyr::mutate(
    samp = as.character(samp),
    season = factor(
      season,
      levels = c("Spring", "Summer", "Autumn", "Winter")
    )
  ) %>%
  dplyr::arrange(match(samp, labels(robust_aitch_dist)))

stopifnot(
  identical(
    metadata_season_pcoa$samp,
    labels(robust_aitch_dist)
  )
)

pcoa_season <- make_pcoa_df(
  dist_obj = robust_aitch_dist,
  metadata_df = metadata_season_pcoa,
  grouping_var = "season"
)

p_pcoa_season <- plot_pcoa(
  pcoa_list = pcoa_season,
  grouping_var = season,
  colors = season_colors,
  legend_title = "Season",
  title = "Robust Aitchison PCoA by season",
  subtitle = "PERMANOVA: R² = 0.04563, p = 0.001"
)

p_pcoa_season



# Optional centroid overlay for host plot
season_centroids <- pcoa_season$data %>%
  dplyr::group_by(season) %>%
  dplyr::summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

p_pcoa_season <- p_pcoa_season +
  geom_point(
    data = season_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = season
    ),
    inherit.aes = FALSE,
    shape = 23,
    colour = "black",
    stroke = 0.9,
    size = 7
  )

p_pcoa_season <- p_pcoa_season&
  theme(
    text = element_text(size = 34),
    axis.title = element_text(size = 34),
    axis.text = element_text(size = 32),
    plot.title = element_text(size = 35),
    plot.subtitle = element_text(size = 31),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 31)
  )


ggsave(
  "Deliverables/ALL/betadiv/AITCH_SEASON-major_PCoA.png",
  plot = p_pcoa_season,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)

# ==============================================================================
# PCoA BY POD
# ==============================================================================

metadata_pod_pcoa <- metadata %>%
  dplyr::filter(samp %in% labels(robust_aitch_dist)) %>%
  dplyr::mutate(
    samp = as.character(samp),
    pod_from_id = factor(
      pod_from_id,
      levels = c("J", "K", "L")
    )
  ) %>%
  dplyr::arrange(match(samp, labels(robust_aitch_dist)))

stopifnot(
  identical(
    metadata_pod_pcoa$samp,
    labels(robust_aitch_dist)
  )
)

pcoa_pod <- make_pcoa_df(
  dist_obj = robust_aitch_dist,
  metadata_df = metadata_pod_pcoa,
  grouping_var = "pod_from_id"
)

p_pcoa_pod <- plot_pcoa(
  pcoa_list = pcoa_pod,
  grouping_var = pod_from_id,
  colors = pod_colors,
  legend_title = "Pod",
  title = "Robust Aitchison PCoA by pod",
  subtitle = "PERMANOVA: R² = 0.01186, p = 0.033"
)

p_pcoa_pod

# Optional centroid overlay for host plot
pod_centroids <- pcoa_pod$data %>%
  dplyr::group_by(pod_from_id) %>%
  dplyr::summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

p_pcoa_pod <- p_pcoa_pod +
  geom_point(
    data = pod_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = pod_from_id
    ),
    inherit.aes = FALSE,
    shape = 23,
    colour = "black",
    stroke = 0.9,
    size = 7
  )

p_pcoa_pod <- p_pcoa_pod&
  theme(
    text = element_text(size = 34),
    axis.title = element_text(size = 34),
    axis.text = element_text(size = 32),
    plot.title = element_text(size = 35),
    plot.subtitle = element_text(size = 31),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 31)
  )


ggsave(
  "Deliverables/ALL/betadiv/AITCH_POD-major_PCoA.png",
  plot = p_pcoa_pod,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)


# ==============================================================================
# FINAL COMBINED PCoA FIGURE
# ==============================================================================

all_pcoa_plots <- (
  p_pcoa_pod |
    p_pcoa_season
) &
  theme(
    text = element_text(size = 34),
    axis.title = element_text(size = 34),
    axis.text = element_text(size = 32),
    plot.title = element_text(size = 35),
    plot.subtitle = element_text(size = 31),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 31)
  )

all_pcoa_plots

ggsave(
  "Deliverables/ALL/betadiv/AITCH_major_PCoA.png",
  plot = all_pcoa_plots,
  width = 20,
  height = 24,
  units = "in",
  dpi = 300
)

# ==============================================================================
# END PCoA SECTION
# 
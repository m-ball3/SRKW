# plots

## Packages ----
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggh4x)

## Load phyloseq object for all plates ----
ps.16s     <- readRDS("srkw-ps.16s.allplates")


# names to drop
control.neg <- c("NegCon_S44", "PosCon_S43")

# for 16S
ps.16s <- prune_samples(!(sample_names(ps.16s) %in% control.neg), ps.16s)

## 12S preprocessing for FACET plot ----
ps.16s <- tax_glom(ps.16s, "Species", NArm = FALSE)
ps.16s <- subset_taxa(ps.16s, Class == "Actinopteri")
ps.16s <- prune_samples(sample_sums(ps.16s) >= 100, ps.16s)

ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  x_rel
})

## Family palette for barplot ----
families_16s  <- sort(unique(as.character(tax_table(ps16s.rel)[, "Species"])))
n_fam_16s     <- length(families_16s)
palette_fam_16s <- rep(
  c(brewer.pal(4,"Pastel2"),
    brewer.pal(8,"Paired")),
  length.out = n_fam_16s
)
names(palette_fam_16s) <- families_16s

## Pod color palette (analogous to predator_colors) ----
pod_colors <- c(
  "J"  = "#006D8F",  # deep teal
  "K"  = "#D9A441",  # ochre
  "L"  = "#8B4E80"   # plum
)


## Alpha diversity data (12S + 16S) ----
pr16s_data <- plot_richness(
  ps.16s, x = "Year",
  measures = c("Shannon", "Simpson"),
  color = "Pod"
)$data

pr16s_data$Year <- factor(pr16s_data$Year, levels = sort(unique(pr16s_data$Year)))

dd_box   <- position_dodge2(width = 0.7, preserve = "single")
dd_point <- position_jitterdodge(jitter.width = 0.01, dodge.width = 0.45)

ylim_shannon <- range(c(
  pr16s_data$value[pr16s_data$variable == "Shannon"]
), na.rm = TRUE)


## 12S alpha-diversity plot (Pod colors) ----
div.16 <- ggplot(pr16s_data, aes(x = Year, y = value, color = pod)) +
  geom_boxplot(
    aes(group = interaction(Year, pod), fill = pod),
    position = dd_box,
    outlier.shape = NA,
    width = 0.45,
    color = "grey40",
    alpha = 0.3
  ) +
  geom_jitter(position = dd_point, size = 3, alpha = 0.85)+
  scale_color_manual(values = pod_colors,, na.value = "lightblue") +
  scale_fill_manual(values = pod_colors, , na.value = "lightblue") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  labs(title = "Alpha Diversity", y = "Diversity Index", x = "Year") +
  scale_y_continuous(limits = ylim_shannon, name = "Diversity Index") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(hjust = 1, angle = 45, size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 22),
    strip.text   = element_text(size = 18),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
div.16

## 12S NMDS (Pod colors) ----
ps.prop.16s <- transform_sample_counts(ps.16s, function(otu) otu / sum(otu))
sample_data(ps.prop.16s) <- sample_data(ps.16s)

ord.nmds.bray <- ordinate(ps.prop.16s, method = "NMDS", distance = "bray")

nmds_scores <- as.data.frame(ord.nmds.bray$points)
nmds_scores$SampleID <- rownames(nmds_scores)

metadata <- as.data.frame(sample_data(ps.prop.16s))
metadata$SampleID <- rownames(metadata)

plot_data <- left_join(nmds_scores, metadata, by = "SampleID")

nmds.16 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = pod)) +
  stat_ellipse(
    aes(group = pod, color = pod),  # outline only, by pod
    level = 0.95,
    type  = "t",
    linewidth = 1    # adjust thickness if you like
  ) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = pod_colors,, na.value = "lightblue") +
  labs(
    title = "Bray Curtis",
    x = "NMDS1",
    y = "NMDS2",
    color = "Pod"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16)
  )

nmds.16


## 12S FACET barplot by Pod (family-level) ----
pod_levels <- levels(sample_data(ps16s.rel)$pod)

# make sure this is run once
library(ggh4x)

pod_levels <- levels(sample_data(ps16s.rel)$pod)
faucet.16s <- plot_bar(ps16s.rel, fill = "Species") +
  scale_fill_manual(values = palette_fam_16s) +
  facet_wrap(
    ~ pod,
    nrow   = 1,
    scales = "free_x"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y  = element_text(size = 14, color = "grey60"),
    axis.title.y = element_text(size = 18, color = "grey40"),
    axis.line.y  = element_line(color = "grey80"),
    axis.ticks.y = element_line(color = "grey80"),
    strip.text   = element_text(size = 16, face = "bold", colour = "black"),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.position = "bottom")

faucet.16s


## Combine and save final SRKW plot ----
final_top <- (div.16 + nmds.16) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_srkw <- final_top / faucet.16s

ggsave(
  filename = "./Deliverables/Beautiful Graphics in R/final_patchwork_16sp1_SRWK.png",
  plot     = final_srkw,
  width    = 22,
  height   = 16,
  dpi      = 350,
  units    = "in",
  bg       = "white"
)

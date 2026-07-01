# ------------------------------------------------------------------
# RAREFACTION
# THIS IS THE FIFTH STEP AFTER DADA2
## rownames-match.r, replicates-contaminated.r,
## Phyloseq-16SALL.r, must be run before this!
# ------------------------------------------------------------------

# Sets up environment
library(phyloseq)
library(tidyverse)
library(vegan)

# loads in data
load("srkw-ps.16s.Rdata")

# checks the depths to choose a cutoff
sums <- sample_sums(ps.16s) 
sums[order(sums)]
summary(sample_sums(ps.16s))


??rarefy_even_depth()

# tries a 5,000 read depth cutoff
min_depth <- 25000
ps_filt <- prune_samples(sample_sums(ps.16s) >= min_depth, ps.16s)

summary(sample_sums(ps_filt))


# creates rarefied phyloseq obj
set.seed(1)
ps_rarefied <- rarefy_even_depth(
  ps_filt,
  sample.size = 25000,
  rngseed = 1,
  replace = FALSE,
  verbose = TRUE
)

# gets relative abundance from rarefied ps obj
ps_rarefied.rel <- transform_sample_counts(ps_rarefied, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

save.image(file = "rarefaction.RData")

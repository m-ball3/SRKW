# ------------------------------------------------------------------
# FORMATS SAMDF TO HAVE ROWNAMES THAT MATCH THE ROWNAMES IN SEQTAB.NOCHIM
# THIS IS THE FIRST STEP AFTER DADA2
# ------------------------------------------------------------------

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(tibble)


# ------------------------------------------------------------------
# Loads in Data
# ------------------------------------------------------------------

# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/DADA2/DADA2 Outputs/SRKW-diet-933.Rdata")
out_df<- as.data.frame(out)


mean(out_df$reads.out)
range(out_df$reads.out)

write.csv(out_df, "./Deliverables/ALL/SRKW_ALL-inout.csv", row.names = TRUE)

# Gets sample metadata
samdf <- read.csv("C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/ALL/SRKW_Diet_SDZWA_MetaData.csv")

# Creates pod column from ID column
samdf <- samdf %>%
  mutate(pod = case_when(
    grepl("^J", ID) ~ "J",
    grepl("^K", ID) ~ "K",
    grepl("^L", ID) ~ "L",
    TRUE            ~ NA_character_
  ))

table(samdf$pod)
# J = 188, K = 51,  = 108

table(samdf$ID)

# ------------------------------------------------------------------
# Ensures rownames are the same
# ------------------------------------------------------------------

# Filter duplucate samples not in seqtab.nochim
samdf <- samdf %>%
  filter(!Sample_name %in% c('2015Jan10-02', '2015Jan10-03', '2015Jan30', '2016Feb25-04'))

# Filter pos and neg controls
samdf <- samdf[!grepl("PosCon", samdf$Sample_name), ]
samdf <- samdf[!grepl("Neg", samdf$Sample_name), ]
# Filter pos and neg controls from sequence table (by rownames)
seqtab.nochim <- seqtab.nochim[!grepl("PosCon", rownames(seqtab.nochim)), ]
seqtab.nochim <- seqtab.nochim[!grepl("NegCon", rownames(seqtab.nochim)), ]

# Sets row names to LabID
rownames(samdf) <- samdf$Sample_name


# Removea _S## from seqtab (keep only base name + run number)
# Replaces all dashes with underscores for consistency
rownames(samdf) <- gsub("_", "-", rownames(samdf))
#rownames(samdf) <- gsub(".", "-", rownames(samdf))
rownames(samdf) <- gsub("\\.", "-", rownames(samdf))
rownames(seqtab.nochim) <- gsub("_", "-", rownames(seqtab.nochim))


#further attempts to make them the same
# gsub(" ", "", x) |>            # remove spaces
# # sub("-S[0-9]+$", "", x) |>     # drop trailing -S##


#rownames(samdf) <- sub("-S\\d+$", "", rownames(samdf))

# # Assuming df1 has the row names you want to match against
# # and df2 has the column you want to filter
# samdf <- samdf %>% 
#   filter(Sample_name %in% rownames(seqtab.nochim))

## These are duplicates (from two extractions)
# rownames(seqtab.nochim) <- sub("F22SEP06_01C_2_S45", "F22SEP06_01C_S45", rownames(seqtab.nochim))

# Checks for identical sample rownames in both
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

### ADDS IN THE -S# FROM SEQTAB.NOCHIM TO SAMDF ----------------------------------------------

# function to drop the trailing -S##
strip_S <- function(x) sub("-S[0-9]+$", "", x)

seq_names      <- rownames(seqtab.nochim)
seq_base_names <- strip_S(seq_names)

# in case some base IDs appear multiple times, keep the first
lookup <- tapply(seq_names, seq_base_names, `[`, 1)

# current metadata rownames (base IDs)
meta_names      <- rownames(samdf)
meta_base_names <- meta_names  # assuming these *are* the base IDs already

# find which metadata samples have a matching entry in the lookup
has_match <- meta_base_names %in% names(lookup)

# create new rownames vector
new_meta_names <- meta_names

# for those with a match, replace with the full ID from lookup (with -S##)
new_meta_names[has_match] <- lookup[meta_base_names[has_match]]

# assign back to samdf
rownames(samdf) <- new_meta_names

### -------------------------------------------------------------------------------------


### CREATES NEW SAMDF AND SEQTAB THAT OVERLAPS
# intersection of sample names
common_samps <- intersect(rownames(samdf), rownames(seqtab.nochim))

samdf_filt        <- samdf[common_samps, , drop = FALSE]
seqtab.nochim_filt <- seqtab.nochim[common_samps, , drop = FALSE]


# Sanity check: row names are the same
rownames(samdf_filt)
rownames(seqtab.nochim_filt)

### Save data
save(samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "rownames-match.RData")

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
  mutate(Pod = case_when(
    grepl("^J", ID) ~ "J",
    grepl("^K", ID) ~ "K",
    grepl("^L", ID) ~ "L",
    TRUE            ~ NA_character_
  ))

table(samdf$ID)

# LINDSEY'S EXTRACTIONS  -------------------------------------------------------------------------------------------------------
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

# SOFIA'S EXTRACTIONS  -------------------------------------------------------------------------------------------------------

sofdf <- read.csv("./metadata/Sofia's samples/SRKW_WO_SDZWA_06182024.csv")

sofdf <- sofdf %>%
  dplyr::select(LabID, FieldID, IndID_SDZWA)

sofdf <- sofdf %>%
  rename(Sample_name = LabID)%>%
  rename(ID = IndID_SDZWA)

# rownames
# Reset row names to default numeric indices
rownames(sofdf) <- NULL

# Removes duplicates and NAs in Sample_name
sofdf <- sofdf[!duplicated(sofdf$Sample_name), ]
sofdf <- sofdf[!is.na(sofdf$Sample_name), ]

# makes Sample_name the rownames
rownames(sofdf) <- sofdf$Sample_name
  
### Fixes issue with extra info from Seqtab that is not in sofdf ----------------------------------------------
samples_with_suffix <- c(
  "WADE-002-014-nc-S40",
  "WADE-002-019-nc-S41",
  "WADE-002-024-nc-S39",
  "WADE-002-028-nc-S42",
  "WADE-002-047-rep-S13",
  "WADE-002-054-rep-nc-S43",
  "WADE-002-056-rep-nc-S44",
  "WADE-002-060-rep-S20",
  "WADE-002-063-rep-nc-S45"
)

core_with_S <- sub("^([^-]+-[^-]+-[^-]+).*(-S[0-9]+)$", "\\1\\2", samples_with_suffix)
core_with_S


# find positions of these samples in seqtab
rn <- rownames(seqtab.nochim)

# 
idx <- match(samples_with_suffix, rn)

# sanity check
cbind(old = rn[idx], new = core_with_S)
# check this looks right before assigning

# 3) replace only those rownames
rownames(seqtab.nochim)[idx] <- core_with_S

### ADDS IN THE -S# FROM SEQTAB.NOCHIM TO SOFDF ----------------------------------------------

# current metadata rownames (base IDs)
meta_names.sof      <- rownames(sofdf)
meta_base_names.sof <- meta_names.sof  # assuming these *are* the base IDs already

# find which metadata samples have a matching entry in the lookup
has_match <- meta_base_names.sof %in% names(lookup)

# create new rownames vector
new_meta_names.sof <- meta_names.sof

# for those with a match, replace with the full ID from lookup (with -S##)
new_meta_names.sof[has_match] <- lookup[meta_base_names.sof[has_match]]

# assign back to samdf
rownames(sofdf) <- new_meta_names.sof



# AMY AND KIM'S EXTRACTIONS  -------------------------------------------------------------------------------------------------

# loads in metadata csv
akdf <- read.csv("./PRJNA1068648/new_prey_meta_10.26.23_pedPod.csv")

akdf <- akdf %>%
  rename(Sample_name = Sample) %>%
  filter(Population == "SRKW")

# loads in csv with sample name and SRA accession #
SRA_list <- read.csv("./PRJNA1068648/srkw_sra_list.csv")

SRA_list <- SRA_list %>%
  rename(Sample_name = SampleName) %>%
  dplyr::select(Sample_name, Run)

# Left joins by Sample_name
akdf <- akdf %>%
  left_join(SRA_list, by = "Sample_name")

# renames columns in akdf to match those in samdf and sofdf
akdf <- akdf %>%
  rename(ID = Individual.ID)%>%
  rename(Sex = Genetic.Sex) %>%
  rename(Year = year) %>%
  rename(Month = month) %>%
  rename(Day = day)

# Changes numeric month to char month abbrev.
akdf$Month <- factor(month.abb[akdf$Month],
                     levels = month.abb)  # keeps calendar order

# cleans NAs
akdf <- akdf[!is.na(akdf$Run), ]

# makes Run the rownames
rownames(akdf) <- akdf$Run

### REMOVES THE -R1 FROM SEQTAB.NOCHIM TO MATCH RUN IN AKDF ----------------------------------------------
rownames(seqtab.nochim) <- sub("-R1$", "", rownames(seqtab.nochim))

# ------------------------------------------------------------------
# Binds all three metadata dfs into one master df
# ------------------------------------------------------------------
samdf_all <- dplyr::bind_rows(samdf, sofdf, akdf)

#---------------------------------------------------------------------------------------------------------

### CREATES NEW SAMDF AND SEQTAB THAT OVERLAPS-------------------------------------------------------------------------------
# intersection of sample names
common_samps <- intersect(rownames(samdf_all), rownames(seqtab.nochim))


# NEED TO FIX THESE:  [7] "WADE-002-024-S39"   "WADE-002-028-S42"   "WADE-002-047-S13"   "WADE-002-054-S43"   "WADE-002-056-S44"   "WADE-002-060-S20"  
# [13] "WADE-002-063-S45" 
# 
# I THINK THEY ARE REPLICATES, NEGATIVE CONTROLS, ETC.

# Samples in metadata but not in OTU table
setdiff(rownames(samdf_all), rownames(seqtab.nochim))

# Samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf_all))


samdf_filt        <- samdf_all[common_samps, , drop = FALSE]
seqtab.nochim_filt <- seqtab.nochim[common_samps, , drop = FALSE]


# Sanity check: row names are the same
rownames(samdf_filt)
rownames(seqtab.nochim_filt)



table(samdf$pod)
# J = 188, K = 51,  = 108

### Save data
save(samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "rownames-match.RData")

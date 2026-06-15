# Loads in output from DADA2 + filtered seqtab.nochim and samdf from rownames-match.r
load("rownames-match.RData")


# Removes Sofia's samples that have great white shark contamination
GWS_contaminated <- c(
  "WADE-002-002-S25",
  "WADE-002-004-S11",
  "WADE-002-005-S32",
  "WADE-002-010-S1",
  "WADE-002-011-S2",
  "WADE-002-014-S27"
)

samdf_filt <- samdf_filt[!rownames(samdf_filt) %in% GWS_contaminated, ]

seqtab.nochim_filt <- seqtab.nochim_filt[!rownames(seqtab.nochim_filt) %in% GWS_contaminated, ]


# REMOVES TISSUE SAMPLE (NOT AN SRKW FECAL SAMPLE)
# DEALS WITH REPLICATES
## Visually compare their absolute and proportional abundance in a separate bar plot. 
## You can also run NMDS on them to see if they cluster. 
## This will tell us a little bit about sequencing/subsampling bias on our data.


save(samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "replicates-contaminated.RData")

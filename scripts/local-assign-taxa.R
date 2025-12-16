# ------------------------------------------------------------------
# ASSIGN TAXONOMY, ASSIGN SPECIES, ADD SPECIES
# more flexibility to easily modify species assignments 
# by doing this step of DADA2 locally
# ------------------------------------------------------------------

# loads in environment

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(tidyverse)

# lOADS IN TAXONOMY REFERENCE FILES

# 12s
taxref <- "./DADA2/Ref-DB/16S_salmon_groundfish_reference_database_2022.fasta"

# lOADS IN RDATA FILE FROM DADA2

# Plate1
load("DADA2/DADA2 Outputs/Plate1/SRKW-diet-16SP1.Rdata")

# # Plate3
# load("DADA2/DADA2 Outputs/Plate3/SRKW-diet-16SP3.Rdata")

# ASSIGNS

### Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95)

# SAVE DATA

# OVERWRITES P1 DADA2 OUTPUT
save(seqtab.nochim, freq.nochim, track, taxa, tax_table, file = "DADA2/DADA2 Outputs/Plate1/SRKW-diet-16SP1.Rdata")

# OVERWRITES P3 DADA2 OUTPUT
# save(seqtab.nochim, freq.nochim, track, taxa, tax_table, file = "DADA2/DADA2 Outputs/Plate3/SRKW-diet-16SP3.Rdata")
# 


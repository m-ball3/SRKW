# save this as test_apptainer.R

# ## dada2 QAQC of Arctic Predator fastq sequences
# ## 1/3/2022 - modified 5/27/25 by mball for WADE-003 Arctic Predator Diets
# ## Amy Van Cise
#
# Set library path to /gscratch (replace with your full path)
.libPaths("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs")

# Create the directory if it doesn't exist
dir.create("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs", recursive = TRUE, showWarnings = FALSE)

# Verify the path is set correctly
.libPaths()

### set up working environment

devtools::install_github("benjjneb/dada2", ref="v1.16", lib = .libPaths()[1])
BiocManager::install("dada2", lib = .libPaths()[1], force = TRUE)
BiocManager::install("S4Vectors")
install.packages("seqinr", lib = .libPaths()[1])

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(dplyr)
# ------------------------------------------------------------------
# THIS IS THE FIFTH STEP AFTER DADA2
## rownames-match.r, replicates-contaminated.r,
## Phyloseq-16SALL.r, must be run before this!
# ------------------------------------------------------------------

##THIS CODE IS HERE TO INVESTIGATE ALL, MAJOR, AND MINOR PREY ITEMS BY POD
# GOAL: LARGE .XLXS WITH TABS FOR THE DIFFERENT PODS

## AFTER THESE CHECKS, MAKE A LIST OF MINOR (6) PREY ITEMS THAT ARE NOT IN THE MAJOR PREY ITEMS LIST
## DROP NAS

# how many samples? avg. proportion (e.g. is 2% in 2 samples, its secondary prey)

library(tidyverse)


# loads in data
load("SRKW-diet-MASTER.Rdata")

otu.prop.major <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv")
otu.prop.minor <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MINOR.csv")
otu.prop.all <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-ALL.csv")

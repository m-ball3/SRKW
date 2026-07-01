# ------------------------------------------------------------------
# THIS IS THE FOURTH STEP AFTER DADA2
## rownames-match.r, replicates-contaminated.r,
## Phyloseq-16SALL.r, must be run before this!
# ------------------------------------------------------------------

##THIS CODE IS HERE TO CHECK ON ANY 'WEIRD' SPECIES ASSIGNMENTS WE GOT FOR SRKW DIET
 # CITHARICHTHYS IS WEIRD
    ## WHAT IS THE PROPORTION OF C. ACROSS ALL SAMPLES?
    ## HOW MANY SAMPLES IS IT IN?
 # SAME FOR SPOTTED RATFISH (COLLEI)
 # AND FOR GADUS, CLUPEA, AND AMODYTES - THESE MAY BE FINE, BUT WE STILL WANT TO DOUBLE CHECK
    ## CHECK AMY'S PAPER AND SEE IF SHE DETECTED ANY OF THESE SPECIES

## AFTER THESE CHECKS, MAKE A LIST OF MINOR (6) PREY ITEMS THAT ARE NOT IN THE MAJOR PREY ITEMS LIST
    ## DROP NAS

# how many samples? avg. proportion (e.g. is 2% in 2 samples, its secondary prey)

library(tidyverse)


# loads in data
otu.prop <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv")

# creates new df with OTU assignments as rows
  out <- otu.prop %>%
  pivot_longer(
    cols = -X,
    names_to  = "species",
    values_to = "abund"
  ) %>%
  group_by(species) %>%
  summarise(
    abund = sum(abund, na.rm = TRUE) / n(),   # n() = number of rows in this species group
    .groups = "drop"
  )

  out <- out %>%
    arrange(desc(abund))

  writexl::write_xlsx(out, "MAJOR-species-assignment-checks.xlsx")
  
  getwd()
  
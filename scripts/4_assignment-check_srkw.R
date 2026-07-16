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
load("SRKW-diet-MASTER.Rdata")
otu.prop.major <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv")
otu.prop.minor <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MINOR.csv")
otu.prop.all <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-ALL.csv")


# checks out the sample data
sampdf<- as.data.frame(ps16s.rel@sam_data)

# creates a vector that gives the common name for each species
# Scientific → common name lookup
common_names <- c(
  "Oncorhynchus.tshawytscha"         = "Chinook salmon",
  "Oncorhynchus.keta"                = "Chum salmon",
  "Oncorhynchus.kisutch"             = "Coho salmon",
  "Ophiodon.elongatus"               = "Lingcod",
  "Oncorhynchus.nerka"               = "Sockeye salmon",
  "Hippoglossus.stenolepis"          = "Pacific halibut",
  "Oncorhynchus.mykiss"              = "Steelhead trout",
  "Anoplopoma.fimbria"               = "Sablefish (black cod)",
  "Raja.binoculata"                  = "Big skate",
  "Atheresthes.stomias"              = "Arrowtooth flounder",
  "Hydrolagus.colliei"               = "Spotted ratfish",
  "Clupea.pallasii"                  = "Pacific herring",
  "Gadus.chalcogrammus"              = "Alaska pollock",
  "Citharichthys.sordidus"           = "Pacific sanddab",
  "Ammodytes.hexapterus"             = "Pacific sand lance", 
  "Citharichthys.stigmaeus"          = "Speckled sanddab",
  "Raja.rhina"                       = "Longnose skate", 
  "Oncorhynchus.gorbuscha"           = "Pink Salmon",
  "Cymatogaster.aggregata"           = "Shiner perch",
  "Myoxocephalus.polyacanthocephalus"= "Great sculpin",
  "Glyptocephalus.zachirus"          = "Rex sole", 
  "Microstomus.pacificus"            = "Pacific Dover Sole",
  "Liparis.pulchellus"               = "Showy snailfish", 
  "Salmo.salar"                      = "Atlantic salmon",
  "Anarrhichthys.ocellatus"          = "Wolf eel"
  
  # NA / NA.1 will stay NA for common_name
)

#MAJOR PREY ITEMS------------------------------------------------------------
# creates new df with OTU assignments as rows
  out.major <- otu.prop.major %>%
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

  out.major <- out.major %>%
    arrange(desc(abund))

  
  # adds common name column to otu table
  # Add common_name column using the lookup vector
  out.major <- out.major %>%
    mutate(common_name = common_names[species])

  writexl::write_xlsx(out.major, "./Deliverables/ALL/MAJOR-species-assignment-checks.xlsx")
  
  getwd()
  
#MINOR PREY ITEMS------------------------------------------------------------
  # creates new df with OTU assignments as rows
  out.minor <- otu.prop.minor %>%
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
  
  out.minor <- out.minor %>%
    arrange(desc(abund))
  
  
  # adds common name column to otu table
  # Add common_name column using the lookup vector
  out.minor <- out.minor %>%
    mutate(common_name = common_names[species])
  
  writexl::write_xlsx(out.minor, "./Deliverables/ALL/MINOR-species-assignment-checks.xlsx")
  
#ALL PREY ITEMS---------------------------------------------------------------------
  # creates new df with OTU assignments as rows
  out.all <- otu.prop.all %>%
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
  
  out.all <- out.all %>%
    arrange(desc(abund))
  
  
  # adds common name column to otu table
  # Add common_name column using the lookup vector
  out.all <- out.all %>%
    mutate(common_name = common_names[species])
  
  writexl::write_xlsx(out.all, "./Deliverables/ALL/ALL-species-assignment-checks.xlsx")
  
  
## BY POD -----------------------------------------------------------------------------------
rownames(otu.prop.major)
  
  
# 2) Build a long data frame: sample, species, abundance
  df_major <- otu.prop.major %>%
    as.data.frame() %>%
    dplyr::rename(Sample_names = X)%>%
    pivot_longer(
      cols = -Sample_names,
      names_to = "species",
      values_to = "abund"
    )
    
  # 3) Add Pod metadata
  
sam_major <- sam_major %>%
  rownames_to_column()%>%
  dplyr::rename(Sample_names = rowname)

  df_major <- df_major %>%
    left_join(
      sam_major %>% select(rowname, pod_from_id),
      by = "Sample_names"
    )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### CODE FOR SAVING AS ONE BIG EXCEL FILE WITH TABS
  # 1. Install and load required packages
  install.packages(c("readxl", "writexl"))
  library(readxl)
  library(writexl)
  
  # 2. Get paths for all Excel files in your target directory
  file_paths <- list.files(path = "path/to/your/files", pattern = "\\.xlsx$", full.names = TRUE)
  
  # 3. Read files into a named list (Tabs will use the filenames)
  file_names <- tools::file_path_sans_ext(basename(file_paths))
  excel_list <- lapply(file_paths, read_excel)
  names(excel_list) <- file_names
  
  # 4. Save all datasets as tabs in one Excel file
  write_xlsx(excel_list, path = "combined_output.xlsx")
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
  
# list of ARKW in 2024 study
## to remove from HYAK
library(tidyverse)


sra    <- read.csv("PRJNA1068648/sra_metadata_with_paths.csv",           
                   stringsAsFactors = FALSE)
meta   <- read.csv("PRJNA1068648/new_prey_meta_10.26.23_pedPod.csv",   
                   stringsAsFactors = FALSE)

meta_sub <- meta %>% 
  select(Sample, Population)

sra_sub  <- sra  %>% 
  select(Run, SampleName)

merged <- sra_sub %>%
  left_join(meta_sub, by = c("SampleName" = "Sample"))

srkw <- merged %>%
  filter(Population == "SRKW")

arkw <- merged %>%
  filter(Population == "ARKW")

write.csv(srkw, "srkw_sra_list.csv", row.names = FALSE)
write.csv(arkw, "arkws_sra_list.csv", row.names = FALSE)
cat("Total ARKW SRA runs found:", nrow(arkw), "\n")

arkw_runs <- arkw |>           # from the previous script
  dplyr::select(Run)

write.csv(arkw_runs, "arkw_runs_only.txt", row.names = FALSE, quote = FALSE)

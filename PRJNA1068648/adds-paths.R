library(tidyverse)
library(dplyr)

sra_final_complete <- read_csv("PRJNA1068648/SraRunInfo.csv")

# use your final metadata object
runs <- sra_final_complete %>%
  distinct(Run) %>%
  arrange(Run)

write.table(runs$Run,
            file = "run_ids.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Also save metadata with paths for later
sra_final_complete %>%
  mutate(
    fastq_R1 = paste0("fastq/", Run, "_1.fastq.gz"),
    fastq_R2 = paste0("fastq/", Run, "_2.fastq.gz")
  ) %>%
  write_csv("sra_metadata_with_paths.csv")

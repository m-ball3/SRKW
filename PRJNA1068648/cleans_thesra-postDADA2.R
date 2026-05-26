# Merges the 
library(tidyverse)

# Amy and Kim's project 
meta   <- read.csv("PRJNA1068648/new_prey_meta_10.26.23_pedPod.csv",   
                   stringsAsFactors = FALSE)

# Gets sample metadata
samdf <- read.csv("C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/ALL/SRKW_Diet_ALL_MetaData.csv")

# Selects the correct columns
meta_sub <- meta %>% 
  select(Sample, Individual.ID, Pod, Genetic.Sex, year, month, day)

sra_sub  <- samdf  %>% 
  select(notes, Sample_name)

# merges the columns of 
merged <- sra_sub %>%
  left_join(meta_sub, by = c("Sample_name" = "Sample"))

# grabs just the ones with SRR in the notes
# base R
df_srr <- merged[grepl("SRR", merged$notes), ]

# dplyr
library(dplyr)
df_srr <- merged %>% 
  filter(grepl("SRR", notes))

df_srr <- df_srr %>%
  dplyr::rename(Year = year)%>%
  dplyr::rename(Month = month)%>%
  dplyr::rename(Day = day)


write.csv(df_srr, "df_srr.csv", row.names = FALSE)

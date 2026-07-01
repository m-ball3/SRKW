# Combines Plate 1-5 metatada into one daasheet

 # Sets up Environment
library(tidyverse)
library(readr)

# Loads in the metadata for each plate
plate_paths <- c(
  "C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate1/SRKW_Diet_Plate1_MetaData.csv",
  "C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate2/SRKW_Diet_Plate2_MetaData.csv",
  "C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate3/SRKW_Diet_Plate3_MetaData.csv",
  "C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate4/SRKW_Diet_Plate4_MetaData.csv",
  "C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate5/SRKW_Diet_Plate5_MetaData.csv"
)


# Counts the rows in each and then prints what the appended df nrow should be
plate_counts <- plate_paths |>
  set_names() |>
  map_df(\(.f) {
    tibble(
      plate_file = .f,
      n_rows = nrow(read_csv(.f))
    )
  })

print(plate_counts)

expected_total <- sum(plate_counts$n_rows)
cat("Expected total rows:", expected_total, "\n")

# Appends all together
ALL <- plate_paths |>
  set_names() |>
  map_dfr(read_csv, .id = "plate_file")

getwd()
write_csv(ALL, "./metadata/ALL/SRKW_Diet_SDZWA_MetaData.csv")


"C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate1/SRKW_Diet_Plate1_MetaData.csv"
"C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate2/SRKW_Diet_Plate2_MetaData.csv"
"C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate3/SRKW_Diet_Plate3_MetaData.csv"
"C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate4/SRKW_Diet_Plate4_MetaData.csv"
"C:/Users/MBall/OneDrive - UW/Documents/WADE LAB/SRKW/metadata/Plate5/SRKW_Diet_Plate5_MetaData.csv"



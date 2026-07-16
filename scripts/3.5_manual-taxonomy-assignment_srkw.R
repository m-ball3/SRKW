# ------------------------------------------------------------------
# MANUALLY ASSIGNS TAXONOMY WHEN NEEDED
# THIS IS AN OPTIONAL INTERMEDIATE STEP BEFORE TAX GLOM IN STEP THREE
## rownames-match.r and replicates-contaminated.r and 
## up to line 63 in 3_phyloseq_srkw.R must be run before this!
# ------------------------------------------------------------------

#loads in data
asv_table <- read.csv("./Deliverables/ALL/asv_table.csv")

asv_nas <- asv_table %>%
  filter(is.na(Species))

na_list <- c(
  ASV32 = "Onchorhynchus spp.", # best match to mykiss
  ASV36 = "Microstomus pacificus", 
  ASV432 = "Thunnus spp.", # best match (+geographic info) to albacares
  #ASV713 = , # best match to salmo salar (>98%), but also matches to chinook (97%)
  ASV1098 = "Gadus chalcogrammus", 
 # ASV1485 = , # best matches to smelt (likely secondary prey item)
  ASV1684 = "Cryptacanthodes spp.", # best match to Cryptacanthodes bergi
  ASV2469 = "Sebastes diploproa", 
  ASV2689 = "Hydrolagus colliei"
)

# assign back in the asv_table
asv_table$Species_fixed <- asv_table$Species
idx <- match(names(na_list), asv_table$ASV)
asv_table$Species_fixed[idx] <- na_list         

write.csv(asv_table, "./Deliverables/ALL/asv_table-fixed.csv", row.names = TRUE)

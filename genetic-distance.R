`# 16S reference database check

AVC December 2025
library(tidyverse)
library(Biostrings)
library(pwalign)
library(DECIPHER)
library(ape)
library(phangorn)

# reads in the .fasta file and tells it that it is a fasta
# saves with each string (?) as fasta header
Get data -------------------------------------------------------------------
  seqs.16s <- readDNAStringSet("C:/Users/Amy Van Cise/Downloads/16S_salmon_groundfish_reference_database_2022.fasta")

check for Ns
has_N <- which(vcountPattern("N", seqs.16s) > 0)
#####19 sequences with Ns

n_Ns <- vcountPattern("N", seqs.16s)[has_N]
summary(n_Ns)

H1: duplicated sequences in the reference database -------------------------
  dup_idx <- duplicated(as.character(seqs.16s))

sum(dup_idx)
#####102 duplicate sequences

duplicated.seqs <- seqs.16s[dup_idx]
dupseq.names <- names(seqs.16s)[dup_idx]
dupseq.row <- which(duplicated(seqs.16s))
alldup.row <- which(duplicated(seqs.16s) | rev(duplicated(rev(seqs.16s))))

check whether duplicate sequences have the same fasta headers
dup_df <- data.frame(row = alldup.row,
                     name = names(seqs.16s)[alldup.row],
                     sequence = as.character(seqs.16s)[alldup.row],
                     stringsAsFactors = FALSE) %>%
  group_by(sequence) %>%
  mutate(group = cur_group_id()) %>%
  relocate(group, .before = name) %>%
  mutate(fasta.row = row*2, .after = row)

name_counts <- dup_df %>%
  group_by(sequence) %>%
  distinct(name, .keep_all = TRUE) %>%
  summarize(nName = n())

identify the sequences with multiple headers
multi_name_seqs <- name_counts %>%
  filter(nName > 1) %>%
  distinct(sequence) %>%
  pull(sequence)
#####32 sequences with more than one fasta header

multi_name_rows <- dup_df %>%
  filter(sequence %in% multi_name_seqs)

H2: salmon sequences are very similar to each other
get Onc seqs only
seqs_onc <- seqs.16s[grepl("Oncorhynchus", names(seqs.16s), ignore.case = TRUE)]

calculate distance
dmat <- as.matrix(pwalign::stringDist(seqs_onc, method = "levenshtein"))

get pairs with 3 or fewer matches
sim_pairs <- which(dmat <= 3 & dmat > 0, arr.ind = TRUE)
exact_pairs <- which(dmat == 0 & row(dmat) != col(dmat), arr.ind = TRUE) #####confirm no exact matches

simpairs_df <- data.frame(row1 = sim_pairs[,1],
                          row2 = sim_pairs[,2],
                          name1 = names(seqs_onc)[sim_pairs[,1]],
                          name2 = names(seqs_onc)[sim_pairs[,2]],
                          mismatches = dmat[sim_pairs]) %>%
  filter(row1 < row2) %>%
  separate(name1, into = c(NA,NA,NA,NA,NA,NA,"species1"), sep = ";") %>%
  separate(name2, into = c(NA,NA,NA,NA,NA,NA,"species2"), sep = ";")

all similar sequences are within the same species, yay.
H3: salmon sequences are not monophyletic
remove --- from sequences
seqs_onc_clean <- DNAStringSet(gsub("-", "", as.character(seqs_onc)))
row_numbers <- match(as.character(seqs_onc), as.character(seqs.16s))

species_names <- as.data.frame(as.character(names(seqs_onc_clean))) %>%
  separate(1, into = c(NA,NA,NA,NA,NA,NA,"species"), sep = ";") %>%
  pull(species)
names(seqs_onc_clean) <- species_names

align sequences
alignment <- AlignSeqs(seqs_onc_clean)

calculate distance (# bp different)
  dist_obj <- as.dist(DistanceMatrix(alignment, type="matrix"))
  
  build a tree
  tree <- nj(dist_obj)
  tree$tip.label <- paste0(species_names, "_row", row_numbers)
  plot(tree, cex = 0.7)
  
  midpoint tree
  tree_mid <- midpoint(tree)
  tree_mid$tip.label <- paste0(species_names, "_row", row_numbers)
  plot(tree_mid, cex = 0.7)
  
  save pdfs (to zoom in)
  pdf(file = "C:/Users/Amy Van Cise/Downloads/Rplot.pdf", height = 15, width = 40)
  plot(tree_mid, cex = 0.7)
  dev.off()
  
  pdf(file = "C:/Users/Amy Van Cise/Downloads/Rplot2.pdf", height = 15, width = 40)
  plot(tree, cex = 0.7)
  dev.off()`
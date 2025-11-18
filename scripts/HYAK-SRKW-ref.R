# ------------------------------------------------------------------
# Uploads .fastq files to HYAK (need to be in the folder containing the files)
# ------------------------------------------------------------------

scp -r * mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/SRKW/rawdata/16SP1/

# ------------------------------------------------------------------
# Uploads reference Db to HYAK
# ------------------------------------------------------------------

  # 16S
  scp -r "C:/Users/Intern/SRKW/DADA2/Ref-DB//.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/SRKW/
  scp -r "C:/Users/Intern/SRKW/DADA2/Ref-DB/SRKW-16S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/SRKW/
  
  # scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S_Arctic_predator_reference_database_07_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  # scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  # 
  # ------------------------------------------------------------------
# Gets output from HYAK
# ------------------------------------------------------------------

#16S
scp -r "mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/SRKW-diet-16SP1.Rdata" "C:/Users/Intern/SRKW/DADA2/DADA2 Outputs"



  
  
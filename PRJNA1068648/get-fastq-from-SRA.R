# 1. navigate to the desired folder

# Run
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz


# Then Run
# Check what's actually there
ls sratoolkit.3.3.0-ubuntu64/
  
  # Should see 'bin/' directory
  ls sratoolkit.3.3.0-ubuntu64/bin/
  
  # Now add to PATH (easiest - no copying needed)
  export PATH="/gscratch/coenv/mball3/530-meta-analysis/sratoolkit.3.3.0-ubuntu64/bin:$PATH"

# Test it works
prefetch --version
fasterq-dump --version

echo 'export PATH="/gscratch/coenv/mball3/SRKW/2024/sratoolkit.3.3.0-ubuntu64/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Then Run
cd /gscratch/coenv/mball3/SRKW/2024

# Make directories
mkdir -p sra fastq logs tmp

# Configure SRA cache (one-time)
vdb-config --interactive
# Set cache to: /gscratch/coenv/mball3/SRKW/2024/sra

# Test 1 run first
head -1 run_ids.txt | prefetch -O sra
ls sra/
  
  # Then
  prefetch --version

#Then
export PATH="/gscratch/coenv/mball3/SRKW/2024/sratoolkit.3.3.0-ubuntu64/bin:$PATH"

#Then
# 1. Make directories
mkdir -p sra fastq logs tmp

# 2. Test download 1 run
head -1 run_ids.txt | prefetch -O sra

# 3. Check it worked
ls sra/
  
  
  # Then
  # Clean run_ids.txt (remove headers/quotes)
  sed '1d' run_ids.txt | sed 's/"//g' > run_ids_clean.txt
mv run_ids_clean.txt run_ids.txt
head run_ids.txt


#Then CONFIGURE THE INTERACTIVE THING WITH VDB??

#then
cd /gscratch/coenv/mball3/SRKW/2024
mkdir -p fastq sra logs tmp

cd /gscratch/coenv/mball3/SRKW/2024

# Remove header and quotes if present, and strip Windows CR
sed '1d' run_ids.txt | sed 's/"//g' | sed 's/\r$//' > run_ids_clean.txt
mv run_ids_clean.txt run_ids.txt

head run_ids.txt   # should just be plain SRR... IDs, one per line


##THEN
cd /gscratch/coenv/mball3/SRKW/2024

# Make sure sratoolkit is on your PATH in this shell
export PATH="/gscratch/coenv/mball3/SRKW/2024/sratoolkit.3.3.0-ubuntu64/bin:$PATH"
prefetch --version

# Background download of ALL SRR IDs into ./sra
prefetch --option-file run_ids.txt -O sra > logs/prefetch.log 2>&1 &
  tail -f logs/prefetch.log

##AND, FINALLY, AFTER THAT
cd /gscratch/coenv/mball3/SRKW/2024

while read -r run_id; do
echo "=== fasterq: $run_id ==="
/gscratch/coenv/mball3/SRKW/2024/sratoolkit.3.4.1-ubuntu64/bin/fasterq-dump \
"sra/${run_id}" \
--outdir fastq \
--split-files \
-e 8 \
-t tmp
done < run_ids.txt


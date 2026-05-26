#!/bin/bash
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2          # Partition/queue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=60G                   # Add units (e.g. G)
#SBATCH --time=48:00:00             # hh:mm:ss

#SBATCH --job-name=SRKW-ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mball3@uw.edu
#SBATCH --output=runRscript_%j.log  # %j = job ID
#SBATCH --error=runRscript_error_%j.log
#SBATCH --chdir=/mmfs1/home/mball3/
#SBATCH --export=ALL

module load apptainer

apptainer run \
/gscratch/coenv/containerstidyverse_4.0.1.sif \
Rscript /mmfs1/home/mball3/SRKW-dada2.R
#!/bin/bash
#SBATCH --job-name=download
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=download.out
#SBTACH --error=download.err
# Title: 
# Author: Amanda Zacharias
# Date: 2024-10-21
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
#
# Code -------------------------------------------
echo Job started at $(date +%T)
# Dependencies
module load StdEnv/2023

# Download data for 1_singleSet and 2_multipleSets tutorials
wget https://figshare.com/ndownloader/files/42997198

# Download data for 3_multipleSetsDiffCells
wget https://figshare.com/ndownloader/files/38830080
wget https://figshare.com/ndownloader/files/38830077

# Fix names of files
mv 42997198 data_humanSkin_CellChat.rda
mv 38830080 cellchat_embryonic_E13.rds
mv 38830077 cellchat_embryonic_E14.rds

echo Job ended at $(date +%T)
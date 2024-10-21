#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Prepare data for CellChat analysis
# Author: Amanda Zacharias
# Date: 2024-10-21
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# Threading and the presto R package help with compute efficiency.

# Options -----------------------------------------------
options(stringsAsFactors = FALSE)

# Packages -----------------------------------------------
library(CellChat) # 2.1.2

# Source -----------------------------------------------

# Pathways -----------------------------------------------
## Input ===========
dataDir <- "0_data"

## Output ===========
baseDir <- "1_singleSet"
rDataDir <- file.path(baseDir, "rDataDir")

dir.create(rDataDir, showWarnings = FALSE)

# Load data -----------------------------------------------
# Rows = genes, Columns = cells
load(file.path(dataDir, "data_humanSkin_CellChat.rda"))
# data_humanSkin is a combined dataset w/ normal and diseased samples.
data.input <- data_humanSkin$data # normalized gene count matrix
meta <- data_humanSkin$meta # cell metadata
cell.use <- rownames(meta)[meta$condition == "LS"] # cell names of disease data

# Subset for cellchat -----------------------------------------------
# We only want diseased cells for this tutorial.
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$labels)
# 1] Inflam. FIB  FBN1+ FIB    APOE+ FIB    COL11A1+ FIB cDC2
#  [6] LC           Inflam. DC   cDC1         CD40LG+ TC   Inflam. TC
# [11] TC           NKT
# 12 Levels: APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 ... NKT

# Create cellchat object -----------------------------------------------
cellchat <- createCellChat(object = data.input,  meta = meta, group.by = "labels")
# `addMeta()` could be used to add metadata after the fact

## Number of cells in each cell group
groupSize <- as.numeric(table(cellchat@idents))
# > groupSize
#  [1] 1228  813  181  484  121  294   67   81  765  266  630   81

# Set ligand-receptor interaction database -----------------------------------------------
# CellChatDB is manually curated db for humans and mice (& zebrafish available).
# CellChatDB2 contains molecular interactions, autocrine/signaling signaling,
#   ecm-receptor, cell-cell contact, and non-protein signaling. More interactions than v1.
#   ALso extra annotations of l-r pairs s.a. UniProtKB keyworks.
# Can add your own curated l-r pairs.
# Non-protein signaling pairs aren't used by default.
CellChatDB <- CellChatDB.human
# Explore the DB
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

## Subset to only secreted signaling pairs
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
# Other subsetting options available...

## Set the used database in the cellchat object
cellchat@DB <- CellChatDB.use

# Preprocess expression data -----------------------------------------------
# Cellchat compares over-expressed ligands in one group with over-expressed receptors in another.
# Interactions are present if either ligand or receptor is over-expressed.
# Can also project gene expresion data to a protein-protein interaction network.
# CellChat uses raw data by default, not projected data.
## Subset data by signaling genes to save compute power ===========
cellchat <- subsetData(cellchat) # necessary even if you are using the entire DB
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 693

## Optional PPI projection ===========
# cellchat <- projectData(cellchat, PPI.human) # nolint
# when running `computeCommunProb()`, set `raw.use = FALSE` so projected values are used.

# Save -----------------------------------------------
saveRDS(cellchat, file.path(rDataDir, "cellchat_part1.rds"))

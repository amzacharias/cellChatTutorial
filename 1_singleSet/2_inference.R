#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Part 2: Infernce of cell-cell communication network
# Author: Amanda Zacharias
# Date: 2024-10-21
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0

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
plotsDir <- file.path(baseDir, "plots")

dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
cellchat <- readRDS(file.path(rDataDir, "cellchat_part1.rds"))

# Compute communication probablility -----------------------------------------------
# Cellchat assigns each interaction w/ a probability value and permutation test.
#   Integrates gene expression with prior known knowledge between ligands, receptors,
#   and cofactors using the law of mass action.
# Law of mass action: rate of chemical reaction is directly proportional to the product
#   of the activities/concentrations of the reactants.
# The # of inferred paires depends on the method for calculating average gene expression
#   per cell group. Cellchat uses "trimean" by default, which produces fewer interactions but
#   those interactions tend to be strong.
# `trimean` approximates 25% truncated mean, meaning avg gene expression is 0 if % of expressed
#   cells in a group is less than 25%.
# Other options used include 5% and 10%.
# `computeAvgExpr` lets you calculate the avg expression of a gene of interest; useful
#   for optimizing the truncation threshold.
## Compute ===========
cellchat <- computeCommunProb(cellchat, type = "triMean")
# 1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-10-21 15:16:02.583159]"
#   |======================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-10-21 15:18:55.365808]" # nolint

## Filter ===========
# Minimum number of cells required for each cell group. Default = 10
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract network as a df -----------------------------------------------
# All inferred communications at ligand/receptor level
df.net <- subsetCommunication(
  cellchat,
  slot.name = "net" # "netP" if you want signaling pathways
)
# source.use = c(1,2) and targets.use = c(4,5) to get communications sending from
#   cell groups 1 and 2, and recieving from groups 4 and 5.
# signaling = c("WNT", "TGFb") to get communications mediated by WNT and TGFb

# Infer communications at pathway level -----------------------------------------------
# Summarize communication probabilities of interactions associated with each signaling
#   pathway.
cellchat <- computeCommunProbPathway(cellchat)

# Calculate aggregated CCC network -----------------------------------------------
# Count number of links / summarize communication probability.
# Can also perform for only a subset of groups via `source.use` and `targets.use`.
cellchat <- aggregateNet(cellchat)

## Visualize aggregated network ===========
groupSize <- as.numeric(table(cellchat@idents))
pdf(file.path(plotsDir, "allCirclePlot.pdf"), width = 7, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)
dev.off()

## Visualize network per cell group ===========
# Also controlling edge weight max, so they're all relative to the overall maximum
pdf(file.path(plotsDir, "indivGroups_allCirclePlot.pdf"), width = 10, height = 10)
mat <- cellchat@net$weight # get overall max edge weight
par(mfrow = c(3, 4), xpd = TRUE)
for (idx in seq_len(nrow(mat))) {
  # Make a new matrix for the cell-group subset
  mat2 <- matrix(0, nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[idx, ] <- mat[idx, ]
  # Plot
  netVisual_circle(
    mat2,
    vertex.weight =  groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[idx]
  )
}
dev.off()

# Save -----------------------------------------------
saveRDS(cellchat, file.path(rDataDir, "cellchat_part2.rds"))

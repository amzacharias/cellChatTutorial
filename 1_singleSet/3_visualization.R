#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Part 3: Visualization of CCC network
# Author: Amanda Zacharias
# Date: 2024-10-21
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# Multiple ways to visualize the network
# Tools to extract and visualize high-order info e.g., prediction of major signaling inputs & outputs
# Quantitatively characterize CCC networks combining social network analysis, pattern recognition,
#   and manifold learning.

# Options -----------------------------------------------
options(stringsAsFactors = FALSE)

# Packages -----------------------------------------------
library(CellChat) # 2.1.2

# Source -----------------------------------------------

# Pathways -----------------------------------------------
## Input ===========
dataDir <- "0_data"

## Output ===========
baseDir <- file.path(getwd(), "1_singleSet")
rDataDir <- file.path(baseDir, "rDataDir")
plotsDir <- file.path(baseDir, "plots")

dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
cellchat <- readRDS(file.path(rDataDir, "cellchat_part2.rds"))

# Visualize each pathway -----------------------------------------------
# edge colours = source / sender
# edge weights = interaction strength. thicker = stronger signal
# circle size (hierarchy & circle) = # of cells in each group
# solid and open circles (hierarchy) = source and target respectively
# thin bar colours (chord) = targets that recieve signal from outer bar
# inner bar size (chord) proportion to signal strength recieved by target

# netVisual_aggregate = signaling pathways level
# netVisual_individual = L-R pairs associated with pathway

## See all sig. pathways ===========
cellchat@netP$pathways
#  [1] "MIF"        "ANNEXIN"    "CypA"       "GALECTIN"   "CXCL"
#  [6] "COMPLEMENT" "FGF"        "TNF"        "CCL"        "GAS"
# [11] "IL4"        "CD40"       "IGFBP"      "LIGHT"      "CSF"
# [16] "VEGF"

## Hierarchical Plot ===========
# `vertex.receiver`: numeric, index of cell groups as target on left
# right of plot shows signaling to remaining groups.
pathways.show <- c("CXCL") # pathway of interest
# signaling to fibroblast on left, right = signaling to immune cells
vertex.reciever <- seq(1, 4) # numeric vector
pdf(file.path(plotsDir, "CXCL_fibroImmune_hierPlot.pdf"), width = 7, height = 7)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver  = vertex.receiver)
dev.off()

## Chord diagram ===========
pdf(file.path(plotsDir, "CXCL_chordPlot.pdf"), width = 7, height = 7)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

## Heatmap ===========
pdf(file.path(plotsDir, "CXCL_heatmapPlot.pdf"), width = 7, height = 7)
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

## Chord custom ===========
# Customize using the `circlize` package
# Group cell clusters into different cell types
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
names(group.cellType) <- levels(cellchat@idents)
pdf(file.path(plotsDir, "CXCL_groupCells_chordPlot.pdf"), width = 7, height = 7)
netVisual_chord_cell(
  cellchat,
  signaling = pathways.show,
  group = group.cellType,
  title.name = paste0(pathways.show, " signaling network")
)
dev.off()

# Contrib of each pair to pathway -----------------------------------------------
netAnalysis_contribution(cellchat, signaling = pathways.show)
# extractEnrichedLR to extract all sig. L-R pairs and related genes in a pathway
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1, ] # get the first L-R pair
# Can plot as before except with the netVisual_individual() function...

## Save all of network with a for loop ===========
# ! setwd for saving plots in next step
setwd(plotsDir)
# All sig pathways
pathways.show.all <- cellchat@netP$pathways
# Make sure order of cells is correct
levels(cellchat@idents)
vertex.receiver <- seq(1, 4) # interested in APOE+ FIB and Inflam. FIB
for (idx in seq_along(pathways.show.all)) {
  # Visualize signaling pathway; will output PDFs in the cwd()
  netVisual(
    cellchat,
    signaling = pathways.show.all[idx],
    vertex.receiver = vertex.receiver,
    layout = "hierarchy",
    out.format = "pdf"
  )
  # Compute & visualize L-R contribution
  tmpGgPlot <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[idx])
  ggsave(
    plot = tmpGgPlot,
    filename = paste0(pathways.show.all[idx], "_L-R_contribution.pdf"), path = plotsDir,
    width = 90, height = 90, units = "mm"
  )
}

# ! resetting wd
setwd(file.path("..", baseDir))

# CCC mediated by multiple L-Rs or pathways -----------------------------------------------
## Bubble plot ===========
# L-R pairs from some cell groups (sources.use) to other groups (targets.use)
# remove.isolate = whether to remove entire empty columns
bubblePlot_4vs5to11 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
ggsave(
  plot = bubblePlot_4vs5to11,
  filename = "bubblePlot_4vs5to11.pdf", path = plotsDir,
  width = 180, height = 180, units = "mm"
)

# In addition, look at certain signaling pathways
bubblePlot_4vs5to11_CCLtoCXCL <- netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = c(5:11),
  remove.isolate = FALSE,
  signaling = c("CCL", "CXCL")
)
ggsave(
  plot = bubblePlot_4vs5to11_CCLtoCXCL,
  filename = "bubblePlot_4vs5to11_CCLtoCXCL.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

# Also, only certain L-R pairs
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL", "CXCL", "FGF"))
netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = c(5:11),
  remove.isolate = FALSE,
  pairLR.use = pairLR.use
)
# Not saving b.c. I don't think it's necessary
# Other customizations are available... e.g., sorting of cell pairs

## Chord Diagram ===========
# Show all L-R pairs from cell groups to other cell groups
# E.g. all interactions sending from Inflam.FIB
netVisual_chord_gene(
  cellchat,
  sources.use = 4,
  targets.use = c(5:11),
  lab.cex = 0.5,
  legend.pos.y = 30
)
dev.off()
# All recieved by Inflam.DC
netVisual_chord_gene(
  cellchat,
  sources.use = c(1, 2, 3, 4),
  targets.use = 8,
  legend.pos.x = 15
)
dev.off()

# All L-R pairs assoc with CCL and CXCL
netVisual_chord_gene(
  cellchat,
  sources.use = c(1, 2, 3, 4),
  targets.use = c(5:11),
  signaling = c("CCL", "CXCL"),
  legend.pos.x = 8
)
dev.off()

# All sig pathways from cell groups to cell groups
netVisual_chord_gene(
  cellchat,
  sources.use = c(1, 2, 3, 4),
  targets.use = c(5:11),
  slot.name = "netP",
  legend.pos.x = 8
)
dev.off()

# Plot signaling gene expression using violin/dot plot -----------------------------------------------
# Uses Seurat's plotGeneExpression function. 
# Three visualization types: violin, dot, or bar.
# Can also extract genes using extractEnrichedLR, and do it yourself.
violinPlot_CXCL <- plotGeneExpression(
  cellchat,
  signaling = "CXCL",
  enriched.only = TRUE,
  type = "violin"
)
ggsave(
  plot = violinPlot_CXCL,
  filename = "violinPlot_CXCL.pdf", path = plotsDir,
  width = 180, height = 180, units  = "mm"
)

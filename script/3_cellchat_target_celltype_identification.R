###########################################################
# CELLCHAT TARGET-CELLTYPE IDENTIFICATION
###########################################################

# A cell type with high levels of intracellular communication with the macrophage subset must be selected
# this is a screening step, so steps are taken to increase speed while maintaining decent accuracy


##############################
# Load Libraries
##############################

library(RColorBrewer)
library(tidyverse)
library(Seurat)
library(CellChat)
library(rstudioapi)
library(reticulate)
library(pheatmap)
library(ggpubr)
library(cowplot)

# Colours!
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")
# set the random seed as 42 for consistent colour generation
set.seed(42)


##############################
# Set Working Directory
##############################

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


##############################
# Load in the Dataset
##############################

# load pbmc
pbmc <- readRDS("../data/pbmc_final.rds")

# retrieve only the cells from the tumour samples
pbmc_T <- subset(pbmc, SPECIMEN_TYPE %in% c("T"))

# save 
# saveRDS(pbmc_T, file = "../data/pbmc_T.rds")
# pbmc_T <- readRDS("../data/pbmc_T.rds")


##############################
# Dataset 1 preparation
##############################

# Make a new dataset with alternative labels
df1 <- pbmc_T

# check the label distribution
table(df1@meta.data$clMidwayPr[df1@meta.data$clTopLevel == "Myeloid"])
table(df1@meta.data$cl295v11SubFull[df1@meta.data$clTopLevel == "Myeloid"])
table(df1@meta.data$cell_type_detailed[df1@meta.data$clTopLevel == "Myeloid"])

# use the top level labels, but with myeloid cells removed and replaced with SPP1 and C1QC level labels
# Add CellChat labels first
df1@meta.data$cclabels <- ifelse(
  df1@meta.data$clTopLevel == "Myeloid",
  df1@meta.data$cell_type_detailed,
  df1@meta.data$clTopLevel
)
table(df1@meta.data$cclabels)

# Keep non-myeloids + selected macrophage subsets
keep_cells <- rownames(df1@meta.data)[
  df1@meta.data$clTopLevel != "Myeloid" |
    df1@meta.data$cell_type_detailed %in% c(
      "Myeloid: Macro-SPP1",
      "Myeloid: Macro-C1QC")]

df1 <- subset(df1, cells = keep_cells)
table(df1@meta.data$cclabels)

# downsample the dataset to reduce running time for preliminary stage
print(ncol(df1))

df1 <- subset(
  df1,
  downsample = 3000 # change this number depending on the time it takes for the next steps
)
print(ncol(df1))

# add samples and change identities
df1$cclabels <- droplevels(factor(df1$cclabels))
Idents(df1) <- df1@meta.data$cclabels
df1$samples <- factor(df1$orig.ident)


###########################################################
# Cellchat
###########################################################
##############################
# # Part I: Data input & processing and initialization of CellChat object
##############################

# create the cellchat object
cellchat <- createCellChat(object = df1, 
                           group.by = "cclabels",
                           assay = "RNA")

# load database
CellChatDB <- readRDS("../data/cellchat_db.rds")

# use all CellChatDB for cell-cell communication analysis
cellchat@DB <- CellChatDB

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 

# identify over-expressed signalling genes in each group
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = TRUE)

# identify over expressed ligand-receptor interactions
cellchat <- identifyOverExpressedInteractions(cellchat)


##############################
# # Part II: Inference of cell-cell communication network
##############################

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, nboot = 20)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)


##############################
# Circle plot
##############################

groupSize <- as.numeric(table(cellchat@idents))

# to decrease noise, keep only interactions and directions of interest for count and weight
# count 
mat.count.out <- cellchat@net$count
# Keep only SPP1/C1QC as senders
mat.count.out[!rownames(mat.count.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0
# Remove SPP1/C1QC as receivers
mat.count.out[, colnames(mat.count.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC")] <- 0

# weight
mat.weight.out <- cellchat@net$weight
# Keep only SPP1/C1QC as senders
mat.weight.out[!rownames(mat.weight.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0
# Remove SPP1/C1QC as receivers
mat.weight.out[, colnames(mat.weight.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC")] <- 0

# show circle plot
# interaction counts
netVisual_circle(
  mat.count.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interactions from Macro-SPP1/C1QC"
)

# show circle plot
# interaction weights
netVisual_circle(
  mat.weight.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interaction strength from Macro-SPP1/C1QC"
)


##############################
# heatmap
##############################

# reformat names for figures
mat.weight.out <- cellchat@net$weight[c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), 
                                      !grepl("^Myeloid", colnames(cellchat@net$weight))]
mat.count.out <- cellchat@net$count[c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"),
                                    !grepl("^Myeloid", colnames(cellchat@net$count))]

a4 <- pheatmap::pheatmap(
  mat.weight.out,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  main = "Interaction Weight",
  color = brewer.pal(9, "YlOrRd"))
a4

b4 <- pheatmap::pheatmap(
  mat.count.out,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  main = "Interaction Count",
  brewer.pal(9, "YlOrRd")
)
b4

# save the cellchat object for later visualisation
# saveRDS(cellchat, file = "../data/cellchat_MtoAll.rds")
# cellchat <- readRDS("../data/cellchat_MtoAll.rds")


##############################
# Analysis of SPP1/C1QC interactions with T cell subtypes
##############################

# T cells were of the strongest outgoing interactions
# initially load in all T cells to investigate the strongest interacting groups

# check the number of each specific cell type are sufficient
table(pbmc_T@meta.data$cell_type_detailed[pbmc_T@meta.data$clTopLevel %in% c("TNKILC")])


# Make a new dataset with alternative labels
df2 <- pbmc_T

# use the top level labels, but with myeloid cells removed and replaced with SPP1 and C1QC level labels
# Add CellChat labels first
df2@meta.data$cclabels <- ifelse(
  df2@meta.data$clTopLevel == "Myeloid",
  df2@meta.data$cell_type_detailed,
  df2@meta.data$clTopLevel
)

# subset into my second dataframe
# Keep non-myeloids + selected macrophage subsets
keep_cells <- rownames(df2@meta.data)[
  df2@meta.data$clTopLevel == "TNKILC" |
    df2@meta.data$cell_type_detailed %in% c(
      "Myeloid: Macro-SPP1",
      "Myeloid: Macro-C1QC")]

df2 <- subset(df2, cells = keep_cells)

# T cell breakdown
  # Tgd cells
  # "TNKILC: gd-like T" "TNKILC: gd-like T PDCD1+"
  # "TNKILC: gd-like T prolif"
  # TZBTB16 cells
  # "TNKILC: PLZF+ T"  "TNKILC: PLZF+ T prolif"  
  # Macrophages
  # "Myeloid: Macro-C1QC"  "Myeloid: Macro-SPP1" 
  # CD4 cells
  # "TNKILC: CD4+ CXCL13+"  "TNKILC: CD4+ IL17+" 
  # "TNKILC: CD4+ IL7R+"  "TNKILC: CD4+ IL7R+CCL5+" 
  # "TNKILC: CD4+ IL7R+HSP+"  "TNKILC: CD4+ IL7R+SELL+" 
  # "TNKILC: CD4+ TFH"  "TNKILC: CD4+ Treg" 
  # "TNKILC: CD4+ Treg prolif" 
  # CD8 cells
  # "TNKILC: CD8+ CXCL13+"  "TNKILC: CD8+ CXCL13+ HSP+"  
  # "TNKILC: CD8+ CXCL13+ prolif"  "TNKILC: CD8+ IL7R+"  
  # "TNKILC: CD8+ T IL17+" "TNKILC: CD8+GZMK+"  
  # ILC cells
  # "TNKILC: ILC3"  
  # NK cells
  # "TNKILC: cTNI22"  "TNKILC: NK CD16A+"  
  # "TNKILC: NK GZMK+"  "TNKILC: NK XCL1+"

# add the midway annotations of the T cells 
# and the detailed annotations of macrophages
df2@meta.data$cclabels <- ifelse(
  df2@meta.data$clTopLevel == "Myeloid",
  df2@meta.data$cell_type_detailed,
  df2@meta.data$clMidwayPr
)

# set these labels as the identities
Idents(df2) <- df2@meta.data$cclabels

# downsample to reduce computational cost
table(df2@meta.data$cclabels)
df2 <- subset(
  df2,
  downsample = 500 
)
table(df2@meta.data$cclabels)

# add a samples column for cellchat 
df2$samples <- factor(df2$orig.ident)


##############################
# Data input & processing and initialization of CellChat object
##############################

# create the cellchat object
cellchat <- createCellChat(object = df2, 
                           group.by = "cclabels",
                           assay = "RNA") # use lognormalised assay data

# load database and select interactions
CellChatDB <- readRDS("../data/cellchat_db.rds")

# use all CellChatDB for cell-cell communication analysis
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)

# identify over-expressed signalling genes in each group
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = TRUE)

# identify over expressed ligand-receptor interactions
cellchat <- identifyOverExpressedInteractions(cellchat)


##############################
# Inference of cell-cell communication network
##############################

# Compute the communication probability and infer cellular communication network
# nboot = 20 as a low power preliminary screen
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, nboot = 20)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)


##############################
# Circle Plots
##############################

# Visualise
groupSize <- as.numeric(table(cellchat@idents))

# keep only interactions and directions of interest four count and weight
# count  
mat.count.out <- cellchat@net$count
# Keep only SPP1/C1QC as senders
mat.count.out[!rownames(mat.count.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0
# Remove SPP1/C1QC as receivers
mat.count.out[, colnames(mat.count.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC")] <- 0

# weight
mat.weight.out <- cellchat@net$weight
# Keep only SPP1/C1QC as senders
mat.weight.out[!rownames(mat.weight.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0
# Remove SPP1/C1QC as receivers
mat.weight.out[, colnames(mat.weight.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC")] <- 0

# interaction counts
# show circle plot
netVisual_circle(
  mat.count.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interactions from Macro-SPP1/C1QC")

# interaction weights
# show circle plot
netVisual_circle(
  mat.weight.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interaction strength from Macro-SPP1/C1QC")


##############################
# heatmap
##############################

# reformat titles
mat.weight.out <- cellchat@net$weight[c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"),
                                      !grepl("^Myeloid", colnames(cellchat@net$weight))]
mat.count.out <- cellchat@net$count[c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"),
                                    !grepl("^Myeloid", colnames(cellchat@net$count))]

# create heatmap 
c4 <- pheatmap::pheatmap(mat.weight.out,
               cluster_rows = FALSE,
               cluster_cols = FALSE, 
               main = "Interaction Weight (T Cells)",
               brewer.pal(9, "YlOrRd"))
c4

d4 <- pheatmap::pheatmap(mat.count.out,
               cluster_rows = FALSE,
               cluster_cols = FALSE, main = "Interaction Count (T Cells)",
               brewer.pal(9, "YlOrRd"))
d4


##############################
# Make Figures
##############################

# Figure 4
ggdraw() +
  draw_plot(a4$gtable, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(b4$gtable, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(c4$gtable, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(d4$gtable, x = .5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15, 
                  x = c(0, .5, .0, .5), 
                  y = c(1, 1, .5, .5))
# save figure
ggsave("../figures/Figure_4_cc_general.jpg", width = 30, height = 15, units = c("cm"), dpi = 300)


##############################
# Final Save
##############################

# save the cellchat object for later visualisation
# saveRDS(cellchat, file = "../data/cellchat_MtoTall.rds")
# cellchat <- readRDS("../data/cellchat_MtoTall.rds")

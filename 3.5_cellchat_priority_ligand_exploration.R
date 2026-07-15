###########################################################
# CELLCHAT SPP1/C1QC - T PRIORITY LIGAND EXPLORATION
###########################################################

# identification of priority interacting ligands for nichenet analysis
# CD8 CD4 and Tgd were origionally selected and Tgd was later excluded from analysis


##############################
# Load libraries
##############################

library(BiocManager)
library(RColorBrewer)
library(future)
library(dplyr)
library(Seurat)
library(CellChat)
library(rstudioapi)
library(reticulate)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# Colours!
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")


##############################
# Set Working Directory
##############################

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


##############################
# Further subsetting
##############################

# as TCD8 and tgd have the strongest and most interactions with the macrophages
# their subtypes will be chosen for further downstream analysis
# in addition, TCD4 will will also be chosen for comparison with preliminary results


##############################
# HPC friendly cellchat data processing
##############################

# this step will be run on the HPC as the cell number is to computationally costly

print("running r script")

library(BiocManager)
library(RColorBrewer)
library(future)
library(dplyr)
library(Seurat)
library(CellChat)
library(rstudioapi)
library(reticulate)

print("libraries loaded")

# load in the dataset 
pbmc_T <- readRDS("../data/pbmc_T.rds")
# subset the Tgd and CD8 cells
df4 <- subset(
  pbmc_T,
  cell_type_detailed %in% c("Myeloid: Macro-SPP1",
                            "Myeloid: Macro-C1QC") | 
    clMidwayPr %in% c("TCD8", "TCD4", "Tgd"))

# add cclabels column
df4@meta.data$cclabels <- df4@meta.data$cell_type_detailed
# add make these the identity
Idents(df4) <- df4@meta.data$cclabels
# add patient numbers
df4$samples <- factor(df4$orig.ident)

print("dataset prepared")

# create the cellchat object
cellchat <- createCellChat(object = df4, 
                           group.by = "cclabels",
                           assay = "RNA")
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

print("cellchat object created")

# Compute the communication probability and infer cellular communication network
# nboot is increased for higher significance
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, nboot = 100)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

print("cellchat object processed")

saveRDS(cellchat, file = "../data/cellchat_MtoT_exploration.rds")

print("cellchat object saved. FINITO")


##############################
# End of HPC
##############################
##############################
# load HPC file
##############################

# read in the the file from the HPC
# cellchat <- readRDS("../data/cellchat_MtoT_exploration.rds")


##############################
# Circle Plots
##############################

# screening of cell subtypes 
# Visualise
groupSize <- as.numeric(table(cellchat@idents))

# remove all interactions other than the ones of interest
# Interaction counts
mat.count.out <- cellchat@net$count
# Remove all senders except the macrophages
mat.count.out[!rownames(mat.count.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0

# Interaction strengths
mat.weight.out <- cellchat@net$weight
# Remove all senders except the macrophages
mat.weight.out[!rownames(mat.weight.out) %in% c("Myeloid: Macro-SPP1","Myeloid: Macro-C1QC"), ] <- 0

# interaction counts
netVisual_circle(
  mat.count.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interaction counts from Macro-SPP1/C1QC"
)

# interaction weights
netVisual_circle(
  mat.weight.out,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Outgoing interaction strength from Macro-SPP1/C1QC"
)


##############################
# Source and target preparation
##############################

# classification of targets and source cells
sources = c("Myeloid: Macro-C1QC", "Myeloid: Macro-SPP1")

targets.cd8 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD8"])
targets.cd4 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD4"])
targets.gd <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "Tgd"])


##############################
# Visualization of cell-cell communication network
##############################
##############################
# Bubble plot
##############################

# CD8
netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = targets.cd8,
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 # customisation
                 color.heatmap = "viridis",
                 n.colors = 2, 
                 direction = 1,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()

# Tgd
netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = targets.gd,
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 # customisation
                 color.heatmap = "viridis",
                 n.colors = 2, 
                 direction = 1,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()

# CD4
netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = targets.cd4,
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 # customisation
                 color.heatmap = "viridis",
                 n.colors = 2, 
                 direction = 1,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()


##############################
# Chord diagram
##############################

# CD8
netVisual_chord_gene(cellchat, sources.use = c("Myeloid: Macro-C1QC", 
                                               "Myeloid: Macro-SPP1"), 
                     targets.use = targets.cd8,
                     slot.name = "netP",
                     show.legend = FALSE
                     )

# Tgd
netVisual_chord_gene(cellchat, sources.use = c("Myeloid: Macro-C1QC", 
                                               "Myeloid: Macro-SPP1"), 
                     targets.use = targets.tgd,
                     slot.name = "netP",
                     show.legend = FALSE
                     )
# CD4
netVisual_chord_gene(cellchat, sources.use = c("Myeloid: Macro-C1QC", 
                                               "Myeloid: Macro-SPP1"), 
                     targets.use = targets.tgd,
                     slot.name = "netP",
                     show.legend = FALSE
)


##############################
# Heat map
##############################

#	Compute and visualize the network centrality scores
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10)

# Visualize dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)


# Identify signals contributing the most to outgoing or incoming signaling of certain cell groups	
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")


##############################
# Global Patterns
##############################

# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
# Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# save the top half of this
netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of receiving cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# save the bottom half of this
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")


##############################
# Functional clustering
##############################

# install python
Sys.which("python")
use_python("C:/Users/mileh/MINICO~1/python.exe", required = TRUE)
py_config()
py_install("umap-learn", pip = TRUE)
import("umap")
py_module_available("umap")
options(future.globals.maxSize = 4 * 1024^3) 
plan(sequential)


### Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


##############################
# Saves
##############################

# save the final cellchat object
# saveRDS(cellchat, file = "../data/cellchat_MtoT_exploration.rds")
# cellchat <- readRDS("../data/cellchat_MtoT_exploration.rds")

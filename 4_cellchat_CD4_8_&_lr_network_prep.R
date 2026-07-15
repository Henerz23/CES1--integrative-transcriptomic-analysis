###########################################################
# CELLCHAT COMBINED FIGURE CREATION AND NICHENET SUBSET PREPARATION
###########################################################

# Part 1: final cellchat object and L-R selection based on CD4 CD8 interactions
# Part 2: L-R pairs distributed into datasets depending on the amount of interaction per sybtype pair
  # this helps guide more informed nichenet analysis of distinct macrophage- T- cell interaction
# Part 3: bubble plot orgarnised into interaction types created


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


###########################################################
# Part 1:
###########################################################
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
df5 <- subset(
  pbmc_T,
  cell_type_detailed %in% c("Myeloid: Macro-SPP1",
                            "Myeloid: Macro-C1QC") | 
    clMidwayPr %in% c("TCD8", "TCD4"))

# use the top level labels
df5@meta.data$cclabels <- df5@meta.data$cell_type_detailed

# add make these the identity
Idents(df5) <- df5@meta.data$cclabels
# add patient numbers
df5$samples <- factor(df5$orig.ident)

print("dataset prepared")

# create the cellchat object
cellchat <- createCellChat(object = df5, 
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

# save cellchat object
# saveRDS(cellchat, file = "../data/cellchat_MtoCD8_4.rds")

print("cellchat object saved. FINITO")


##############################
# End of HPC
##############################


###########################################################
# Part 2:
###########################################################
##############################
# load HPC file
##############################

# read in the the file from the HPC
cellchat <- readRDS("../data/cellchat_MtoCD8_4.rds")
pbmc <- readRDS("../data/pbmc_final.rds")

# redefine sources and targets
sources = c("Myeloid: Macro-C1QC", "Myeloid: Macro-SPP1")
targets.cd8 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD8"])
targets.cd4 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD4"])


##############################
# manual distinction of LR pairs
##############################

# LR pairs will be distinguished into 6 dfs, depending on:
  # source cell weighted more towards SPP1, C1QC or neither (both)
  # whether the interaction is present in CD8 or CD4
    # the dsitinction of a both CD8 and CD4 isnt neccesary for nichenet as 
    # the genes upregulated are distinct and therefore a "both CD8 and CD4" input
    # would incorrectly assume that upregulated genes would be identical for identical L-R pairs
  # however, for an aesthetic distinction in the bubble plot, a both CD8 and CD4 section
  # of LR pairs equally present in both combinations will be used

## T-cell targets for easier function input
targets <- list(TCD8 = targets.cd8, TCD4 = targets.cd4)

# Function 1:
# retrives the LR pairs that are expressed in that source to target combination
# makes 4 total dfs
   # e.g. all LR interactions between SPP1 and CD4
get_lri <- function(source, targets) {
  # subset communication from the cellchat object, defining all sources and targets
  subsetCommunication(
    cellchat,
    sources.use = source,
    targets.use = targets) %>%
    # select communication columns needed (prob being communication probability)
    select(interaction_name, ligand, receptor, prob) %>%
    # rename for nichenet compatability
    rename(from = ligand, to = receptor) %>%
    
    # collapse rows with the same interactions into one 
    group_by(interaction_name, from, to) %>%
    summarise(
      # mean the probs of LR interactions of each T cell types
      prob = mean(prob),
      .groups = "drop"
    )
}

# Function 2: 
# define function that compares lr interactions from SPP1 and C1QC macrophages
# fold is the threshold to decide which interaction is stronger than the other
classify_lr <- function(spp1, c1qc, fold = 2) {
  
  # combines the two interaction tables by matching interactions
  both <- full_join(spp1, c1qc,
    by = c("from", "to", "interaction_name"),
    suffix = c("_spp1", "_c1qc"))
  # reactions not found in one celltype are left as NA
  
  # list returned with 3 dataframes with only the LR columns
  list(
    spp1 = both %>%
      # keep interactions that are absent in C1QC or
      # the fold times stronger in SPP1 than C1QC
      filter(is.na(prob_c1qc) | prob_spp1 >= fold * prob_c1qc) %>%
      select(interaction_name, from, to),
    
    c1qc = both %>%
      # keep interactions that are absent in SPP1 or
      # the fold times stronger in C1QC than SPP1
      filter(is.na(prob_spp1) | prob_c1qc >= fold * prob_spp1) %>%
      select(interaction_name, from, to),
    
    both = both %>%
      # keep interactions in both macrophage types
      filter(!is.na(prob_spp1), !is.na(prob_c1qc),
             # where neither is stronger than the other
        prob_spp1 < fold * prob_c1qc,
        prob_c1qc < fold * prob_spp1
      ) %>%
      select(interaction_name, from, to))}


# Create 6 dataframes. loop over targets
for (nm in names(targets)) {
  
  # use function 1 to find all interactions per cell type pair
  spp1 <- get_lri("Myeloid: Macro-SPP1", targets[[nm]])
  c1qc <- get_lri("Myeloid: Macro-C1QC", targets[[nm]])
  
  # assign a specific df name
  assign(paste0("lri_spp1_", nm), spp1)
  assign(paste0("lri_c1qc_", nm), c1qc)
  
  # use function 2 to classify into the 3 groups per T cell type
  cls <- classify_lr(spp1, c1qc, fold = 2)
  
  # assign specific names to these- 6 total
  assign(paste0("lr_spp1_", nm), cls$spp1)
  assign(paste0("lr_c1qc_", nm), cls$c1qc)
  assign(paste0("lr_both_", nm), cls$both)
}

# add these dataframes to a list
lr_networks <- list(
  lr_spp1_TCD8 = lr_spp1_TCD8,
  lr_spp1_TCD4 = lr_spp1_TCD4,
  lr_both_TCD8 = lr_both_TCD8,
  lr_both_TCD4 = lr_both_TCD4,
  lr_c1qc_TCD8 = lr_c1qc_TCD8,
  lr_c1qc_TCD4 = lr_c1qc_TCD4
)

# then create a cellchat specific one which has overlaps in CD4 and CD8- 9 groups total
lr_networks9 <- list(
  # spp1
  spp1_TCD8 = anti_join(lr_spp1_TCD8, lr_spp1_TCD4, by = c("interaction_name", "from", "to")),
  spp1_both = inner_join(lr_spp1_TCD8, lr_spp1_TCD4, by = c("interaction_name", "from", "to")),
  spp1_TCD4 = anti_join(lr_spp1_TCD4, lr_spp1_TCD8, by = c("interaction_name", "from", "to")),
  # both 
  both_TCD8 = anti_join(lr_both_TCD8, lr_both_TCD4, by = c("interaction_name", "from", "to")),
  both_both = inner_join(lr_both_TCD8, lr_both_TCD4, by = c("interaction_name", "from", "to")),
  both_TCD4 = anti_join(lr_both_TCD4, lr_both_TCD8, by = c("interaction_name", "from", "to")),
  # c1qc
  c1qc_TCD8 = anti_join(lr_c1qc_TCD8, lr_c1qc_TCD4, by = c("interaction_name", "from", "to")),
  c1qc_both = inner_join(lr_c1qc_TCD8, lr_c1qc_TCD4, by = c("interaction_name", "from", "to")),
  c1qc_TCD4 = anti_join(lr_c1qc_TCD4, lr_c1qc_TCD8, by = c("interaction_name", "from", "to")))


###########################################################
# Part 3:
###########################################################
##############################
# Figure Preparation
##############################

# 1. Extract enriched L-R pairs for your signaling pathways of interest
ccnn_order <- data.frame(
  interaction_name = unlist(
    lapply(lr_networks9, function(x) x$interaction_name),
    use.names = FALSE
  )
)

# create figure ready organised LR interaction bubble plot
a5 <- netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = c(targets.cd4, targets.cd8),
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 pairLR.use = ccnn_order,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()
a5


##############################
# Saves
##############################

# save the final cellchat object
# saveRDS(cellchat, file = "../data/cellchat_MtoCD8_4.rds")
# cellchat <- readRDS("../data/cellchat_MtoCD8_4.rds")

# saveRDS(lr_networks, file = "../data/lr_networks.rds")
# lr_networks <-  readRDS("../data/lr_networks.rds")

# saveRDS(lr_networks9, file = "../data/lr_networks9.rds")
# lr_networks9 <- readRDS("../data/lr_networks9.rds")

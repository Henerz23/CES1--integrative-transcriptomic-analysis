###########################################################
# STROMAL FULL ANALYSIS
###########################################################

# Part 1: 
# Part 2: 


##############################
# load libraries
##############################

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrastr)
library(Hmisc)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(stringr)
library(cowplot)

# Colours!
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")


###########################################################
# Part 1
###########################################################
##############################
# Analysis of CES1
##############################

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load in pbmc
pbmc <- readRDS("../data/pbmc_final.rds")

# create stromal subset 
table(pbmc@meta.data$clTopLevel)
pbmc_stromal <- subset(pbmc, subset = clTopLevel %in% "Strom")

# analyse CES1
Idents(pbmc_stromal) <- "clMidwayPr"
DotPlot(pbmc_stromal, features = 'CES1', cols = cont_2) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# change names to only contain the cell types
pbmc_stromal@meta.data$cl295v11SubFull <- 
  gsub(".*\\((.*)\\).*", "\\1", 
       pbmc_stromal@meta.data$cl295v11SubFull)

# analyse Fibro
pbmc_fibro <- subset(pbmc_stromal, subset = clMidwayPr %in% "Fibro")
Idents(pbmc_fibro) <- "cell_type_detailed"
DotPlot(pbmc_fibro, features = 'CES1', cols = cont_2) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


# violin plot
VlnPlot(pbmc_fibro, features = c("CES1"), ncol = 3)


##############################
# Data Preprocessing and Clustering
##############################

## normalization
pbmc_fibro <- NormalizeData(pbmc_fibro, normalization.method = "LogNormalize", scale.factor = 10000) 
pbmc_fibro <- FindVariableFeatures(pbmc_fibro, selection.method = 'vst', nfeatures = 2000)

# scaling
pbmc_fibro <- ScaleData(pbmc_fibro, vars.to.regress = "percent.mt")
pbmc_fibro <- RunPCA(pbmc_fibro, features = VariableFeatures(object = pbmc_fibro)) 

# dim use identification
ElbowPlot(pbmc_fibro, reduction = "pca", ndims = 50)

# dimensional reduction
dim.use<-1:20
pbmc_fibro <- FindNeighbors(pbmc_fibro, dims = dim.use)
pbmc_fibro <- FindClusters(pbmc_fibro, resolution = 0.8 )
pbmc_fibro <- RunUMAP(pbmc_fibro, dims = dim.use)

# plot clusters
a8 <- DimPlot(pbmc_fibro, reduction = 'umap',label=F, 
              group.by = 'cl295v11SubFull', cols = distinctColorPalette(11)) +
  ggtitle(NULL)
a8
Idents(pbmc_fibro) <- "cl295v11SubFull"


##############################
# Finding differentially expressed genes
##############################

# find the markers for the fibroblast clusters
pbmc_fibro.markers <- FindAllMarkers(pbmc_fibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_fibro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# saveRDS(pbmc_fibro.markers, file = "../data/pbmc_fibro.markers.rds")
# pbmc_fibro.markers <- readRDS("../data/pbmc_fibro.markers.rds")


##############################
# fibroblast/Meta data Check
##############################

# patient distribution
DimPlot(object = pbmc_fibro, group.by="orig.ident", reduction='umap')

# metastasis status
# one group is called:
# pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum.
# change to pM1c
pbmc_fibro@meta.data$MetastasisStatus[pbmc_fibro@meta.data$MetastasisStatus == 
                                             "pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum."] <- "pM1c"
DimPlot(object = pbmc_fibro, group.by="MetastasisStatus",reduction='umap')

# tissue site
DimPlot(object = pbmc_fibro, group.by="TissueSiteSimple",reduction='umap')

# sex
DimPlot(object = pbmc_fibro, group.by="Sex",reduction='umap')

# N or T
DimPlot(object = pbmc_fibro, group.by="SPECIMEN_TYPE",
              reduction='umap', cols = distinctColorPalette(4))

# mmrstatus
DimPlot(object = pbmc_fibro, group.by="MMRStatus",reduction='umap')


##############################
# Analysis of CES1 in clusters
##############################

# dot plot of the fibroblast cells
DotPlot(pbmc_fibro, features = 'CES1', cols = cont_2, group.by = "cl295v11SubFull") + RotatedAxis()
# fibro stem cell niche is high but these are alll normal cells so make CAF only calss

# dotplot of just the CAFs
pbmc_CAF <- subset(pbmc_fibro, cl295v11SubFull %in% c("CXCL14+ CAF", "GREM1+ CAF", "MMP3+ CAF", "Myofibro"))
DotPlot(pbmc_CAF, features = 'CES1', cols = cont_2, group.by = "cl295v11SubFull") + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())



# bar chart showing the distribution of normal and tumour cells between fibroblasts 

b8 <- ggplot(
  pbmc_CAF@meta.data %>% filter(CES1_isExpressed),
  aes(x = cl295v11SubFull, fill = SPECIMEN_TYPE)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = disc_10) +
  labs(
    x = " ",
    y = "CES1 Expression Distribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x=element_blank())

b8


##############################
# CES1 average expression
##############################

# make a box plot of average CES1 expression for each patient and group by cl295v11SubFull
df_avg_CES1_exp_f <-  pbmc_CAF@meta.data %>%
  # group by
  group_by(orig.ident, cl295v11SubFull, MMRStatus) %>%
  # replace the NAs with Normal
  mutate(MMRStatus = ifelse(is.na(MMRStatus), "Normal", MMRStatus)) %>%
  # the summarise function
  summarise(avg_CES1_exp = mean(CES1_expression))

# define comparisons for statistical test
f_comparisons <- list(c("Myofibro", "CXCL14+ CAF"), 
                      c("Myofibro", "GREM1+ CAF"), 
                      c("Myofibro", "MMP3+ CAF"))

# make a boxplot depicting the average CES1 expression per patient
c8 <- ggplot(df_avg_CES1_exp_f, aes(x = cl295v11SubFull, y = avg_CES1_exp, fill = cl295v11SubFull)) +
  # 
  geom_boxplot(position="dodge", notch=TRUE, notchwidth = 0.8) +
  # 
  scale_fill_manual(values = disc_10) +
  labs(x = "Cell type", y = "Average Expression") +
  stat_compare_means(comparisons = f_comparisons, 
                     aes(label = after_stat(p.signif)),
                     method = "t.test") +
  theme_minimal() +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank())
c8

# make a boxplot depicting the avg ces1 exp for just the myofib class
ggplot(df_avg_CES1_exp_f %>% filter(cl295v11SubFull %in% c("Myofibro")),
             aes(x = cl295v11SubFull, y = avg_CES1_exp, fill = MMRStatus)) +
  # 
  geom_boxplot(notch=TRUE, notchwidth = 0.8) +
  # 
  scale_fill_manual(values = disc_10) +
  labs(x = "Cell type", y = "Average Expression") +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())



##############################
# CES1 percentage expressed
##############################

# make a box plot of Percentage Expressed for each patient and group by Caf type
# first make a summary dataframe
exp_pct_S_C_f <-  pbmc_CAF@meta.data %>%
  ##group by
  group_by(orig.ident, cl295v11SubFull, MMRStatus) %>%
  # replace the NAs with Normal
  mutate(MMRStatus = ifelse(is.na(MMRStatus), "Normal", MMRStatus)) %>%
  ##the summarise function
  summarise(percent_CES1_exp = mean(CES1_isExpressed) * 100)

# make a boxplot depicting the percentage CES1 expressed per patient
d8 <- ggplot(exp_pct_S_C_f, aes(x = cl295v11SubFull, y = percent_CES1_exp, fill = cl295v11SubFull)) +
  # 
  geom_boxplot(notch=TRUE, notchwidth = 0.8) +
  scale_fill_manual(values = disc_10) +
  labs(x = "Cell type", y = "Percent Expressed") +
  # stats test
  stat_compare_means(comparisons = f_comparisons, 
                     aes(label = after_stat(p.signif)),
                     method = "t.test") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())
d8

# make a bar chart depicting the MMRstatus and percentage of cells expresing CES1
ggplot(exp_pct_S_C, aes(x = myeloid_type, y = percent_CES1_exp, fill = MMRStatus)) +
  scale_fill_discrete("MMRStatus") +
  # 
  geom_boxplot(position="dodge", notch=TRUE, notchwidth = 0.8) +
  # 
  scale_fill_manual(values = disc_10) +
  labs(x = "Cell type", y = "% of cells expressing CES1") +
  # stats test
  stat_compare_means(method = "t.test", paired = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


##############################
# Final Saves
##############################

# saveRDS(pbmc_stromal, file = "../data/pbmc_stromal.rds")
# pbmc_stromal <- readRDS("../data/pbmc_stromal.rds")

# saveRDS(pbmc_fibro.markers, file = "../data/pbmc_fibro.markers.rds")
# pbmc_fibro.markers <- readRDS("../data/pbmc_fibro.markers.rds")

# saveRDS(pbmc_fibro, file = "../data/pbmc_fibro_final.rds")
# pbmc_fibro <- readRDS("../data/pbmc_fibro_final.rds")

# saveRDS(pbmc_CAF, file = "../data/pbmc_CAF.rds")
# pbmc_CAF <- readRDS("../data/pbmc_CAF.rds")


###########################################################
# Part 1
###########################################################
##############################
# Cellchat
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
# subset the Fibro
df6 <- subset(
  pbmc_T,
  cell_type_detailed %in% c("Myeloid: Macro-SPP1",
                            "Myeloid: Macro-C1QC") | 
    clMidwayPr %in% c("Fibro"))

# add cclabels column
df6@meta.data$cclabels <- df6@meta.data$cell_type_detailed
# add make these the identity
Idents(df6) <- df6@meta.data$cclabels
# add patient numbers
df6$samples <- factor(df6$orig.ident)

print("dataset prepared")

# create the cellchat object
cellchat <- createCellChat(object = df6, 
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

saveRDS(cellchat, file = "../data/cellchat_MtoF_exploration.rds")

print("cellchat object saved. FINITO")


##############################
# End of HPC
##############################
##############################
# load HPC file
##############################

# read in the the file from the HPC
# cellchat <- readRDS("../data/cellchat_MtoF_exploration.rds")


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
targets.fibro <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "Fibro"])


##############################
# Visualization of cell-cell communication network
##############################
##############################
# Bubble plot
##############################

# Fibro
netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = targets.fibro,
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()

table(pbmc_fibro@meta.data$SPECIMEN_TYPE[pbmc_fibro@meta.data$clMidwayPr == "Fibro"])


##############################
# manual distinction of LR pairs
##############################

## T-cell targets for easier function input
targets <- list(Fibro = targets.fibro)

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
lr_networks_Fibro <- list(
  lr_spp1_Fibro = lr_spp1_Fibro,
  lr_both_Fibro = lr_both_Fibro,
  lr_c1qc_Fibro = lr_c1qc_Fibro
)

# 1. Extract enriched L-R pairs for your signaling pathways of interest
ccnn_order <- data.frame(
  interaction_name = unlist(
    lapply(lr_networks_Fibro, function(x) x$interaction_name),
    use.names = FALSE
  )
)

# create ordered figure
a7 <- netVisual_bubble(cellchat, sources.use = sources, 
                 targets.use = c(targets.fibro),
                 remove.isolate = TRUE, sort.by.source = T, sort.by.target = T,
                 pairLR.use = ccnn_order,
                 angle.x = 45,
                 vjust.x = NULL,
                 hjust.x = NULL) +
  
  coord_flip() +
  RotatedAxis()
a7
##############################
# Saves
##############################

# save the final cellchat object
# saveRDS(cellchat, file = "../data/cellchat_MtoF_exploration.rds")
# cellchat <- readRDS("../data/cellchat_MtoF_exploration.rds")

# saveRDS(lr_networks_Fibro, file = "../data/lr_networks_Fibro.rds")
# lr_networks_Fibro <- readRDS("../data/lr_networks_Fibro.rds")


###########################################################
# Part 3
###########################################################
##############################
# HPC friendly nichenet processing step
##############################
##############################
# Load Libraries and data
##############################

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(patchwork)
library(cowplot)

print("libraries loaded")

ligand_target_matrix <- readRDS("../data/ligand_target_matrix_human.rds")
weighted_networks  <- readRDS("../data/weighted_networks_human.rds")
lr_networks_Fibro <-  readRDS("../data/lr_networks_Fibro.rds")
pbmc <- readRDS("../data/pbmc_final.rds")

print("data loaded")

# define targets
targets.fibro <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "Fibro"])

# read in the seurat object
seurat_obj <- subset(pbmc, cell_type_detailed %in% c(targets.fibro, "Myeloid: Macro-C1QC", "Myeloid: Macro-SPP1"))
seurat_obj <- subset(seurat_obj, !(SPECIMEN_TYPE == "N" & cell_type_detailed %in% c("Myeloid: Macro-SPP1", "Myeloid: Macro-C1QC")))

# need to make the idents the spp1/ c1qc labels then just CD8 and CD4
seurat_obj@meta.data$cclabels <- ifelse(
  seurat_obj@meta.data$clTopLevel == "Myeloid",
  seurat_obj@meta.data$cell_type_detailed,
  seurat_obj@meta.data$clMidwayPr)

# Set cell type labels as the active identity for NicheNet
Idents(seurat_obj) <- seurat_obj$cclabels

# Verify cell types and conditions
table(seurat_obj$cclabels)
table(seurat_obj$SPECIMEN_TYPE)

print("seurat object created")


##############################
# Define nichenet function
##############################

# nichnet function must be defined to loop over all lr_networks
run_nichenet <- function(lr_network_use,
                         # the name of the network used for later classification
                         network_name,
                         # the up/downregulation cutoff
                         logfc_cutoff = 0.25,
                         # the default number of ligands used for further analysis
                         top_n = 31, 
                         # the default number of gene targets displayed
                         ngene_targets_use = 300) {
  
  # define sender and reciever based on the lr_network name
  receiver = "Fibro"
  
  if (network_name == "lr_spp1_Fibro") {
    sender <- "Myeloid: Macro-SPP1"
  } else if (network_name == "lr_c1qc_Fibro") {
    sender <- "Myeloid: Macro-C1QC"
  } else if (network_name == "lr_both_Fibro") {
    sender <- c("Myeloid: Macro-SPP1",
                "Myeloid: Macro-C1QC")
  } else {
    stop("Unknown network_name: ", network_name)
  }
  
  # add the sender to sender cell types
  sender_celltypes <- sender
  
  # Define the condition comparison
  condition_oi        <- "T"  # condition of interest
  condition_reference <- "N"         # reference condition
  
  
  ##############################
  # Define the gene set of interest
  ##############################
  
  # Extract the count matrix and cell identity vector once for efficiency.
  expr_mat    <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  cell_idents <- Idents(seurat_obj)
  
  # Identify expressed genes in the receiver population
  # A gene is "expressed" if detected in >=10% of receiver cells
  expressed_genes_receiver <- get_expressed_genes(
    celltype_oi    = receiver,
    object         = expr_mat,
    celltype_annot = cell_idents,
    pct            = 0.10)
  
  # Define background genes: expressed receiver genes present in the NicheNet model
  background_expressed_genes <- expressed_genes_receiver %>%
    .[. %in% rownames(ligand_target_matrix)]
  
  # Identify expressed genes in each sender cell type (reuse the same matrix/idents)
  list_expressed_genes_sender <- sender_celltypes %>%
    lapply(function(ct) get_expressed_genes(
      celltype_oi    = ct,
      object         = expr_mat,
      celltype_annot = cell_idents,
      pct            = 0.10))
  
  expressed_genes_sender <- list_expressed_genes_sender %>%
    unlist() %>%
    unique()
  
  ##############################
  # Run ligand activity analysis
  ##############################
  
  # Subset to receiver cell type
  seurat_receiver <- subset(seurat_obj, idents = receiver)
  
  # Set condition as identity for DE testing within receiver cells
  Idents(seurat_receiver) <- seurat_receiver$SPECIMEN_TYPE
  
  # Find DE genes: Post_Treatment vs Healthy in CD4+ T cells
  DE_table_receiver <- FindMarkers(
    object   = seurat_receiver,
    ident.1  = condition_oi,
    ident.2  = condition_reference,
    min.pct  = 0.10)
  
  # Filter to significant, biologically meaningful DE genes
  # and restrict to genes present in the NicheNet ligand-target model
  if (logfc_cutoff > 0) {
    geneset_oi <- DE_table_receiver %>%
      filter(p_val_adj <= 0.05,
             avg_log2FC >= logfc_cutoff) %>%
      rownames()
  } else {
    geneset_oi <- DE_table_receiver %>%
      filter(p_val_adj <= 0.05,
             avg_log2FC <= logfc_cutoff) %>%
      rownames()
  } %>%  # upregulated in Post_Treatment
    rownames() 
  
  geneset_oi <- geneset_oi %>%
    intersect(rownames(ligand_target_matrix))
  length(geneset_oi)  # how many genes passed filtering
  
  
  ##############################
  # Downstream analysis
  ##############################
  
  # Identify all receptors expressed by the receiver population
  expressed_receptors <- intersect(
    lr_network_use$to,
    expressed_genes_receiver
  )
  
  # --- Sender-focused approach ---
  expressed_ligands_sender <- intersect(
    lr_network_use$from,
    expressed_genes_sender
  )
  
  # Potential ligands: expressed in sender AND have expressed receptor in receiver
  potential_ligands_focused <- lr_network_use %>%
    filter(from %in% expressed_ligands_sender & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()
  
  
  ##############################
  # Run Ligand Activity Analysis
  ##############################
  
  print(potential_ligands_focused)
  # line to remove if omnipath???
  potential_ligands_focused <- intersect(
    potential_ligands_focused,
    colnames(ligand_target_matrix))
  print(potential_ligands_focused)
  
  if (is.null(potential_ligands_focused) || length(potential_ligands_focused) == 0) {
    return(NULL)
  }
  
  
  # Run ligand activity analysis — sender-focused
  ligand_activities_focused <- predict_ligand_activities(
    geneset                  = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix     = ligand_target_matrix,
    potential_ligands        = potential_ligands_focused
  )
  
  
  # Sort by corrected AUPR and assign rank using row_number() on the sorted result.
  
  ligand_activities_focused <- ligand_activities_focused %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank = row_number())
  
  
  ##############################
  # Select Ligands for downstream analysis
  ##############################
  
  best_upstream_ligands_focused <- ligand_activities_focused %>%
    top_n(top_n, aupr_corrected) %>%
    pull(test_ligand) %>%
    unique()
  
  # Combine (union of both approaches for comprehensive follow-up)
  best_upstream_ligands <- best_upstream_ligands_focused
  
  #length(best_upstream_ligands)
  # remove TNF as it dominates results
  # best_upstream_ligands <- best_upstream_ligands[! best_upstream_ligands %in% c("TNFSF12")]
  # length(best_upstream_ligands)
  
  
  ##############################
  # Run Ligand Activity Analysis
  ##############################
  
  ligand_activities_df <- ligand_activities_focused
  
  top_ligands <- ligand_activities_df %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = top_n) %>%
    pull(test_ligand)
  
  # Helper to build a bar plot for one set of ligand activities
  plot_ligand_activity <- function(ligand_activities_df, top_ligands, subtitle) {
    ligand_activities_df %>%
      filter(test_ligand %in% top_ligands) %>%
      mutate(test_ligand = factor(test_ligand, levels = rev(top_ligands))) %>%
      ggplot(aes(x = test_ligand, y = aupr_corrected)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title    = "NicheNet Ligand Activity Scores",
        subtitle = subtitle,
        x        = "Ligand",
        y        = "Corrected AUPR"
      ) +
      theme_classic()
  }
  
  
  p_focused <- plot_ligand_activity(
    ligand_activities_focused,
    best_upstream_ligands_focused,
    "Top ligands — sender-focused")
  
  ##############################
  # Predict Target Genes of Top Ligands
  ##############################
  
  # Define a regulatory potential cutoff
  # get_weighted_ligand_target_links extracts the top target genes
  # for each ligand above a given threshold
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset              = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n                    = ngene_targets_use) %>%  
    bind_rows() %>%
    drop_na()
  
  # Convert to a wide matrix for heatmap visualization.
  # cutoff = 1/3 retains only weights above the 33rd percentile;
  # prepare_ligand_target_visualization handles this filtering internally.
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df     = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff               = 0.33333
  )
  
  # Order ligands by activity rank for consistent heatmap layout
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
    rev()
  order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links))
  
  # Generate the ligand-target heatmap
  if (length(order_targets) == 0 || length(order_ligands) == 0) {
    warning("No ligand-target interactions for this comparison.")
    vis_ligand_target <- NULL
    vis_ligand_receptor_network <- NULL
  } else {
    vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands, drop = FALSE]
    
    ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
      best_upstream_ligands,
      expressed_receptors,
      lr_network_use,
      weighted_networks$lr_sig
    )
    
    # SAFE CHECK
    if (is.null(ligand_receptor_links_df) ||
        nrow(ligand_receptor_links_df) == 0 ||
        length(unique(ligand_receptor_links_df$from)) < 2 ||
        length(unique(ligand_receptor_links_df$to)) < 2) {
      
      vis_ligand_receptor_network <- NULL
      
    } else {
      
      vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
        ligand_receptor_links_df,
        best_upstream_ligands,
        order_hclust = "both"
      )
      
    }
  }
  
  return(list(
    # ligand_activity = ligand_activities_focused,
    # ligand_heatmap = p_ligand_target_network,
    # ligand_barplot = p_focused,
    # ligands = best_upstream_ligands,
    # geneset = geneset_oi,
    gene_expression = vis_ligand_target,
    receptor_expression = vis_ligand_receptor_network
  ))
}


##############################
# run nichenet function for all 6 dataframes
##############################

# find up and down regulated ligands
logfcs <- c(0.25, -0.25)

# create a list for results
results <- list()

# loop over differen networks
for (net_name in names(lr_networks_Fibro)) {
  net <- lr_networks_Fibro[[net_name]]
  
  # for up and down regulated genes
  for (lfc in c(0.25, -0.25)) {
    # add up or down to the df name
    run_name <- paste0(net_name, ifelse(lfc > 0, "_up", "_down"))
    
    cat("Running:", run_name, "\n")
    
    
    results[[run_name]] <- run_nichenet(
      lr_network_use = net,
      network_name = net_name,
      logfc_cutoff = lfc,
      ngene_targets_use = 300, # top 200 genes for later figure
      top_n = 31 # use all ligands from cellchat
    )
  }
}

saveRDS(results, file = "../data/nichenet_results_fibro.rds")

print("nichenet object saved. FINITO")


##############################
# End of HPC
##############################
###########################################################
# Part 2
###########################################################
##############################
# Load data and libraries
##############################

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(patchwork)

pbmc <- readRDS("../data/pbmc_final.rds")
lr_networks_Fibro <- readRDS("../data/lr_networks_Fibro.rds")

# load in the results from nichenet
results <- readRDS("../data/nichenet_results_fibro.rds")


##############################
# visualise results- differences between priority genes
##############################

# create a list of all up/downregulated data to combine into a single figure
View(results)
# upregulated 
mat_list <- list(
  lr_spp1_Fibro = results$lr_spp1_Fibro_up$gene_expression,
  lr_both_Fibro = results$lr_both_Fibro_up$gene_expression,
  lr_c1qc_Fibro = results$lr_c1qc_Fibro_up$gene_expression
)

# downregulated
# fix for ones with no columns
mat_list <- list(
  lr_spp1_Fibro = results$lr_spp1_Fibro_down$gene_expression,
  lr_both_Fibro = results$lr_both_Fibro_down$gene_expression,
  lr_c1qc_Fibro = results$lr_c1qc_Fibro_down$gene_expression
)

# create a single matrix from the list of matrices
mat_list <- lapply(mat_list, as.matrix)
# get all unique genes across matrices
all_genes <- Reduce(union, lapply(mat_list, rownames))
# align them into the same gene set
mat_list <- lapply(mat_list, function(m) {
  # find missing genes
  missing <- setdiff(all_genes, rownames(m))
  # set missing genes to 0
  if (length(missing) > 0) {
    pad <- matrix(0,
                  nrow = length(missing),
                  ncol = ncol(m),
                  dimnames = list(missing, colnames(m)))
    # add missing genes as 0 expression rows
    m <- rbind(m, pad)
  }
  # reorder rows consistently
  m[all_genes, , drop = FALSE]
})

# add suffix to column names
mat_list <- Map(function(m, i) {colnames(m) <- paste0(colnames(m), "_", i)
m}, mat_list, seq_along(mat_list))

# combine all matrices into one big matrix
mega_heatmap_up <- do.call(cbind, mat_list)

# too many genes are present
# remove genes that have a lower regulatory potential
nrow(mega_heatmap_up)
keep <- apply(mega_heatmap_up, 1, function(x) max(x) - min(x) >= 0.15)
mega_heatmap_up <- mega_heatmap_up[keep, , drop = FALSE]
nrow(mega_heatmap_up)


# visualise the mega matrix ligand target network
b7 <- mega_heatmap_up %>%
  make_heatmap_ggplot(
    y_name = paste("Target genes in "),
    x_name = paste("prioritized ligands "),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) +
  coord_flip()
b7


##############################
# identify genes expressed through specific types of interaction
##############################
lr_networks_Fibro <- readRDS("../data/lr_networks_Fibro.rds")

# load in the results from nichenet
results <- readRDS("../data/nichenet_results_fibro.rds")

# load appropriate library
library(ComplexUpset)

geneset_lists <- list()

# loop version
# loop over the different interaction conditions
for (nm in names(results)) {
  # include gene expression
  mat <- results[[nm]]$gene_expression
  # skip if NULL or empty
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    next
  }
  # ensure matrix
  mat <- as.matrix(mat)
  # keep genes with any expression signal
  expressed_genes <- rownames(mat)[rowSums(!is.na(mat) & mat > 0) > 0.03]
  # add to geneset lists
  geneset_lists[[nm]] <- expressed_genes
}


# sort
all_geneset <- sort(unique(unlist(geneset_lists)))
# covert to dataframe
geneset_df <- data.frame(geneset = all_geneset)

# loop over dataframes
for (nm in names(geneset_lists)) {
  geneset_df[[nm]] <- all_geneset %in% geneset_lists[[nm]]
}

# add row names
rownames(geneset_df) <- geneset_df$geneset
geneset_df$geneset <- NULL

# plot
ComplexUpset::upset(
  geneset_df,
  intersect = colnames(geneset_df),
  sort_intersections_by = "degree",
  sort_sets = "descending",
  keep_empty_groups = TRUE,
  base_annotations = list("Intersection size" = intersection_size()))

# get a list of each gene in each group
geneset_logical <- geneset_df > 0
group_key <- apply(
  geneset_logical,
  1,
  function(x) paste(colnames(geneset_logical)[x], collapse = "///")
)
groups <- split(rownames(geneset_logical), group_key)
group_table <- data.frame(
  group = names(groups),
  genes = sapply(groups, paste, collapse = ", ")
)

# display groups
groups
# saveRDS(groups, file = "../data/nichenet_gene_groups_Fibro.rds")

##############################
# Make Figure
##############################

# Figure 7
ggdraw() +
  draw_plot(a7, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(b7, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B"), size = 15, 
                  x = c(0, 0), 
                  y = c(1, .5))
# save figure
ggsave("../figures/Figure_7_nichenet_general_Fibro.jpg", width = 40, height = 30, units = c("cm"), dpi = 300)


##############################
# Make Figure
##############################

# load in the results from nichenet
results <- readRDS("../data/nichenet_results_fibro.rds")
groups <- readRDS("../data/nichenet_gene_groups_Fibro.rds")

# include only gene which are exclusive to SPP1
spp1_Fibro_up_only_conditions <- c(groups$lr_spp1_Fibro_up)

heatmap_spp1_Fibro_up <- results$lr_spp1_Fibro_up$gene_expression[, -c(1, 3, 4, 5, 6, 8, 9, 10)]
heatmap_spp1_Fibro_up <- heatmap_spp1_Fibro_up[rownames(heatmap_spp1_Fibro_up) %in% spp1_Fibro_up_only_conditions, ]


nrow(heatmap_spp1_Fibro_up)
keep <- apply(heatmap_spp1_Fibro_up, 1, function(x) max(x) >= 0.04)
heatmap_spp1_Fibro_up <- heatmap_spp1_Fibro_up[keep, , drop = FALSE]
nrow(heatmap_spp1_Fibro_up)


e8 <- heatmap_spp1_Fibro_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
e8

# C1QC
# include only gene which are exclusive to c1qc
c1qc_Fibro_up_only_conditions <- c(groups$lr_c1qc_Fibro_up)

heatmap_c1qc_Fibro_up <- results$lr_c1qc_Fibro_up$gene_expression# [, -c(2)]
heatmap_c1qc_Fibro_up <- heatmap_c1qc_Fibro_up[rownames(heatmap_c1qc_Fibro_up) %in% c1qc_Fibro_up_only_conditions, ]


nrow(heatmap_c1qc_Fibro_up)
keep <- apply(heatmap_c1qc_Fibro_up, 1, function(x) max(x) >= 0.015)
heatmap_c1qc_Fibro_up <- heatmap_c1qc_Fibro_up[keep, , drop = FALSE]
nrow(heatmap_c1qc_Fibro_up)

f8 <- heatmap_c1qc_Fibro_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
f8



ggdraw() +
  draw_plot(a8, x = 0, y = .5, width = .33, height = .5) +
  draw_plot(b8, x = .33, y = .5, width = .33, height = .5) +
  draw_plot(c8, x = .66, y = .5, width = .33, height = .5) +
  draw_plot(d8, x = 0, y = 0, width = .33, height = .5) +
  draw_plot(e8, x = .33, y = 0, width = .33, height = .5) +
  draw_plot(f8, x = .66, y = 0, width = .33, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 15, 
                  x = c(0, .33, .66, 0, .33, .66), 
                  y = c(1, 1, 1, .5, .5, .5))

ggsave("../figures/Figure_8_fibroblast.jpg", width = 30, height = 20, units = c("cm"), dpi = 300)


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

##### DATASET SCREENING AND SUBSET IDENTIFICATION
# Step 1: Preproccessing and preliminary analysis (following Xie et al)
# Step 2: Clustering and analysis of CES1 expression


##############################
# load in the Pelka Dataset
##############################

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load in dataset into a dataframe
pbmc_data <- Read10X_h5("../data/GSE178341_crc10x_full_c295v4_submit.h5")

# Initialise the Seurat object with the raw non-normalized data
pbmc <- CreateSeuratObject(counts=pbmc_data, 
                           project = "pelka_sc", # project name 
                           min.cells = 10, # include features detected with at least 10 cells
                           min.features = 500) #include cells where at least this many features are detected


##############################
# Preprocessing
##############################

# create Mitochondrial gene column
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# visualise QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# one line filter all cells with below 6000 feature counts, above 50, and less than 20 mitochondrial
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)


##############################
# Normalizing the data
##############################

# normalize feature expression measurements for each cell by total expression
# use default values (scale factor 10000, lognormalize function)
pbmc <- NormalizeData(pbmc)


##############################
# feature selection
##############################

# the top 2000 are identified based on mean and dispersion
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# save progress
# saveRDS(pbmc, file = "../data/pbmc_post_feature_selection.rds")
# pbmc <- readRDS("../data/pbmc_post_feature_selection.rds")


##############################
# Data Scaling 
##############################

# shifts and scales gene expression so the mean is 0 and the variance is 1. 
# gives equal weight to downstream analysis
# stored in pbmc[["RNA"]]$scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt"))


##############################
# Linear dimensional reduction
##############################

# principal component analysis
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# visualise PCAs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()

# dim heatmaps allow for easy exploration of primary sources of heterogeneity
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

##############################
# Determine dimentionality
##############################

# to overcome extensive technical noise
ElbowPlot(pbmc, reduction = "pca",ndims = 50)
# most of the dimensionality is in the first 30 PCs


##############################
# cluster the cells
##############################

# the distance metric that drives clustering (based on PCs) stays the same
dim.use <- 1:30
pbmc <- FindNeighbors(pbmc, dims = dim.use)
pbmc <- FindClusters(pbmc, resolution = 0.1)
# head(Idents(pbmc), 5)

##############################
# Run non-linear dimensional reduction
##############################

# visualisation of the dimensionality
pbmc <- RunUMAP(pbmc, dims = dim.use)
DimPlot(pbmc, reduction = "umap", label = T)
# 13 clusters identified

# save progress
# saveRDS(pbmc, file = "../data/pbmc_post_UMAP_30PCs.rds")
# pbmc <- readRDS("../data/pbmc_pbmc_post_UMAP_30PCs.rds")


##############################
# Finding differentially expressed features (cluster biomarkers)
##############################

# Nine major cell types were identified in these datasets, 
# including T cells (CD3D), natural killer cells (KLRF1), B cells (MS4A1),
# plasma cells (MZB1), myeloid cells (LYZ), mast cells (TPSAB1), 
# neutrophils (FCGR3B), epithelial cells (EPCAM), and stromal cells (DCN)

# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# save progress
# saveRDS(pbmc.markers, file = "../data/pbmc.markers.rds")
# pbmc.markers <- readRDS("../data/pbmc.markers.rds")

# identification of marker genes
# cluster 0: (CD3D) TRBC2- T cell
head(pbmc.markers[pbmc.markers$cluster == 0, ], 5)
# cluster 1: KRT18 KRT8 (EPCAM)- Epithelial
head(pbmc.markers[pbmc.markers$cluster == 1, ], 5)
# cluster 2: (MZB1) DERL3 (CD79A)- Plasma
head(pbmc.markers[pbmc.markers$cluster == 2, ], 5)
# cluster 3: AIF1 TYR0BP (LYZ)- Macrophage
head(pbmc.markers[pbmc.markers$cluster == 3, ], 5)
# cluster 4: (MS4A1) CD79A- B cell
head(pbmc.markers[pbmc.markers$cluster == 4, ], 5)
# cluster 5: PIGR LGALS4 (FABP1)- Enterocyte
head(pbmc.markers[pbmc.markers$cluster == 5, ], 10)
# cluster 6: CALD1 IGFBP7 (COL6A2)- Fibroblasts
head(pbmc.markers[pbmc.markers$cluster == 6, ], 5)
# cluster 7: (GNG11) RAMP2- Endothelial cell
head(pbmc.markers[pbmc.markers$cluster == 7, ], 5)
# cluster 8: IGF2 SMIM22 (EPCAM)- Epithelial
head(pbmc.markers[pbmc.markers$cluster == 8, ], 5)
# cluster 9: TNNC2 CXCL14 Malignant epithelial
head(pbmc.markers[pbmc.markers$cluster == 9, ], 10)
# cluster 10: (TPSAB1) CPA3- Mast Cell
head(pbmc.markers[pbmc.markers$cluster == 10, ], 5)
# cluster 11: AZGP1 WFDC2- Tumour epithelial
head(pbmc.markers[pbmc.markers$cluster == 11, ], 10)
# cluster 12: (LAMP3) CCR7- Dendritic cells
head(pbmc.markers[pbmc.markers$cluster == 12, ], 5)
# cluster 13: LILRA4 GZMB (PLAC8)- Dendritic cells
head(pbmc.markers[pbmc.markers$cluster == 13, ], 10)


# use violin plot to visualise identified marker expression
VlnPlot(pbmc, features = c("CES1"), raster = TRUE, split.by = "celltype", log = T, pt.size = 0)

features= c("EPCAM","CD3D","MZB1","LYZ","MS4A1","FABP1","COL6A2","GNG11",
            "TPSAB1","LAMP3","PLAC8")
DotPlot(pbmc, features = unique(features)) + RotatedAxis()

# feature plots 
FeaturePlot(pbmc, features = features,
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")

# build a cluster tree
pbmc <- BuildClusterTree(
  pbmc,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
PlotClusterTree(pbmc)

# label clusters
cluster2celltype <- c( "0"="T",
                       "1"="Epithelial", 
                       "2"="Plasma", 
                       "3"= "Myeloid", 
                       "4"= "B", 
                       "5"= "Epithelial",
                       "6"= "Stromal", 
                       "7"= "Endothelial", 
                       "8"= "Epithelial",
                       "9"= "Epithelial",
                       "10"="Mast",
                       "11"="Epithelial", 
                       "12"="Dendritic", 
                       "13"="Dendritic"
                       
)

# add cluster labels to pbmc
pbmc[['celltype']] = unname(cluster2celltype[pbmc@meta.data$seurat_clusters])
Idents(pbmc) <- "celltype"
View(pbmc@meta.data)


##############################
# Meta analysis
##############################

# import meta data
meta_data <- read.csv("../data/GSE178341_crc10x_full_c295v4_submit_metatables_.csv")
# check the col names align
pbmc[[]]
head(meta_data)

# reformat row names 
rownames(meta_data) <- meta_data$cellID
meta_data$cellID <- NULL
meta_data <- meta_data[colnames(pbmc), ]

# then add the metadata
pbmc <- AddMetaData(pbmc, metadata = meta_data)

#meta data check
describe(pbmc@meta.data)

# visualisation
#patient
DimPlot(object = pbmc, group.by="orig.ident", reduction='umap')
# some patinets have unique epithelial mutations

# methylatedor nonmethylated (inherited or not) MMRd
DimPlot(object = pbmc, group.by="MMRMLH1Tumor", reduction='umap')

# tissue site
DimPlot(object = pbmc, group.by="TissueSiteSimple", reduction='umap')
dev.off()

#Histological grade
DimPlot(object = pbmc, group.by="HistologicGradeSimple", reduction='umap')

#Sex
DimPlot(object = pbmc, group.by="Sex", reduction='umap')

# Tumour stage
DimPlot(object = pbmc, group.by="TumorStage", reduction='umap')

# add extra meta data values
# add ces1 expression to the metadata
ces1 <- FetchData(pbmc, vars = "CES1")
pbmc$CES1_expression <- ces1$CES1

# add CES1 expression TRUE/FALSE column to metadata
pbmc_macrophage@meta.data$CES1_isExpressed <- FetchData(pbmc_macrophage, vars = "CES1")[,1] > 0


##############################
# Cluster analysis
##############################

cluster_data <- read.csv("../data/crc10x_full_c295v4_submit_cluster.csv")
# check the col names align
pbmc[[]]
head(cluster_data)
View(cluster_data)

# reformat row names 
rownames(cluster_data) <- cluster_data$sampleID
cluster_data$sampleID <- NULL

# check row names match
all(colnames(pbmc) %in% rownames(cluster_data))
all(rownames(cluster_data) == colnames(pbmc))

# check the row names are identical
head(colnames(pbmc))
head(rownames(cluster_data))

# reorder the cluster colnames
cluster_data <- cluster_data[colnames(pbmc), ]

# then add the cluster data 
pbmc <- AddMetaData(pbmc, metadata = cluster_data)

#cluster data check
View(pbmc@meta.data)
rm(pbmc@meta.data$sampleID)


##############################
# Make Figures
##############################

# save featureplot of CES1
pdf(file = "../figures/pbmc_CES1_featureplot.pdf", width = 10, height = 6)
FeaturePlot(pbmc, features = c("CES1"), label = T,
            reduction = "umap", pt.size = 1, cols = cont_2) +
  ggtitle(NULL)
dev.off()

# save dotplot of CES1
pdf(file = "../figures/pbmc_CES1_dotplot.pdf", width = 5, height = 6)
c1 <- DotPlot(pbmc, features = 'CES1', cols = cont_2) + RotatedAxis() + NoLegend() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
c1
dev.off()

# fibroblasts have high levels of CES1
# check with just immune cells
pbmc_key_cells <- subset(pbmc, celltype %in% c("T","Plasma", "Mast","B", "Dendritic", "Myeloid"))
# save dotplot of CES1 with only key cells
pdf(file = "../figures/pbmc_key_cells_CES1_dotplot.pdf", width = 5, height = 6)
c2 <- DotPlot(pbmc_key_cells, features = 'CES1', cols = cont_2) + RotatedAxis() + NoLegend() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
c2
dev.off()

# save umap with my labelled cells
pdf(file = "../figures/pbmc_labelled_umap.pdf", width = 10, height = 6)
DimPlot(pbmc, reduction = "umap", group.by = 'celltype',label = T, cols = disc_10) +
  ggtitle(NULL)
dev.off()

# make figure from previous clusters
pdf(file = "../figures/pbmc_labelled_umap_clTopLevel.pdf", width = 10, height = 6)
a1 <- DimPlot(pbmc, reduction = "umap", group.by = 'clTopLevel' ,label = F, cols = disc_10) +
  ggtitle(NULL) + NoLegend() + clean_theme() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), plot.margin=grid::unit(c(0,0,0,0), "mm"))
a1
dev.off()

# save umap with Normal/Tumour labels
pdf(file = "../figures/pbmc_labelled_umap_specimen_type.pdf", width = 10, height = 6)
DimPlot(pbmc, reduction = "umap", group.by = 'SPECIMEN_TYPE' ,label = F, cols = disc_10) +
  ggtitle(NULL)
dev.off()


##############################
# Final Saves
##############################

# these are the necessary files with all metadata and processing

# saveRDS(pbmc, file = "../data/pbmc_final.rds")
# pbmc <- readRDS("../data/pbmc_final.rds")

# saveRDS(pbmc_key_cells, file = "../data/pbmc_key_cells.rds")
# pbmc <- readRDS("../data/pbmc_key_cells.rds")

# saveRDS(pbmc.markers, file = "../data/pbmc.markers.rds")
# pbmc.markers <- readRDS("../data/pbmc.markers.rds")

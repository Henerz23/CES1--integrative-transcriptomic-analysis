###########################################################
# DATASET SCREENING AND SUBSET IDENTIFICATION
###########################################################

# Part 1: The Pelka dataset screened for appropriate granularity
# Part 2: Cell subtypes analysed for CES1 expression


##############################
# load libraries
##############################

library(Seurat)
library(patchwork)
library(Hmisc)
library(RColorBrewer)
library(tidyverse)
library(cowplot)


# Colours!
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")


###########################################################
# Part 1
###########################################################
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

# one line filter all cells with below 6000 feature counts, above 500, and less than 20 mitochondrial
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
# all.genes <- rownames(pbmc)
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
# Determine dimensionality
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

##############################
# Run non-linear dimensional reduction
##############################

# visualisation of the dimensionality
pbmc <- RunUMAP(pbmc, dims = dim.use)
DimPlot(pbmc, reduction = "umap", label = T)
table(Idents(pbmc))


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
# some patients have unique epithelial mutations

# methylated or non-methylated (inherited or not) MMRd
DimPlot(object = pbmc, group.by="MMRMLH1Tumor", reduction='umap')

# tissue site
DimPlot(object = pbmc, group.by="TissueSiteSimple", reduction='umap')

#Histological grade
DimPlot(object = pbmc, group.by="HistologicGradeSimple", reduction='umap')

#Sex
DimPlot(object = pbmc, group.by="Sex", reduction='umap')

# Tumour stage
DimPlot(object = pbmc, group.by="TumorStage", reduction='umap')

# add extra meta data values
ces1 <- FetchData(pbmc, vars = "CES1")
# add ces1 expression to the metadata
pbmc@meta.data$CES1_expression <- ces1$CES1
# add CES1 expression TRUE/FALSE column to metadata
pbmc@meta.data$CES1_isExpressed <- ces1$CES1 > 0


##############################
# Cluster analysis
##############################

cluster_data <- read.csv("../data/crc10x_full_c295v4_submit_cluster.csv")
# check the col names align
head(pbmc[[]])
# head(cluster_data)

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

###########################################################
# Part 2
###########################################################
##############################
# Make Figures
##############################

# make the top level annotations the identities
Idents(pbmc) <- "clTopLevel"

# save featureplot of CES1
c1 <- FeaturePlot(pbmc, features = c("CES1"), label = F,
            reduction = "umap", cols = cont_2, order = T, pt.size = 0.1)+
            NoLegend() +
            ggtitle(NULL)

# add figure legend to clearly distinguish between cells expressing and not expressing CES1
c1 <- c1 +
  annotate("point", x = -15, y = -15, colour = cont_2[2], size = 4) +
  annotate("text",  x = -14, y = -15, label = "Expressing CES1", hjust = 0, size = 4) +
  annotate("point", x = -15, y = -17, colour = cont_2[1], size = 4) +
  annotate("text",  x = -14, y = -17, label = "Not Expressing CES1", hjust = 0, size = 4) 
c1

# save dotplot of CES1
d1 <- DotPlot(pbmc, features = 'CES1', cols = cont_2) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
d1
d1$data

# Stromal have high levels of CES1
# check with just immune cells
pbmc_immune_cells <- subset(pbmc, clTopLevel %in% c("TNKILC","Plasma", "Mast","B", "Dendritic", "Myeloid"))
# save dotplot of CES1 with only immune cells
DotPlot(pbmc_immune_cells, features = 'CES1', cols = cont_2) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# make figure from previous clusters
a1 <- DimPlot(pbmc, reduction = "umap", group.by = 'clTopLevel' ,label = F, cols = disc_10) +
  ggtitle(NULL)
a1

# save umap with Normal/Tumour labels
b1 <- DimPlot(pbmc, reduction = "umap", group.by = 'SPECIMEN_TYPE' ,label = F, cols = disc_10) +
  ggtitle(NULL)
b1

# Figure 1
ggdraw() +
  draw_plot(a1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(b1, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(c1, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(d1, x = .5, y = 0, width = .3, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15, 
                  x = c(0, .5, 0, .5), 
                  y = c(1, 1, .5, .5))
# save
ggsave("../figures/Figure_1_general.jpg", width = 30, height = 20, units = c("cm"), dpi = 300)

##############################
# Final Saves
##############################

# saveRDS(pbmc, file = "../data/pbmc_final.rds")
# pbmc <- readRDS("../data/pbmc_final.rds")

# saveRDS(pbmc.markers, file = "../data/pbmc.markers.rds")
# pbmc.markers <- readRDS("../data/pbmc.markers.rds")

# saveRDS(pbmc_immune_cells, file = "../data/pbmc_immune_cells.rds")
# pbmc_immune_cells <- readRDS("../data/pbmc_immune_cells.rds")

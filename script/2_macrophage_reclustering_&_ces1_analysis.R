###########################################################
# SPP1/C1QC MACROPHAGE RECLUSTERING AND CES1 ANALYSIS
###########################################################

# Part 1: As Pelka has insufficient macrophage labels, they will be reclustered
# Part 2: New macrophage subsets will be analysed for CES1 levels


##############################
# Load Libraries
##############################

library(Seurat)
library(patchwork)
library(ggrastr)
library(Hmisc)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(cowplot)
library(randomcoloR)
library(scales)

# Establish some colours
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")
# set the random seed as 42 for consistent colour generation
set.seed(42)

# as the macrophage class has been identified to have high expression of CES1, 
# they will be further analysed

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###########################################################
# Part 1
###########################################################
##############################
# Myeloid subset analysis
##############################

# load in pbmc
# pbmc <- readRDS("../data/pbmc_final.rds")

# subset the myeloid group
# done using the previous annotations to avoid spending days clustering the wrong thing before realising it isnt going to work for cellchat
pbmc_myeloid <- subset(pbmc, clTopLevel %in% c("Myeloid"))
dim(pbmc_myeloid)
table(pbmc_myeloid@meta.data$celltype)


##############################
# Data Preprocessing
##############################

# normalization
pbmc_myeloid <- NormalizeData(pbmc_myeloid, normalization.method = "LogNormalize", scale.factor = 10000) 
pbmc_myeloid <- FindVariableFeatures(pbmc_myeloid, selection.method = 'vst', nfeatures = 2000)

# scaling
pbmc_myeloid <- ScaleData(pbmc_myeloid, vars.to.regress = "percent.mt")
pbmc_myeloid <- RunPCA(pbmc_myeloid, features = VariableFeatures(object = pbmc_myeloid)) 

# dim use identification
ElbowPlot(pbmc_myeloid, reduction = "pca", ndims = 50)

# dimensional reduction
dim.use<-1:30
pbmc_myeloid <- FindNeighbors(pbmc_myeloid, dims = dim.use)
pbmc_myeloid <- FindClusters(pbmc_myeloid, resolution = 0.8 )
pbmc_myeloid <- RunUMAP(pbmc_myeloid, dims = dim.use)

# plot clusters
DimPlot(pbmc_myeloid, reduction = 'umap',label=T)


##############################
# Finding differentially expressed genes
##############################

# find the markers for the macrophage clusters
pbmc_myeloid.markers <- FindAllMarkers(pbmc_myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_myeloid.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# saveRDS(pbmc_myeloid.markers, file = "../data/pbmc_myeloid.markers_PC30.rds")
# pbmc_myeloid.markers <- readRDS("../data/pbmc_myeloid.markers_PC30.rds")


##########################################################################################
# Clustering 
##########################################################################################

# 3 steps:
# Plan A: check Pelka dataset for useful clustering annotations
# Plan B: check Xie et al for useful clustering annotations
# Plan C: use Xie et al markers to recluster manually


##############################
# Clustering using Pelka clusters
##############################

# analysis using clustering information
a2 <- DimPlot(pbmc_myeloid, reduction = 'umap',label=F, 
              group.by = 'cl295v11SubFull', cols = distinctColorPalette(10)) +
  ggtitle(NULL) 
a2
# no annotation of SPP1 and C1QC+ Macrophages


##############################
# Checking Xie at al markers
##############################

# no clustering data is available from Xie et al. 
# marker genes were therefore used to identify clusters using a multistep analysis

# utilise dotplots and their markers to identify subgroups
# pDC-LILRA4
pDC_LILRA4_features = c("LILRA4", "GZMB", "PTGDS", "CLIC3", "IRF7", "PLAC8", "TSPAN13", "C12orf75")
DotPlot(pbmc_myeloid, features = pDC_LILRA4_features, cols = cont_2) + RotatedAxis()
# 15

# cDC1-CLEC9A
cDC1_CLEC9A_features = c("CLEC9A", "DNASE1L3", "IDO1", "TACSTD2", "ASB2", "XCR1", "HLA-DRs")
DotPlot(pbmc_myeloid, features = cDC1_CLEC9A_features, cols = cont_2) + RotatedAxis()
# 19

# cDC2-CD1C
cDC2_CD1C_features = c("CD1C", "CD1E", "FCER1A", "CLEC10A", "AREG", "PPA1", "LGALS2", "HLA-DRs")
DotPlot(pbmc_myeloid, features = cDC2_CD1C_features, cols = cont_2) + RotatedAxis()
# 9

# cDC3-LAMP3
cDC3_LAMP3_features = c("LAMP3", "CCL17", "CCL19", "CCL22","BIRC3", "CCR7", "FSCN1", "TXN", "IDO1")
DotPlot(pbmc_myeloid, features = cDC3_LAMP3_features, cols = cont_2) + RotatedAxis()
# 13

# cDC4-LTB
cDC4_LTB_features = c("LTB", "CDH17", "PRDM16", "TLR10", "SFTPD", "AFF3", "SUSD3", "IGFBP6")
DotPlot(pbmc_myeloid, features = cDC4_LTB_features, cols = cont_2) + RotatedAxis()
#

# Macro-C1QC
Macro_C1QC_features = c("C1QA", "C1QB", "C1QC", "SEPP1", "SLC40A1", "APOE", "RNASE1", "CTSC")
DotPlot(pbmc_myeloid, features = Macro_C1QC_features, cols = cont_2) + RotatedAxis()
# 0, 12 10 6 14 18

# Macro-FCN1
Macro_FCN1_features = c("FCN1", "S100A8", "S100A9", "S100A12","VCAN", "EREG", "THBS1", "TIMP1", "IL1B", "VEGFA")
DotPlot(pbmc_myeloid, features = Macro_FCN1_features, cols = cont_2) + RotatedAxis()
# 5 1 2 16 4 

# Macro-MKI67
Macro_MKI67_features = c("MKI67", "HIST1H4C", "STMN1", "HMGB2", "PCLAF", "TYMS", "TUBB", "TK1")
DotPlot(pbmc_myeloid, features = Macro_MKI67_features, cols = cont_2) + RotatedAxis()
# 8

# Macro-SPP1
Macro_SPP1_features = c("SPP1", "CCL2", "CLEC5A", "CSTB", "APOC1", "CXCL8", "SDC2", "FN1", "INHBA")
DotPlot(pbmc_myeloid, features = Macro_SPP1_features, cols = cont_2) + RotatedAxis()
# 3 4 10 11 16 

# Doublets-B
Doublets_B_features = c("CD79A", "CD79B","MZB1", "MS4A1", "JCHAIN", "DERL3", "FKBP11", "TNFRSF17")
DotPlot(pbmc_myeloid, features = Doublets_B_features, cols = cont_2) + RotatedAxis()
# 17

# Doublets-Epi
Doublets_Epi_features = c("EPCAM", "SMIM22", "PRSS3", "PERP", "KLF5", "C19orf33", "CEACAM5")
DotPlot(pbmc_myeloid, features = Doublets_Epi_features, cols = cont_2) + RotatedAxis()
# 20

# Doublets-T
Doublets_T_features = c("CD3D", "CD3E", "CD3G", "GNLY", "IL32", "CCL5", "CD2", "GZMA", "KLRB1", "TBRC2")
DotPlot(pbmc_myeloid, features = Doublets_T_features, cols = cont_2) + RotatedAxis()
# 7


##############################
# Checking Xie at al markers part 2
##############################

# analysis of feature distribution across different macrophages
FeaturePlot(pbmc_myeloid, features = Macro_C1QC_features,
            cols = cont_2 ,reduction = "umap", label = T)
FeaturePlot(pbmc_myeloid, features = Macro_SPP1_features,
            cols = cont_2 ,reduction = "umap", label = T)
FeaturePlot(pbmc_myeloid, features = Macro_FCN1_features,
            cols = cont_2 ,reduction = "umap", label = T)

#build a cluster tree
pbmc_myeloid <- BuildClusterTree(
  pbmc_myeloid,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
PlotClusterTree(pbmc_myeloid)


##############################
# Checking Xie at al markers part 3
##############################

# CDC identification
FeaturePlot(pbmc_myeloid, features = c("LILRA4","CLEC9A","CD1C","LAMP3"),
            cols = cont_2 ,reduction = "umap")
# Macrophage identification
FeaturePlot(pbmc_myeloid, features = c("SPP1","C1QC","FCN1","MKI67"),
            cols = cont_2 ,reduction = "umap")

# Final comparison to Xie et al
features = rev(c(
  "CD68", "CD163", "CD14", "CD3D", "EPCAM", "CD79A", "MKI67",
  "SPP1", "CSTB", "SDC2", "FCN1", "S100A9", "S100A8",
  "C1QC", "C1QB", "C1QA", "LTB", "LAMP3", "CD1C", "CLEC9A", "LILRA4", "CES1"
))

# Dot plot group comparison with Xie et al paper
Idents(pbmc_myeloid) <- "seurat_clusters"
Idents(pbmc_myeloid) <- factor(
  Idents(pbmc_myeloid),
  levels = c(
    "1","2","5",             # Macro-FCN1
    "11", "3","4","10","16", # Macro-SPP1
    "0","6","12","14","18",  # Macro-C1QC
    "8",                     # Macro-MKI67
    "19","9","13","15",      # cDC1 → cDC2 → cDC3 → pDC
    "7","17","20"
  )
)

# analysis using clustering information
b2 <- DimPlot(pbmc_myeloid, reduction = 'umap',label=F, 
              group.by = 'seurat_clusters', cols = distinctColorPalette(21)) +
  ggtitle(NULL) 
b2

# dotplot all main features
c2 <- DotPlot(pbmc_myeloid, features = unique(features), cols = cont_2) + 
  RotatedAxis() + coord_flip() + 
  labs(y = "Cluster Number")
c2


##############################
# Label and visualise
##############################

# monocytes and FCN1 have a pretty overlapping definition
# label the cells
cluster5celltype <-c("0"="Macro-C1QC", #
                     "1"="Macro-FCN1", #
                     "2"="Macro-FCN1", #
                     "3"="Macro-SPP1", #
                     "4"="Macro-SPP1", #
                     "5"="Macro-FCN1", #
                     "6"="Macro-C1QC", #
                     "7"="Doublets-T", #
                     "8"="Macro-MKI67", #
                     "9"="cDC2-CD1C", #
                     "10"="Macro-SPP1", #
                     "11"="Macro-SPP1", #
                     "12"="Macro-C1QC", #
                     "13"="cDC3-LAMP3", #
                     "14"="Macro-C1QC", #
                     "15"="pDC-LILRA4", #
                     "16"="Macro-SPP1", #
                     "17"="Doublets-B", #
                     "18"="Macro-C1QC", #
                     "19"="cDC1-CLEC9A", #
                     "20"="Doublets-Epi" #
)

# add cluster labels to pbmc
pbmc_myeloid[['myeloid_type']] = unname(cluster5celltype[pbmc_myeloid@meta.data$seurat_clusters])

# switch between setting Idents and myeloid
Idents(pbmc_myeloid) <- "myeloid_type"

# save labelled myeloid dimplot
d2 <-DimPlot(pbmc_myeloid, reduction = "umap", 
        group.by = 'myeloid_type',label = F, cols = distinctColorPalette(11)) +
        ggtitle(NULL)
d2


##############################
# Myeloid/Meta data Check
##############################

# patient distribution
DimPlot(object = pbmc_myeloid, group.by="orig.ident", reduction='umap')
# interestingly, the entire MMRP MLH hypermethylation group of SPP1 are from one patient
pbmc_myeloid_C165 <-  subset(pbmc_myeloid, subset =  startsWith(as.character(orig.ident), "C165"))
DimPlot(object = pbmc_myeloid_C165, group.by="orig.ident", reduction='umap')

# metastasis status
# one group is called:
# pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum.
# change to pM1c
pbmc_myeloid@meta.data$MetastasisStatus[pbmc_myeloid@meta.data$MetastasisStatus == 
                                          "pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum."] <- "pM1c"
DimPlot(object = pbmc_myeloid, group.by="MetastasisStatus",reduction='umap')

# tissue site
DimPlot(object = pbmc_myeloid, group.by="TissueSiteSimple",reduction='umap')

# sex
DimPlot(object = pbmc_myeloid, group.by="Sex",reduction='umap')

# N or T
DimPlot(object = pbmc_myeloid, group.by="SPECIMEN_TYPE",reduction='umap')

# mmrstatus
DimPlot(object = pbmc_myeloid, group.by="MMRStatus",reduction='umap', cols = disc_10)


###########################################################
# Part 2
###########################################################
##############################
# Analysis of CES1 in clusters
##############################

# separate into a macrophage subset
pbmc_macrophage<-subset(pbmc_myeloid,myeloid_type %in% c("Macro-C1QC","Macro-FCN1",
                                              "Macro-MKI67","Macro-SPP1"))

# CES1 expression dot plot
a3 <- DotPlot(pbmc_macrophage, features = 'CES1', cols = cont_2) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
a3


##############################
# Further CES1 expression analysis and Statistics
##############################

# utilise box plots, and violin plots to investigate CES1: Percent Expressed and Average Expression
# also analyse MMRp MMRd distribution in groups and test significance.

##############################
# Initial Exploratory Figures
##############################


# Firstly plot the total number of macrophages expressing CES1
ggplot(
  pbmc_macrophage@meta.data,
  aes(x = CES1_isExpressed, fill = CES1_isExpressed)
) +
  geom_bar() +
  scale_fill_manual(values = cont_2) +
  scale_x_discrete(labels = c("Not expressing", "Expressing")) +
  labs(x = "CES1 expression", y = "Number of cells") +
  theme_minimal() +
  theme(legend.position = "none")

# plot the Percentage of Normal and Tumor cells in each macrophage type that are expressing CES1
b3 <- ggplot(
  pbmc_macrophage@meta.data %>% filter(CES1_isExpressed),
  aes(x = myeloid_type, fill = SPECIMEN_TYPE)
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
    axis.title.x=element_blank(), 
    legend.position = "none")

b3

# violin plot of CES1 expression
VlnPlot(pbmc_macrophage, features = c("CES1"), ncol = 3)
# violins arent visible 

# Violin plot no NAs
# to view the violin plots remove all 0 values
pbmc_macrophage_CES1exp <- pbmc_macrophage[,pbmc_macrophage$CES1_isExpressed == TRUE]
# violin plot
VlnPlot(pbmc_macrophage_CES1exp, features = c("CES1"), ncol = 3, cols = disc_10)


##############################
# CES1 average expression
##############################

# make a box plot of average CES1 expression for each patient and group by myeloid_type
# first make a summary dataframe
df_avg_CES1_exp <-  pbmc_macrophage@meta.data %>%
  # group by
  group_by(orig.ident, myeloid_type, SPECIMEN_TYPE) %>%
  # the summarise function
  summarise(avg_CES1_exp = mean(CES1_expression))

# add comparisons values for later stats test
m_comparisons <- list( c("Macro-C1QC", "Macro-SPP1"))

# make a boxplot depicting the average CES1 expression per patient
d3 <- ggplot(df_avg_CES1_exp, 
             aes(x = myeloid_type, y = avg_CES1_exp, fill = myeloid_type)) +
  # add notches to the boxplot
  geom_boxplot( notch=TRUE, notchwidth = 0.8) +
  scale_fill_manual(values = distinctColorPalette(4)) +
  labs(x = "Cell type", y = "Average Expression") +
  # add statistical test between C1QC and SPP1
  stat_compare_means(comparisons = m_comparisons, 
                     aes(label = after_stat(p.signif)),
                     method = "t.test") +
  theme_minimal() +
  # clean up for figure
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank())
d3

# add comparisons values for later stats test
m2_comparisons <- list( c("T", "N"))

# make a boxplot depicting the N vs T
c3 <- ggplot(df_avg_CES1_exp %>% filter(myeloid_type %in% c("Macro-SPP1", "Macro-C1QC")), 
             aes(x = SPECIMEN_TYPE, y = avg_CES1_exp, fill = SPECIMEN_TYPE)) +
  scale_fill_discrete("SPECIMEN_TYPE") +
  geom_boxplot(position="dodge", notch=TRUE, notchwidth = 0.8) +
  scale_fill_manual(values = disc_10) +
  labs(x = "Cell type", y = "Average Expression") +
  facet_wrap(~myeloid_type) +
  # do a stats test between the MMRstatuses
  stat_compare_means(comparisons = m2_comparisons,
                     aes(label = after_stat(p.signif)),
                     method = "t.test", paired = FALSE) +
  theme_minimal()  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
c3

##############################
# CES1 percentage expressed
##############################

# make a box plot of Percentage Expressed for each patient and group by myeloid_type
# first make a summary dataframe
exp_pct_S_C <-  pbmc_macrophage@meta.data %>%
  # group by
  group_by(orig.ident, myeloid_type) %>%
  ##the summarise function
  summarise(percent_CES1_exp = mean(CES1_isExpressed) * 100)

# make a boxplot depicting the percentage CES1 expressed per patient
e3 <- ggplot(exp_pct_S_C, aes(x = myeloid_type, y = percent_CES1_exp, fill = myeloid_type)) +
  # add notches
  geom_boxplot(notch=TRUE, notchwidth = 0.8) +
  scale_fill_manual(values = distinctColorPalette(4)) +
  labs(x = "Cell type", y = "Percent Expressed") +
  # add stats test 
  stat_compare_means(comparisons = m_comparisons, 
                     aes(label = after_stat(p.signif)),
                     method = "t.test") +
  theme_minimal() +
  # make clean for figure
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())
e3


# Figure 2
ggdraw() +
  draw_plot(a2, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(b2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(c2, x = 0, y = 0, width = .5, height = .4) +
  draw_plot(d2, x = .5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15, 
                  x = c(0, .5, 0, .5), 
                  y = c(1, 1, .5, .5))
ggsave("../figures/Figure_2_myeloid_general.jpg", width = 40, height = 25, units = c("cm"), dpi = 300)

# Figure 3
ggdraw() +
  draw_plot(a3, x = 0, y = .5, width = .33, height = .5) +
  draw_plot(b3, x = .33, y = .5, width = .23, height = .5) +
  draw_plot(c3, x = .56, y = .5, width = .43, height = .5) +
  draw_plot(d3, x = .1, y = 0, width = .33, height = .5) +
  draw_plot(e3, x = .44, y = 0, width = .48, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15, 
                  x = c(0, .33, .56, .1, .44), 
                  y = c(1, 1, 1, .5, .5))
ggsave("../figures/Figure_3_myeloid_ces1.jpg", width = 30, height = 20, units = c("cm"), dpi = 300)


##############################
# Pre cell chat cell type labelling
##############################

# CellChat requires a single group with the cell information
# my information is split between a few variables so i will change it to one column

# change the cl295... column for myeloid_type in pbmc only for myeloid cells
pbmc@meta.data$cell_type_detailed <- pbmc@meta.data$cl295v11SubFull
head(pbmc@meta.data)

cells <- rownames(pbmc_myeloid@meta.data)
pbmc@meta.data[cells, "cell_type_detailed"] <- paste0("(", pbmc_myeloid@meta.data[cells, "myeloid_type"], ")")

# make a new label column for cellchat that 
pbmc@meta.data$cell_type_detailed <-
  paste0(pbmc@meta.data$clTopLevel, ": ", 
         str_extract(pbmc@meta.data$cell_type_detailed, "(?<=\\().*(?=\\))"))


##############################
# Final Saves
##############################

# saveRDS(pbmc_myeloid, file = "../data/pbmc_myeloid_final.rds")
# pbmc_myeloid <- readRDS("../data/pbmc_myeloid_final.rds")

# saveRDS(pbmc_myeloid.markers, file = "../data/pbmc_myeloid.markers_PC30.rds")
# pbmc_myeloid.markers <- readRDS("../data/pbmc_myeloid.markers_PC30.rds")

# saveRDS(pbmc_macrophage, file = "../data/pbmc_macrophage_final.rds")
# pbmc_macrophage <- readRDS("../data/pbmc_macrophage_final.rds")

# saveRDS(pbmc, file = "../data/pbmc_final.rds")
# pbmc <- readRDS("../data/pbmc_final.rds")

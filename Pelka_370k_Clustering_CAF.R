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

# Establish some colours
cont_2 <- brewer.pal(9, "YlOrRd")[c(1, 9)]
grey_red <- c("lightgrey", "#b81f25")
disc_10 <- brewer.pal(10, "Set3")


##############################
# Subset analysis
##############################

# as the CAF class has also been identified to have high expression of CES1, 
# they will be further analysed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################
# Stromal subset
##############################

# initially separate my cells identified as stromal
pbmc_stromal <- subset(pbmc,celltype %in% c("Stromal"))
# identify how many cells in subset
nrow(pbmc_stromal@meta.data)

# then use the clustering annotations to only keep cells identified as stromal by pelka
pbmc_stromal <- subset(pbmc_stromal, subset =  startsWith(cl295v11SubFull, "cS"))
# identify how many are lost
nrow(pbmc_stromal@meta.data)
# only 41 cells lost

# list the unique cell type names
unique(pbmc_stromal@meta.data$cl295v11SubFull)

# change names to only contain the cell types
pbmc_stromal@meta.data$cl295v11SubFull <- 
  gsub(".*\\((.*)\\).*", "\\1", 
       pbmc_stromal@meta.data$cl295v11SubFull)

##############################
# Fibroblast subset
##############################

# this note from pelka needs to be considered
### Note: cS30 and cS31 are overwhelmingly from two tumors which grew below non-neoplastic tissue and 
# may not be purely tumor-derived

# keep only fibroblasts
pbmc_fibroblast <- subset(pbmc_stromal, 
                           cl295v11SubFull %in% c("Fibro stem cell niche",
                                                  "Fibro BMP-producing",
                                                  "Fibro CCL8+",
                                                  "Myofibro",
                                                  "CXCL14+ CAF",
                                                  "GREM1+ CAF",
                                                  "MMP3+ CAF"))
# check how many cells
nrow(pbmc_fibroblast@meta.data)
# over half


##############################
# Data Preprocessing and Clustering
##############################

## normalization
pbmc_fibroblast <- NormalizeData(pbmc_fibroblast, normalization.method = "LogNormalize", scale.factor = 10000) 
pbmc_fibroblast <- FindVariableFeatures(pbmc_fibroblast, selection.method = 'vst', nfeatures = 2000)

# scaling
pbmc_fibroblast <- ScaleData(pbmc_fibroblast, vars.to.regress = "percent.mt")
pbmc_fibroblast <- RunPCA(pbmc_fibroblast, features = VariableFeatures(object = pbmc_fibroblast)) 

# dim use identification
ElbowPlot(pbmc_fibroblast, reduction = "pca", ndims = 50)

# dimensional reduction
dim.use<-1:20
pbmc_fibroblast <- FindNeighbors(pbmc_fibroblast, dims = dim.use)
pbmc_fibroblast <- FindClusters(pbmc_fibroblast, resolution = 0.8 )
pbmc_fibroblast <- RunUMAP(pbmc_fibroblast, dims = dim.use)

# plot clusters
b1 <- DimPlot(pbmc_fibroblast, reduction = 'umap',label=F, group.by = 'cl295v11SubFull') +
  ggtitle(NULL)
b1
Idents(pbmc_fibroblast) <- "cl295v11SubFull"


##############################
# Finding differentially expressed genes
##############################

# find the markers for the fibroblast clusters
pbmc_fibroblast.markers <- FindAllMarkers(pbmc_fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_fibroblast.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# saveRDS(pbmc_fibroblast.markers, file = "../data/pbmc_fibroblast.markers.rds")
# pbmc_fibroblast.markers <- readRDS("../data/pbmc_fibroblast.markers.rds")


##############################
# fibroblast/Meta data Check
##############################

# patient distribution
DimPlot(object = pbmc_fibroblast, group.by="orig.ident", reduction='umap')

# metastasis status
# one group is called:
# pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum.
# change to pM1c
pbmc_fibroblast@meta.data$MetastasisStatus[pbmc_fibroblast@meta.data$MetastasisStatus == 
                                          "pM1c (Metastases the peritoneal surface, alone or with other site or organ metastases): Sites involved: Liver and peritoneum."] <- "pM1c"
DimPlot(object = pbmc_fibroblast, group.by="MetastasisStatus",reduction='umap')

# tissue site
DimPlot(object = pbmc_fibroblast, group.by="TissueSiteSimple",reduction='umap')

# sex
DimPlot(object = pbmc_fibroblast, group.by="Sex",reduction='umap')

# N or T
DimPlot(object = pbmc_fibroblast, group.by="SPECIMEN_TYPE",reduction='umap')

# mmrstatus
DimPlot(object = pbmc_fibroblast, group.by="MMRStatus",reduction='umap')


##############################
# Analysis of CES1 in clusters
##############################

# dot plot of the fibroblast cells
DotPlot(pbmc_fibroblast, features = 'CES1', cols = cont_2, group.by = "cl295v11SubFull") + RotatedAxis()
# fibro stem cell niche is high but these are alll normal cells so make CAF only calss

# dotplot of just the CAFs
pbmc_CAF <- subset(pbmc_fibroblast, cl295v11SubFull %in% c("CXCL14+ CAF", "GREM1+ CAF", "MMP3+ CAF", "Myofibro"))
d1 <- DotPlot(pbmc_CAF, features = 'CES1', cols = cont_2, group.by = "cl295v11SubFull") + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
d1


# bar chart showing the distribution of normal and tumour cells between fibroblasts 
e1 <- ggplot(
  pbmc_fibroblast@meta.data,
  aes(x = cl295v11SubFull, fill = SPECIMEN_TYPE)) +
  geom_bar(position="dodge") +
  scale_fill_manual(values = disc_10) +
  labs(y = "Number of cells") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x=element_blank())
e1


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
f1 <- ggplot(df_avg_CES1_exp_f, aes(x = cl295v11SubFull, y = avg_CES1_exp, fill = cl295v11SubFull)) +
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
f1

# make a boxplot depicting the avg ces1 exp for just the myofib class
h1 <- ggplot(df_avg_CES1_exp_f %>% filter(cl295v11SubFull %in% c("Myofibro")),
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
h1


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
g1 <- ggplot(exp_pct_S_C_f, aes(x = cl295v11SubFull, y = percent_CES1_exp, fill = cl295v11SubFull)) +
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
g1

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

# saveRDS(pbmc_fibroblast.markers, file = "../data/pbmc_fibroblast.markers.rds")
# pbmc_fibroblast.markers <- readRDS("../data/pbmc_fibroblast.markers.rds")

# saveRDS(pbmc_fibroblast, file = "../data/pbmc_fibroblast_final.rds")
# pbmc_fibroblast <- readRDS("../data/pbmc_fibroblast_final.rds")


##############################
# Build final figures (Figure 1)
##############################


# Figure 1
ggdraw() +
  draw_plot(b1, x = 0, y = .33, width = .5, height = .66) +
  draw_plot(a1, x = .06, y = .52, width = .1, height = .2) +
  draw_plot(c1, x = .5, y = .33, width = .17, height = .66) +
  draw_plot(d1, x = .67, y = .33, width = .33, height = .66) +
  draw_plot(e1, x = 0, y = 0, width = .25, height = .33) +
  draw_plot(f1, x = .25, y = 0, width = .18, height = .33) +
  draw_plot(g1, x = .43, y = 0, width = .32, height = .33) +
  draw_plot(h1, x = .75, y = 0, width = .25, height = .33) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H"), size = 15, 
                  x = c(0, .062, .5, .65, 0, .25, .43, .75), 
                  y = c(1, .71, 1, 1, .37, .37, .37, .37))

ggsave("../figures/Figure_1_CAF.jpg", width = 45, height = 20, units = c("cm"), dpi = 300)

# Figure 2
ggdraw() +
  draw_plot(b2, x = 0, y = .33, width = .5, height = .66) +
  draw_plot(a1, x = .23, y = .50, width = .1, height = .2) +
  draw_plot(c2, x = .5, y = .33, width = .17, height = .66) +
  draw_plot(d2, x = .67, y = .33, width = .33, height = .66) +
  draw_plot(e2, x = 0, y = 0, width = .25, height = .33) +
  draw_plot(f2, x = .25, y = 0, width = .18, height = .33) +
  draw_plot(g2, x = .43, y = 0, width = .32, height = .33) +
  draw_plot(h2, x = .75, y = 0, width = .25, height = .33) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H"), size = 15, 
                  x = c(0, .232, .5, .65, 0, .25, .43, .75), 
                  y = c(1, .69, 1, 1, .37, .37, .37, .37))

ggsave("../figures/Figure_2_Macrophage.jpg", width = 45, height = 20, units = c("cm"), dpi = 300)


    
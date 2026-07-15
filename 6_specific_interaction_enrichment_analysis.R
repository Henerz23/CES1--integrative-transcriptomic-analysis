###########################################################
# SPECIFIC INTERACTION ANALYSIS AND ENRICHMENT
###########################################################

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(patchwork)

pbmc <- readRDS("../data/pbmc_final.rds")
results <- readRDS("../data/nichenet_results_300.rds")
lr_networks <- readRDS("../data/lr_networks.rds")
groups <- readRDS("../data/nichenet_gene_groups.rds")

# focus on SPP1- CD8
# individual heat maps for specific interaction analysis
# remove the first column

# I REMOVED APP!!!!!!


# include only gene which are exclusive to SPP1
spp1_up_only_conditions <- c(groups$`lr_spp1_TCD8_up///lr_spp1_TCD4_up`, 
                          groups$lr_spp1_TCD8_up, 
                          groups$lr_spp1_TCD4_up)

heatmap_spp1_cd8_up <- results$lr_spp1_TCD8_up$gene_expression[, -c(1, 4, 5, 6)]
heatmap_spp1_cd8_up <- heatmap_spp1_cd8_up[rownames(heatmap_spp1_cd8_up) %in% spp1_up_only_conditions, ]


nrow(heatmap_spp1_cd8_up)
keep <- apply(heatmap_spp1_cd8_up, 1, function(x) max(x) >= 0.04)
heatmap_spp1_cd8_up <- heatmap_spp1_cd8_up[keep, , drop = FALSE]
nrow(heatmap_spp1_cd8_up)


a6 <- heatmap_spp1_cd8_up %>% make_heatmap_ggplot(
    y_name = paste("Target genes in "),
    x_name = paste("prioritized ligands "),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) 
a6



# NISH
spp1_down_only_conditions <- c(groups$lr_spp1_TCD4_down, 
                             groups$lr_spp1_TCD8_down)
heatmap_spp1_cd8_down <- results$lr_spp1_TCD8_down$gene_expression# [, -c(1, 2, 3, 6, 8)]
heatmap_spp1_cd8_down <- heatmap_spp1_cd8_down[rownames(heatmap_spp1_cd8_down) %in% spp1_down_only_conditions, ]
"""
nrow(heatmap_spp1_cd8_down)
keep <- apply(heatmap_spp1_cd8_down, 1, function(x) max(x) - min(x) >= 0.04)
heatmap_spp1_cd8_down <- heatmap_spp1_cd8_down[keep, , drop = FALSE]
nrow(heatmap_spp1_cd8_down)
"""
heatmap_spp1_cd8_down %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 


# C1QC- CD8
# NISH



# C1QC- CD4
# include only gene which are exclusive to c1qc
c1qc_up_only_conditions <- c(groups$lr_c1qc_TCD4_up)

heatmap_c1qc_cd4_up <- results$lr_c1qc_TCD4_up$gene_expression[, -c(1, 4)]
heatmap_c1qc_cd4_up <- heatmap_c1qc_cd4_up[rownames(heatmap_c1qc_cd4_up) %in% c1qc_up_only_conditions, ]


nrow(heatmap_c1qc_cd4_up)
keep <- apply(heatmap_c1qc_cd4_up, 1, function(x) max(x) >= 0.015)
heatmap_c1qc_cd4_up <- heatmap_c1qc_cd4_up[keep, , drop = FALSE]
nrow(heatmap_c1qc_cd4_up)

c6 <- heatmap_c1qc_cd4_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
c6

# cd4 up
heatmap_spp1_cd4_up <- results$lr_spp1_TCD4_up$gene_expression[, -c(1, 2, 3, 6, 8)]
heatmap_spp1_cd4_up <- heatmap_spp1_cd4_up[rownames(heatmap_spp1_cd4_up) %in% spp1_up_only_conditions, ]


nrow(heatmap_spp1_cd4_up)
keep <- apply(heatmap_spp1_cd4_up, 1, function(x) max(x) >= 0.04)
heatmap_spp1_cd4_up <- heatmap_spp1_cd4_up[keep, , drop = FALSE]
nrow(heatmap_spp1_cd4_up)

b6 <- heatmap_spp1_cd4_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
b6


# SPP1 CD4 down
# NISH
heatmap_spp1_cd4_down <- results$lr_spp1_TCD4_down$gene_expression# [, -c(1, 2, 3, 6, 8)]
heatmap_spp1_cd4_down <- heatmap_spp1_cd4_down[rownames(heatmap_spp1_cd4_down) %in% spp1_down_only_conditions, ]
"""
nrow(heatmap_spp1_cd4_down)
keep <- apply(heatmap_spp1_cd4_down, 1, function(x) max(x) - min(x) >= 0.04)
heatmap_spp1_cd4_down <- heatmap_spp1_cd4_down[keep, , drop = FALSE]
nrow(heatmap_spp1_cd4_down)
"""
heatmap_spp1_cd4_down %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
















# receptors
# visualise receptor interactions for specific interactions
receptors_spp1_cd4_up <- results$lr_spp1_TCD4_up$receptor_expression[, -c(1, 2, 4, 6)]

# d6 <- t(receptors_spp1_cd8_up) %>%
  make_heatmap_ggplot(
    y_name       = "Ligands",
    x_name       = "Receptors expressed in receiver (CD4+ T cells)",
    color        = "mediumvioletred",
    legend_title = "Prior interaction\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) +
  coord_flip()
# d6




# panel D will be some kind of key referencing where the genes 
# are uniquely up/downreuglated in that condition

# Panel E some kind of table showing what the effects are of each ligand
# Panel F some kind of diagram showing pathways
d6 <- 0
e6 <- 0
f6 <- 0
ggdraw() +
  draw_plot(a6, x = 0, y = .5, width = .33, height = .5) +
  draw_plot(b6, x = .33, y = .5, width = .33, height = .5) +
  draw_plot(c6, x = .66, y = .5, width = .33, height = .5) +
  draw_plot(d6, x = 0, y = 0, width = .33, height = .5) +
  draw_plot(e6, x = .33, y = 0, width = .33, height = .5) +
  draw_plot(f6, x = .66, y = 0, width = .33, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 15, 
                  x = c(0, .33, .66, 0, .33, .66), 
                  y = c(1, 1, 1, .5, .5, .5))
ggsave("../figures/Figure_6_nichenet_enrichment.jpg", width = 30, height = 20, units = c("cm"), dpi = 300)



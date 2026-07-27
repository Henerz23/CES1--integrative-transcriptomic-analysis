###########################################################
# SPECIFIC INTERACTION ANALYSIS AND ENRICHMENT
###########################################################

# Part 1: Analysis of specific ligands and upregulated genes
# Part 2: GSEA enrichment analysis
# Warning: due to iterative enrichment investigation style, this script has been hardcoded 
# Thresholds were handtuned for figure clarity


##############################
# Load Libraries
##############################

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(patchwork)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load in dependencies
pbmc <- readRDS("../data/pbmc_final.rds")
results <- readRDS("../data/nichenet_results_300.rds")
lr_networks <- readRDS("../data/lr_networks.rds")
groups <- readRDS("../data/nichenet_gene_groups.rds")
weighted_networks  <- readRDS("../data/weighted_networks_human.rds")


###########################################################
# Part 1
###########################################################
##############################
# Exclusive gene upregulation analysis
##############################

# include only gene which are exclusive to SPP1
spp1_up_only_conditions <- c(groups$`lr_spp1_TCD8_up///lr_spp1_TCD4_up`, 
                          groups$lr_spp1_TCD8_up, 
                          groups$lr_spp1_TCD4_up)

# SPP1 CD8
# keep only matrix columns of interest
heatmap_spp1_cd8_up <- results$lr_spp1_TCD8_up$gene_expression[, colnames(results$lr_spp1_TCD8_up$gene_expression) %in% c("RETN", "APP")]
heatmap_spp1_cd8_up <- heatmap_spp1_cd8_up[rownames(heatmap_spp1_cd8_up) %in% spp1_up_only_conditions, ]

# only keep genes where the regulatory potential is high
nrow(heatmap_spp1_cd8_up)
keep <- apply(heatmap_spp1_cd8_up, 1, function(x) max(x) >= 0.04)
heatmap_spp1_cd8_up <- heatmap_spp1_cd8_up[keep, , drop = FALSE]
nrow(heatmap_spp1_cd8_up)

# display heatmap
a6 <- heatmap_spp1_cd8_up %>% make_heatmap_ggplot(
    y_name = paste("Target genes in "),
    x_name = paste("prioritized ligands "),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) 
a6

# SPP1 CD4
# keep only matrix columns of interest
heatmap_spp1_cd4_up <- results$lr_spp1_TCD4_up$gene_expression[, colnames(results$lr_spp1_TCD4_up$gene_expression) %in% c("RETN", "APP", "CCL20")]
heatmap_spp1_cd4_up <- heatmap_spp1_cd4_up[rownames(heatmap_spp1_cd4_up) %in% spp1_up_only_conditions, ]

# only keep genes where the regulatory potential is high
nrow(heatmap_spp1_cd4_up)
keep <- apply(heatmap_spp1_cd4_up, 1, function(x) max(x) >= 0.04)
heatmap_spp1_cd4_up <- heatmap_spp1_cd4_up[keep, , drop = FALSE]
nrow(heatmap_spp1_cd4_up)

# display heatmap
b6 <- heatmap_spp1_cd4_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
b6

# C1QC- CD8
# No genes present

# C1QC- CD4
# include only gene which are exclusive to c1qc
c1qc_up_only_conditions <- c(groups$lr_c1qc_TCD4_up)

# keep only matrix columns of interest
heatmap_c1qc_cd4_up <- results$lr_c1qc_TCD4_up$gene_expression[, colnames(results$lr_c1qc_TCD4_up$gene_expression) %in% c("HLA-DMB", "HLA-DQA2")]
heatmap_c1qc_cd4_up <- heatmap_c1qc_cd4_up[rownames(heatmap_c1qc_cd4_up) %in% c1qc_up_only_conditions, ]

# only keep genes where the regulatory potential is high
nrow(heatmap_c1qc_cd4_up)
keep <- apply(heatmap_c1qc_cd4_up, 1, function(x) max(x) >= 0.015)
heatmap_c1qc_cd4_up <- heatmap_c1qc_cd4_up[keep, , drop = FALSE]
nrow(heatmap_c1qc_cd4_up)

# display heatmap
c6 <- heatmap_c1qc_cd4_up %>% make_heatmap_ggplot(
  y_name = paste("Target genes in "),
  x_name = paste("prioritized ligands "),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic")) 
c6


##############################
# LR-gene pathway analysis
##############################

# Pathway analysis for individual ligand - receptor - gene sequences

# Infer Signaling Paths from Top Ligands to Target Genes
# Load the ligand-TF matrix — a prebuilt matrix (ligands x TFs) encoding
ligand_tf_matrix <- readRDS("../data/ligand_tf_matrix_human.rds")

# Choose a top ligand to trace the signaling path
top_ligand <- "RETN"

# Choose the top predicted target gene for that ligand (from Step 9)
top_target <- "AHSA1"

# Infer the signaling path from top_ligand to top_target
sig_path <- get_ligand_signaling_path(
  ligand_tf_matrix  = ligand_tf_matrix,
  ligands_all       = top_ligand,
  targets_all       = top_target,
  top_n_regulators  = 4,
  weighted_networks = weighted_networks,
  minmax_scaling    = TRUE  # normalise edge weights for cleaner visualisation
)

# Combine both edge tables and build an igraph for visualisation
sig_path_combined <- bind_rows(sig_path$sig, sig_path$gr)

# Build igraph and visualise with ggraph
sig_graph <- graph_from_data_frame(sig_path_combined, directed = TRUE)

#plot
ggraph(sig_graph, layout = "kk") +
  geom_edge_link(
    aes(edge_width = weight, edge_alpha = weight),
    arrow    = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap  = circle(3, "mm"),
    color    = "grey50"
  ) +
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = paste("Signaling path:", top_ligand, "->", top_target))

# limitation- Nichenet is a bit of a blackbox in that inferring exact pathways is hard


###########################################################
# Part 2
###########################################################
##############################
# Load Libraries
##############################

library(clusterProfiler)  
library(org.Hs.eg.db)
library(tidyverse)
library(DOSE)
library(ReactomePA)
library(enrichplot)

##############################
# GSEA enrichment analysis
##############################

# read in the heatmap file
mega_heatmap_up <- readRDS("../data/mega_heatmap_up.rds")

# run a loop for only SPP1-CD8/CD4 only and C1QC-CD4 only
x = c("_1", "_2", "_5")
plots <- list()

# loop over all DE genes for each group exclusive to SPP1 or C1QC
# identifies the main biological processes that the genes are driving 
for (numb in x) {
  genes <- mega_heatmap_up %>%
    as.data.frame() %>%
    mutate(
      gene = rownames(.),
      max_rp = do.call(pmax, c(dplyr::select(., ends_with(numb)), na.rm = TRUE))
    ) %>%
    filter(max_rp > 0) %>%
    arrange(desc(max_rp)) %>%
    pull(gene)
  
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # specify background expressed gene set (still gene symbols at this point)
  background_expressed_genes <- switch(
    numb,
    "_1" = results$lr_spp1_TCD8_up$background_expressed_genes,
    "_2" = results$lr_spp1_TCD4_up$background_expressed_genes,
    "_5" = results$lr_c1qc_TCD4_up$background_expressed_genes
  )
  
  # convert the universe to ENTREZ IDs too, since keyType = "ENTREZID" below
  universe_ids <- bitr(background_expressed_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  ego <- enrichGO(
    gene = ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    universe = universe_ids
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    message("No enriched GO terms for ", numb, " — skipping")
    next
  }
  
  s_ego<-clusterProfiler::simplify(ego)
  
  # plot results and add to list
  plots[[numb]] <-  s_ego %>% filter(p.adjust < 0.03) %>%
    ggplot(showCategory = 20,
           aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range=c(1, 7)) +
    theme_minimal() + 
    xlab("Gene Ratio") +
    ylab(NULL) + 
    ggtitle("GO Enrichment of up-regulated genes")
}


# add list to figure variables
d6 <- plots[["_1"]]
e6 <- plots[["_2"]]
f6 <- plots[["_5"]]


##############################
# Make figures
##############################

# make figure 6
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
# save figure
ggsave("../figures/Figure_6_nichenet_enrichment.jpg", width = 40, height = 20, units = c("cm"), dpi = 300)


###########################################################
# SPECIFIC INTERACTION ANALYSIS AND ENRICHMENT
###########################################################

# Part 1: Complex heatmap- finsing genes in overlapping conditions
# Part 2: Analysis of specific ligands and upregulated genes
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
library(cowplot)

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
# identify genes expressed through specific types of interaction
##############################


# load appropriate library
library(ComplexUpset)

geneset_lists <- list()

# loop version
# loop over the different interaction conditions
for (nm in names(results)) {
  # include gene expression
  mat <- results[[nm]]$gene_expression
  # filter to include only genes in mega heatmap
  # skip if NULL or empty
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    next
  }
  # ensure matrix
  mat <- as.matrix(mat)
  # keep genes with any expression signal
  expressed_genes <- rownames(mat)[rowSums(!is.na(mat) & mat > 0.05) > 0]
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
b5 <- ComplexUpset::upset(
  geneset_df,
  intersect = colnames(geneset_df),
  sort_intersections_by = "degree",
  sort_sets = "descending",
  keep_empty_groups = TRUE,
  base_annotations = list("Intersection size" = intersection_size()))
b5

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

# save groups
saveRDS(groups, file = "../data/nichenet_gene_groups.rds")


###########################################################
# Part 1
###########################################################
##############################
# Exclusive gene upregulation analysis
##############################


# SPP1 CD8
# keep only matrix columns of interest
heatmap_spp1_cd8_up <- results$lr_spp1_TCD8_up$gene_expression
# filter to only include genes in spp1 cd8 
heatmap_spp1_cd8_up <- heatmap_spp1_cd8_up[rownames(heatmap_spp1_cd8_up) %in% c(groups$lr_spp1_TCD8_up, groups$`lr_spp1_TCD8_up///lr_spp1_TCD4_up`), ]

# display heatmap
c5 <- heatmap_spp1_cd8_up %>% make_heatmap_ggplot(
    y_name = paste("Target genes"),
    x_name = paste("Prioritised Ligands"),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_text(angle = 30)) +
  ggtitle("Unique Genes from \nlr_spp1_TCD8_up")
c5


# SPP1 CD4
# keep only matrix columns of interest
heatmap_spp1_cd4_up <- results$lr_spp1_TCD4_up$gene_expression
heatmap_spp1_cd4_up <- heatmap_spp1_cd4_up[rownames(heatmap_spp1_cd4_up) %in% c(groups$lr_spp1_TCD4_up, groups$`lr_spp1_TCD8_up///lr_spp1_TCD4_up`), ]

# display heatmap
d5 <- heatmap_spp1_cd4_up %>% make_heatmap_ggplot(
  y_name = paste("Target Genes"),
  x_name = paste("Prioritised Ligands"),
  color = "purple",
  legend_title = "Regulatory\npotential"
) +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_text(angle = 30)) +
  ggtitle("Unique Genes from \nlr_spp1_TCD4_up")
d5
?theme()
# C1QC- CD8
# No genes present

# C1QC- CD4
# no genes present

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

x <- c("RETN_1", "RETN_2", "APP_1", "APP_2")
group_labels <- c("RETN_1" = "RETN-CD8", "RETN_2" = "RETN-CD4", "APP_1" = "APP-CD8", "APP_2" = "APP-CD4")

# build per-group gene lists (ENTREZ) + collect backgrounds
gene_list <- list()
background_list <- list()
top_n_genes <- 250

for (numb in x) {
  
  genes <- mega_heatmap_up %>%
    as.data.frame() %>%
    mutate(
      gene = rownames(.),
      max_rp = do.call(pmax, c(dplyr::select(., ends_with(numb)), na.rm = TRUE))
    ) %>%
    filter(max_rp > 0) %>%                 # drop true zeros only
    arrange(desc(max_rp)) %>%
    slice_head(n = top_n_genes) %>%         # relative, per-group cutoff
    pull(gene)
  
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  gene_list[[group_labels[[numb]]]] <- ids$ENTREZID
  
  background_expressed_genes <- switch(
    numb,
    "RETN_1" = results$lr_spp1_TCD8_up$background_expressed_genes,
    "RETN_2" = results$lr_spp1_TCD4_up$background_expressed_genes,
    "APP_1" = results$lr_spp1_TCD8_up$background_expressed_genes,
    "APP_2" = results$lr_spp1_TCD4_up$background_expressed_genes
  )
  
  background_list[[group_labels[[numb]]]] <- background_expressed_genes
}

# shared universe: union of all three backgrounds, converted once
all_background_symbols <- unique(unlist(background_list))

universe_ids <- bitr(
  all_background_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)$ENTREZID

# compareCluster across the three groups
cc <- compareCluster(
  geneCluster   = gene_list,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  universe      = universe_ids,
  pAdjustMethod = "BH")

cc_simplified <- clusterProfiler::simplify(cc)

# plot
dotplot(cc_simplified, showCategory = 20) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) + ylab(NULL)


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


##############################
# Make Figure
##############################

# read in previous figure half
a5 <- readRDS("../figures/a5.rds")

# Figure 5
ggdraw() +
  draw_plot(a5, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(b5, x = 0, y = .25, width = 1, height = .25) +
  draw_plot(c5, x = 0, y = 0, width = .5, height = .25) +
  draw_plot(d5, x = .5, y = 0, width = .5, height = .25) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15, 
                  x = c(0, 0, 0, .5), 
                  y = c(1, .5, .25, .25))
# save figure
ggsave("../figures/Figure_5_nichenet.jpg", width = 40, height = 50, units = c("cm"), dpi = 300)

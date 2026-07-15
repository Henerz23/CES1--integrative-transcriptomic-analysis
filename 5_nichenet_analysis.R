###########################################################
# MACROPHAGE TO CD8/CD4 NICHENET GENE EXPRESSION ANALYSIS
###########################################################

# Part 1: a function will be used to loop over the distinct LR networks to generate distinct gene expression heatmaps for each 
# Part 2 data exploration and figure creation


###########################################################
# Part 1
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
lr_networks <-  readRDS("../data/lr_networks.rds")
pbmc <- readRDS("../data/pbmc_final.rds")

print("data loaded")

# define targets
targets.cd8 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD8"])
targets.cd4 <- unique(pbmc@meta.data$cell_type_detailed[pbmc@meta.data$clMidwayPr == "TCD4"])

# read in the seurat object
seurat_obj <- subset(pbmc, cell_type_detailed %in% c(targets.cd4, targets.cd8, "Myeloid: Macro-C1QC", "Myeloid: Macro-SPP1"))
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
  
  if (network_name == "lr_spp1_TCD4") {
    sender <- "Myeloid: Macro-SPP1"
    receiver <- "TCD4"
  } else if (network_name == "lr_spp1_TCD8") {
    sender <- "Myeloid: Macro-SPP1"
    receiver <- "TCD8"
  } else if (network_name == "lr_c1qc_TCD4") {
    sender <- "Myeloid: Macro-C1QC"
    receiver <- "TCD4"
  } else if (network_name == "lr_c1qc_TCD8") {
    sender <- "Myeloid: Macro-C1QC"
    receiver <- "TCD8"
  } else if (network_name == "lr_both_TCD4") {
    sender <- c("Myeloid: Macro-SPP1",
                "Myeloid: Macro-C1QC")
    receiver <- "TCD4"
  } else if (network_name == "lr_both_TCD8") {
    sender <- c("Myeloid: Macro-SPP1",
                "Myeloid: Macro-C1QC")
    receiver <- "TCD8"
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
for (net_name in names(lr_networks)) {
  net <- lr_networks[[net_name]]
  
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

saveRDS(results, file = "../data/nichenet_results_300.rds")

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
lr_networks <- readRDS("../data/lr_networks.rds")

# load in the results from nichenet
# first run- less genes
results <- readRDS("../data/nichenet_results.rds")
# second run- 250 genes- then remove lower regulating ones
results <- readRDS("../data/nichenet_results_250.rds")
# third and final run
results <- readRDS("../data/nichenet_results_300.rds")


##############################
# visualise results- differences between priority genes
##############################

# create a list of all up/downregulated data to combine into a single figure
View(results)
# upregulated 
mat_list <- list(
  lr_spp1_TCD8 = results$lr_spp1_TCD8_up$gene_expression,
  lr_spp1_TCD4 = results$lr_spp1_TCD4_up$gene_expression,
  lr_both_TCD8 = results$lr_both_TCD8_up$gene_expression,
  lr_both_TCD4 = results$lr_both_TCD4_up$gene_expression,
  # lr_c1qc_TCD8 = results$lr_c1qc_TCD8_up$gene_expression,
  lr_c1qc_TCD4 = results$lr_c1qc_TCD4_up$gene_expression
)

# downregulated
# fix for ones with no columns
mat_list <- list(
   lr_spp1_TCD8 = results$lr_spp1_TCD8_down$gene_expression,
   lr_spp1_TCD4 = results$lr_spp1_TCD4_down$gene_expression,
   lr_both_TCD8 = results$lr_both_TCD8_down$gene_expression,
   lr_both_TCD4 = results$lr_both_TCD4_down$gene_expression,
  # lr_c1qc_TCD8 = results$lr_c1qc_TCD8_down$gene_expression,
   lr_c1qc_TCD4 = results$lr_c1qc_TCD4_down$gene_expression
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
keep <- apply(mega_heatmap_up, 1, function(x) max(x) - min(x) >= 0.05)
mega_heatmap_up <- mega_heatmap_up[keep, , drop = FALSE]
nrow(mega_heatmap_up)


# visualise the mega matrix ligand target network
b5 <- mega_heatmap_up %>%
  make_heatmap_ggplot(
    y_name = paste("Target genes in "),
    x_name = paste("prioritized ligands "),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) +
  coord_flip()
b5
# ggsave("../figures/nichenet_check.jpg", width = 80, height = 40, units = c("cm"), dpi = 300)

# downreg
mega_heatmap_up %>%
  make_heatmap_ggplot(
    y_name = paste("Target genes in "),
    x_name = paste("prioritized ligands "),
    color = "purple",
    legend_title = "Regulatory\npotential"
  ) +
  theme(axis.text.x = element_text(face = "italic")) +
  coord_flip()

##############################
# identify genes expressed through specific typesof interaction
##############################


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
saveRDS(groups, file = "../data/nichenet_gene_groups.rds")

##############################
# Make Figure
##############################
lr_networks9
# Figure 5
ggdraw() +
  draw_plot(a5, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(b5, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B"), size = 15, 
                  x = c(0, 0), 
                  y = c(1, .5))
# save figure
ggsave("../figures/Figure_5_nichenet_general.jpg", width = 40, height = 40, units = c("cm"), dpi = 300)

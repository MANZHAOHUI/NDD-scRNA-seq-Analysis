#!/usr/bin/env Rscript
# ==============================================================================
# Step 5: UMAP Dimensionality Reduction and Clustering
# ==============================================================================
# Description: Perform UMAP and identify cell clusters
# Input: Integrated Seurat objects from Step 4
# Output: Clustered Seurat objects with UMAP coordinates
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# Function: Perform UMAP
# ------------------------------------------------------------------------------
perform_umap <- function(integrated_obj, dims = 1:30,
                        n_neighbors = 30, min_dist = 0.6) {
  message("\n========== Performing UMAP Reduction ==========")
  message("  Parameters: dims=1:", max(dims),
          ", n_neighbors=", n_neighbors,
          ", min_dist=", min_dist)
  
  # Check if PCA exists
  if (!"pca" %in% names(integrated_obj@reductions)) {
    stop("Please run PCA first (Step 4)")
  }
  
  # Run UMAP
  integrated_obj <- RunUMAP(
    integrated_obj,
    dims = dims,
    reduction = "pca",
    assay = DefaultAssay(integrated_obj),
    umap.method = "uwot",
    n.neighbors = as.integer(n_neighbors),
    n.components = 2L,
    metric = "cosine",
    n.epochs = 1000,
    learning.rate = 1,
    min.dist = min_dist,
    spread = 2,
    set.op.mix.ratio = 1,
    local.connectivity = 1L,
    repulsion.strength = 1,
    negative.sample.rate = 5L,
    uwot.sgd = FALSE,
    seed.use = 42L,
    verbose = TRUE,
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  
  message("  UMAP complete")
  
  return(integrated_obj)
}

# ------------------------------------------------------------------------------
# Function: Perform Clustering
# ------------------------------------------------------------------------------
perform_clustering <- function(integrated_obj, resolution_range = seq(0.1, 1.5, 0.1)) {
  message("\n========== Building SNN Graph ==========")
  
  # Build k-nearest neighbor graph
  integrated_obj <- FindNeighbors(
    integrated_obj,
    reduction = "pca",
    dims = 1:30,
    k.param = 30,
    prune.SNN = 1/15,
    nn.method = "annoy",
    n.trees = 50,
    annoy.metric = "euclidean",
    verbose = FALSE
  )
  
  message("========== Multi-resolution Clustering ==========")
  
  # Store clustering results at different resolutions
  clustering_results <- list()
  
  for (res in resolution_range) {
    message("  Resolution: ", res)
    
    integrated_obj <- FindClusters(
      integrated_obj,
      resolution = res,
      verbose = FALSE
    )
    
    # Store results
    cluster_col <- paste0("RNA_snn_res.", res)
    n_clusters <- length(unique(integrated_obj@meta.data[[cluster_col]]))
    
    clustering_results[[as.character(res)]] <- list(
      resolution = res,
      n_clusters = n_clusters,
      clusters = integrated_obj@meta.data[[cluster_col]]
    )
    
    message("    Clusters: ", n_clusters)
  }
  
  # Select optimal resolution
  optimal_res <- select_optimal_resolution(clustering_results)
  
  # Set default clustering
  Idents(integrated_obj) <- integrated_obj@meta.data[[paste0("RNA_snn_res.", optimal_res)]]
  integrated_obj$seurat_clusters <- Idents(integrated_obj)
  
  message("  Selected resolution: ", optimal_res, " (",
          length(unique(integrated_obj$seurat_clusters)), " clusters)")
  
  return(integrated_obj)
}

# ------------------------------------------------------------------------------
# Function: Select Optimal Resolution
# ------------------------------------------------------------------------------
select_optimal_resolution <- function(clustering_results) {
  # Heuristic: prefer 10-20 clusters
  scores <- sapply(clustering_results, function(x) {
    n_clusters <- x$n_clusters
    
    if (n_clusters >= 10 && n_clusters <= 20) {
      score <- 1 / abs(n_clusters - 15)
    } else {
      score <- 0.1
    }
    
    return(score)
  })
  
  optimal_idx <- which.max(scores)
  optimal_res <- as.numeric(names(clustering_results)[optimal_idx])
  
  return(optimal_res)
}

# ------------------------------------------------------------------------------
# Function: Generate UMAP Visualizations
# ------------------------------------------------------------------------------
generate_umap_visualizations <- function(integrated_obj, comparison = "AD.NORMAL") {
  message("\n========== Generating UMAP Visualizations ==========")
  
  output_dir <- file.path(PATHS$figures, "UMAP", comparison)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Basic UMAP (by clusters)
  p1 <- DimPlot(
    integrated_obj,
    reduction = "umap",
    label = TRUE,
    label.size = 6,
    pt.size = 0.1
  ) +
    labs(title = "UMAP by Clusters") +
    theme(legend.position = "none")
  
  ggsave(
    file.path(output_dir, "UMAP_clusters.png"),
    p1, width = 10, height = 8, dpi = 300
  )
  
  # 2. By condition
  p2 <- DimPlot(
    integrated_obj,
    reduction = "umap",
    group.by = "orig.ident",
    pt.size = 0.1
  ) +
    labs(title = "UMAP by Condition")
  
  ggsave(
    file.path(output_dir, "UMAP_condition.png"),
    p2, width = 10, height = 8, dpi = 300
  )
  
  # 3. By sample
  p3 <- DimPlot(
    integrated_obj,
    reduction = "umap",
    group.by = "sampID",
    pt.size = 0.1
  ) +
    labs(title = "UMAP by Sample") +
    theme(legend.position = "bottom")
  
  ggsave(
    file.path(output_dir, "UMAP_sample.png"),
    p3, width = 12, height = 8, dpi = 300
  )
  
  # 4. Split by condition
  p4 <- DimPlot(
    integrated_obj,
    reduction = "umap",
    split.by = "orig.ident",
    pt.size = 0.1
  ) +
    labs(title = "UMAP Split by Condition")
  
  ggsave(
    file.path(output_dir, "UMAP_split_condition.png"),
    p4, width = 15, height = 6, dpi = 300
  )
  
  # 5. Feature plots for marker genes
  DefaultAssay(integrated_obj) <- "SCT"
  
  marker_genes <- c(
    "SLC1A2",   # Astrocyte
    "RBFOX3",   # Neuron (NeuN)
    "SLC17A7",  # Excitatory neuron
    "GAD1",     # Inhibitory neuron
    "AIF1",     # Microglia
    "MOBP"      # Oligodendrocyte
  )
  
  # Check which genes are available
  available_genes <- marker_genes[marker_genes %in% rownames(integrated_obj)]
  
  if (length(available_genes) > 0) {
    p5 <- FeaturePlot(
      integrated_obj,
      features = available_genes,
      reduction = "umap",
      pt.size = 0.1,
      ncol = 3
    )
    
    ggsave(
      file.path(output_dir, "UMAP_markers.png"),
      p5, width = 15, height = 10, dpi = 300
    )
  }
  
  message("  Visualizations saved to: ", output_dir)
  
  DefaultAssay(integrated_obj) <- "integrated"
}

# ------------------------------------------------------------------------------
# Function: Generate Clustering Statistics
# ------------------------------------------------------------------------------
generate_clustering_report <- function(integrated_obj, comparison = "AD.NORMAL") {
  message("\n========== Generating Clustering Report ==========")
  
  # Cluster statistics
  cluster_stats <- table(integrated_obj$seurat_clusters)
  
  # Condition distribution per cluster
  cluster_condition <- table(
    integrated_obj$seurat_clusters,
    integrated_obj$orig.ident
  )
  
  # QC metrics per cluster
  qc_stats <- integrated_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      n_cells = n(),
      mean_genes = mean(nFeature_RNA),
      mean_umi = mean(nCount_RNA),
      mean_mt = mean(percent.mt),
      .groups = 'drop'
    )
  
  # Save statistics
  output_dir <- PATHS$integration_results
  
  write.csv(
    qc_stats,
    file.path(output_dir, paste0(comparison, "_cluster_qc_stats.csv")),
    row.names = FALSE
  )
  
  write.csv(
    cluster_condition,
    file.path(output_dir, paste0(comparison, "_cluster_condition_distribution.csv"))
  )
  
  # Generate text report
  report_file <- file.path(output_dir, paste0(comparison, "_clustering_report.txt"))
  
  report_lines <- c(
    "========== Clustering Analysis Report ==========",
    paste("Comparison:", comparison),
    paste("Analysis time:", Sys.time()),
    "",
    "Clustering Statistics:",
    paste("  Total clusters:", length(unique(integrated_obj$seurat_clusters))),
    paste("  Largest cluster:", max(cluster_stats), "cells"),
    paste("  Smallest cluster:", min(cluster_stats), "cells"),
    paste("  Mean cluster size:", round(mean(cluster_stats), 2), "cells"),
    "",
    "UMAP Parameters:",
    "  n_neighbors: 30",
    "  min_dist: 0.6",
    "  metric: cosine",
    "",
    "Clustering Parameters:",
    "  Algorithm: Louvain",
    "  k.param: 30"
  )
  
  writeLines(report_lines, report_file)
  message("  Clustering report saved: ", report_file)
}

# ------------------------------------------------------------------------------
# Function: Main UMAP/Clustering Pipeline
# ------------------------------------------------------------------------------
run_umap_clustering <- function(comparison = "AD.NORMAL") {
  message("\n", strrep("=", 70))
  message("Starting UMAP and Clustering: ", comparison)
  message(strrep("=", 70))
  
  # Read integrated object
  input_file <- file.path(
    PATHS$integration_results,
    comparison,
    paste0(comparison, ".integrated_PCA.rds")
  )
  
  if (!file.exists(input_file)) {
    stop("Integrated object not found. Please run Step 4 first.")
  }
  
  integrated_obj <- readRDS(input_file)
  
  # Perform UMAP
  integrated_obj <- perform_umap(integrated_obj)
  
  # Perform clustering
  integrated_obj <- perform_clustering(integrated_obj)
  
  # Generate visualizations
  generate_umap_visualizations(integrated_obj, comparison)
  
  # Generate report
  generate_clustering_report(integrated_obj, comparison)
  
  # Save results
  output_dir <- file.path(PATHS$integration_results, comparison)
  output_file <- file.path(output_dir, paste0(comparison, ".integrated_PCA_UMAP.rds"))
  
  saveRDS(integrated_obj, output_file)
  
  message("\n", strrep("=", 70))
  message("✓ Results saved: ", output_file)
  message(strrep("=", 70))
  
  return(integrated_obj)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message(strrep("=", 70))
  message("Step 5: UMAP Dimensionality Reduction and Clustering")
  message(strrep("=", 70))
  
  # Process each comparison
  comparisons <- c("AD.NORMAL", "PD.NORMAL", "DLB.NORMAL")
  
  for (comparison in comparisons) {
    tryCatch({
      run_umap_clustering(comparison)
    }, error = function(e) {
      warning("Failed to process ", comparison, ": ", e$message)
    })
  }
  
  message("\n", strrep("=", 70))
  message("✓ Step 5 completed successfully!")
  message(strrep("=", 70))
}

# Run if executed as script
if (!interactive()) {
  main()
}

#!/usr/bin/env Rscript
# ==============================================================================
# Step 4: Sample Integration and Batch Correction
# ==============================================================================
# Description: Merge samples and correct batch effects using SCTransform
# Input: Filtered Seurat objects from Step 3
# Output: Integrated Seurat objects with PCA
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(sctransform)
  library(glmGamPoi)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Function: Prepare Integration Data
# ------------------------------------------------------------------------------
prepare_integration_data <- function(comparison = "AD.NORMAL") {
  message("\n========== Preparing integration data: ", comparison, " ==========")
  
  # Define sample groups
  groups <- strsplit(comparison, "\\.")[[1]]
  message("Groups: ", paste(groups, collapse = " vs "))
  
  # Read doublet-filtered metadata
  doublet_results <- list.files(
    PATHS$qc_results,
    pattern = "doublet_finder_.*\\.rds",
    full.names = TRUE
  )
  
  # Build sample list
  sample_list <- list()
  
  for (result_file in doublet_results) {
    result <- readRDS(result_file)
    sample_id <- result$sample_id
    
    # Read original Seurat object
    seurat_file <- file.path(PATHS$seurat_objects, paste0(sample_id, ".rds"))
    
    if (!file.exists(seurat_file)) {
      warning("Seurat object not found: ", seurat_file)
      next
    }
    
    seurat_obj <- readRDS(seurat_file)
    
    # Update metadata with doublet information
    metadata_updated <- result$metadata[colnames(seurat_obj), ]
    seurat_obj@meta.data <- metadata_updated
    
    # Filter doublets
    df_col <- grep("^DF.classifications", colnames(metadata_updated), value = TRUE)[1]
    singlet_cells <- rownames(metadata_updated)[metadata_updated[[df_col]] == "Singlet"]
    seurat_obj <- subset(seurat_obj, cells = singlet_cells)
    
    # Check if sample belongs to comparison groups
    if (seurat_obj@meta.data$group[1] %in% groups) {
      sample_list[[sample_id]] <- seurat_obj
      message("  Added sample: ", sample_id, " (", ncol(seurat_obj), " cells)")
    }
  }
  
  if (length(sample_list) < 2) {
    stop("Need at least 2 samples for integration. Found: ", length(sample_list))
  }
  
  message("Total samples for integration: ", length(sample_list))
  
  return(sample_list)
}

# ------------------------------------------------------------------------------
# Function: SCT Integration
# ------------------------------------------------------------------------------
integrate_samples_sct <- function(sample_list, comparison = "AD.NORMAL",
                                  integration_features = 10000) {
  message("\n========== SCT Integration: ", comparison, " ==========")
  
  if (length(sample_list) < 2) {
    stop("At least 2 samples required for integration")
  }
  
  # Merge objects
  message("Merging ", length(sample_list), " samples...")
  combined <- merge(
    x = sample_list[[1]],
    y = sample_list[-1],
    add.cell.ids = names(sample_list)
  )
  
  message("  Merged cells: ", ncol(combined))
  message("  Genes: ", nrow(combined))
  
  # Remove mitochondrial genes to reduce technical noise
  mt_genes <- grep("^MT-", rownames(combined), value = TRUE)
  message("  Removing ", length(mt_genes), " mitochondrial genes")
  combined <- subset(combined, features = setdiff(rownames(combined), mt_genes))
  
  # Split by sample
  combined.list <- SplitObject(combined, split.by = "sampID")
  rm(combined)
  gc()
  
  # SCT normalization for each sample
  message("\nPerforming SCTransform normalization...")
  for (i in seq_along(combined.list)) {
    message("  Processing sample ", i, "/", length(combined.list), ": ",
            names(combined.list)[i])
    
    combined.list[[i]] <- tryCatch({
      SCTransform(
        combined.list[[i]],
        method = 'glmGamPoi',
        return.only.var.genes = FALSE,
        verbose = FALSE,
        vst.flavor = "v2",
        conserve.memory = TRUE
      )
    }, error = function(e) {
      warning("SCTransform failed (", names(combined.list)[i], "): ", e$message)
      return(NULL)
    })
  }
  
  # Remove failed samples
  combined.list <- combined.list[!sapply(combined.list, is.null)]
  
  # Filter small samples
  combined.list <- combined.list[sapply(combined.list, ncol) > 100]
  message("Samples after SCT: ", length(combined.list))
  
  # Select integration features
  message("\nSelecting integration features...")
  features <- SelectIntegrationFeatures(
    object.list = combined.list,
    nfeatures = integration_features
  )
  message("  Selected ", length(features), " features")
  
  # Prepare SCT integration
  combined.list <- PrepSCTIntegration(
    object.list = combined.list,
    anchor.features = features
  )
  
  # Select reference dataset (largest sample)
  cell_counts <- sapply(combined.list, ncol)
  reference_idx <- which.max(cell_counts)
  message("  Reference dataset: ", names(combined.list)[reference_idx],
          " (", cell_counts[reference_idx], " cells)")
  
  # Find integration anchors
  message("\nFinding integration anchors...")
  min_cells <- min(cell_counts)
  k_anchor <- min(5, min_cells - 1)
  
  anchors <- FindIntegrationAnchors(
    object.list = combined.list,
    normalization.method = "SCT",
    anchor.features = features,
    reference = reference_idx,
    k.anchor = k_anchor,
    verbose = FALSE
  )
  
  message("  Found ", length(anchors@anchors[, 1]), " anchors")
  
  # Integrate data
  message("\nIntegrating data...")
  combined.integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    verbose = FALSE
  )
  
  # Clean up memory
  rm(anchors, combined.list)
  gc()
  
  message("Integration complete!")
  
  return(combined.integrated)
}

# ------------------------------------------------------------------------------
# Function: PCA Analysis
# ------------------------------------------------------------------------------
perform_pca_analysis <- function(integrated_obj, npcs = 50) {
  message("\n========== Performing PCA Analysis ==========")
  
  # Run PCA
  integrated_obj <- RunPCA(
    integrated_obj,
    npcs = npcs,
    verbose = FALSE
  )
  
  # Calculate variance explained
  pct_var <- integrated_obj@reductions$pca@stdev^2 /
    sum(integrated_obj@reductions$pca@stdev^2) * 100
  
  # Determine significant PCs (elbow method)
  var_change <- diff(pct_var)
  elbow_point <- which(abs(var_change) < 0.1)[1]
  
  if (is.na(elbow_point)) {
    elbow_point <- 30
  }
  
  message("  Recommended PCs: ", elbow_point)
  message("  Variance explained by PC1-30: ", round(sum(pct_var[1:30]), 2), "%")
  
  return(integrated_obj)
}

# ------------------------------------------------------------------------------
# Function: Generate PCA Diagnostic Plots
# ------------------------------------------------------------------------------
generate_pca_plots <- function(integrated_obj, comparison = "AD.NORMAL") {
  message("\nGenerating PCA diagnostic plots...")
  
  output_dir <- file.path(PATHS$figures, "PCA")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate variance explained
  pct_var <- integrated_obj@reductions$pca@stdev^2 /
    sum(integrated_obj@reductions$pca@stdev^2) * 100
  
  # 1. Elbow plot
  elbow_data <- data.frame(
    PC = 1:length(pct_var),
    Variance = pct_var
  )
  
  p1 <- ggplot(elbow_data, aes(x = PC, y = Variance)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
    labs(
      title = "PCA Elbow Plot",
      x = "Principal Component",
      y = "% Variance Explained"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    file.path(output_dir, paste0(comparison, "_elbow_plot.png")),
    p1, width = 8, height = 6, dpi = 300
  )
  
  # 2. PCA by sample
  p2 <- DimPlot(
    integrated_obj,
    reduction = "pca",
    group.by = "sampID",
    pt.size = 0.1
  ) +
    labs(title = "PCA by Sample") +
    theme(legend.position = "right")
  
  ggsave(
    file.path(output_dir, paste0(comparison, "_pca_by_sample.png")),
    p2, width = 10, height = 8, dpi = 300
  )
  
  # 3. PCA by condition
  p3 <- DimPlot(
    integrated_obj,
    reduction = "pca",
    group.by = "orig.ident",
    pt.size = 0.1
  ) +
    labs(title = "PCA by Condition")
  
  ggsave(
    file.path(output_dir, paste0(comparison, "_pca_by_condition.png")),
    p3, width = 10, height = 8, dpi = 300
  )
  
  message("  PCA plots saved to: ", output_dir)
}

# ------------------------------------------------------------------------------
# Function: Generate Integration Report
# ------------------------------------------------------------------------------
generate_integration_report <- function(integrated_obj, comparison) {
  report_file <- file.path(
    PATHS$integration_results,
    paste0(comparison, "_integration_report.txt")
  )
  
  # Collect statistics
  n_cells <- ncol(integrated_obj)
  n_genes <- nrow(integrated_obj)
  n_samples <- length(unique(integrated_obj$sampID))
  condition_counts <- table(integrated_obj$orig.ident)
  
  # Calculate variance metrics
  pct_var <- integrated_obj@reductions$pca@stdev^2 /
    sum(integrated_obj@reductions$pca@stdev^2) * 100
  
  report_lines <- c(
    "========== Integration Analysis Report ==========",
    paste("Comparison:", comparison),
    paste("Analysis time:", Sys.time()),
    "",
    "Data Overview:",
    paste("  Total cells:", n_cells),
    paste("  Genes:", n_genes),
    paste("  Samples:", n_samples),
    "",
    "Condition Distribution:"
  )
  
  for (condition in names(condition_counts)) {
    report_lines <- c(
      report_lines,
      paste("  ", condition, ":", condition_counts[condition], "cells")
    )
  }
  
  report_lines <- c(
    report_lines,
    "",
    "Integration Parameters:",
    "  Normalization method: SCTransform",
    "  Integration features: 10000",
    "  Anchor algorithm: CCA",
    "",
    "Quality Metrics:",
    paste("  PC1 variance explained:", round(pct_var[1], 2), "%"),
    paste("  PC1-30 cumulative variance:", round(sum(pct_var[1:30]), 2), "%")
  )
  
  writeLines(report_lines, report_file)
  message("  Integration report saved: ", report_file)
}

# ------------------------------------------------------------------------------
# Function: Main Integration Pipeline
# ------------------------------------------------------------------------------
run_integration_pipeline <- function(comparison = "AD.NORMAL") {
  message("\n", strrep("=", 70))
  message("Starting Integration Pipeline: ", comparison)
  message(strrep("=", 70))
  
  # Create output directory
  output_dir <- file.path(PATHS$integration_results, comparison)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data
  sample_list <- prepare_integration_data(comparison)
  
  if (length(sample_list) < 2) {
    stop("Insufficient samples for integration")
  }
  
  # SCT integration
  integrated_obj <- integrate_samples_sct(sample_list, comparison)
  
  # PCA analysis
  integrated_obj <- perform_pca_analysis(integrated_obj)
  
  # Generate plots
  generate_pca_plots(integrated_obj, comparison)
  
  # Generate report
  generate_integration_report(integrated_obj, comparison)
  
  # Save results
  output_file <- file.path(output_dir, paste0(comparison, ".integrated_PCA.rds"))
  saveRDS(integrated_obj, output_file)
  
  message("\n", strrep("=", 70))
  message("✓ Integration complete: ", output_file)
  message(strrep("=", 70))
  
  return(integrated_obj)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message(strrep("=", 70))
  message("Step 4: Sample Integration and Batch Correction")
  message(strrep("=", 70))
  
  # Process each comparison
  comparisons <- c("AD.NORMAL", "PD.NORMAL", "DLB.NORMAL")
  
  for (comparison in comparisons) {
    tryCatch({
      run_integration_pipeline(comparison)
    }, error = function(e) {
      warning("Failed to process ", comparison, ": ", e$message)
    })
  }
  
  message("\n", strrep("=", 70))
  message("✓ Step 4 completed successfully!")
  message(strrep("=", 70))
}

# Run if executed as script
if (!interactive()) {
  main()
}

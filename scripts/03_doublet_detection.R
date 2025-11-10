#!/usr/bin/env Rscript
# ==============================================================================
# Step 3: Doublet Detection and Removal
# ==============================================================================
# Description: Identify and remove doublets using DoubletFinder
# Input: Seurat objects from Step 2
# Output: Filtered Seurat objects and QC reports
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(DoubletFinder)
  library(sctransform)
  library(glmGamPoi)
  library(ggplot2)
})

set.seed(7)  # For reproducibility

# ------------------------------------------------------------------------------
# Function: Optimize DoubletFinder Parameters
# ------------------------------------------------------------------------------
optimize_doubletfinder_params <- function(seurat_obj, sample_id) {
  message('\nOptimizing DoubletFinder parameters for: ', sample_id)
  
  # SCTransform normalization
  seurat_obj <- SCTransform(
    seurat_obj,
    method = 'glmGamPoi',
    return.only.var.genes = TRUE,
    verbose = FALSE
  )
  
  # PCA dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, assay = 'SCT', verbose = FALSE)
  
  # Calculate number of PCs to use (explaining 90% variance)
  eig <- seurat_obj[['pca']]@stdev^2
  pct_var_explained <- eig / sum(eig)
  cumsum_var <- cumsum(pct_var_explained)
  min.pc <- max(20, which(cumsum_var >= 0.9)[1])
  
  message('  Using ', min.pc, ' principal components')
  
  # Initial clustering for homotypic doublet estimation
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min.pc, assay = 'SCT', verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:min.pc, assay = 'SCT', verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, assay = 'SCT', verbose = FALSE)
  
  # Parameter sweep for optimal pK
  message('  Running parameter sweep...')
  sweep_res_list <- list()
  
  for(ix in 1:10) {
    if (ix %% 2 == 0) message('    Iteration ', ix, '/10')
    
    # pK range
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, 0.01))
    
    # Filter valid pK values
    min.cells <- round(nrow(seurat_obj@meta.data)/(1-0.05) - nrow(seurat_obj@meta.data))
    pK <- pK[round(pK * min.cells) >= 1]
    
    # Subsample for speed if needed
    if(nrow(seurat_obj@meta.data) > 10000) {
      sample_cells <- sample(colnames(seurat_obj), 10000)
      data_subset <- seurat_obj@assays$RNA@counts[, sample_cells]
    } else {
      data_subset <- seurat_obj@assays$RNA@counts
    }
    
    # Run parameter sweep
    sweep_res <- paramSweep_v3(seurat_obj, PCs = 1:min.pc, sct = TRUE)
    sweep_res_list[[ix]] <- sweep_res
  }
  
  # Summarize results
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  # Select optimal pK
  optimal.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message('  Optimal pK: ', optimal.pk)
  
  return(list(
    seurat_obj = seurat_obj,
    optimal_pk = optimal.pk,
    min_pc = min.pc
  ))
}

# ------------------------------------------------------------------------------
# Function: Detect Doublets
# ------------------------------------------------------------------------------
detect_doublets <- function(sample_id, expected_doublet_rate = 0.065) {
  message("\n========== Processing sample: ", sample_id, " ==========")
  
  # Read Seurat object
  input_file <- file.path(PATHS$seurat_objects, paste0(sample_id, '.rds'))
  
  if (!file.exists(input_file)) {
    warning("File not found: ", input_file)
    return(NULL)
  }
  
  seurat_obj <- readRDS(input_file)
  cells_before <- ncol(seurat_obj)
  message("Original cell count: ", cells_before)
  
  # Optimize parameters
  opt_result <- optimize_doubletfinder_params(seurat_obj, sample_id)
  seurat_obj <- opt_result$seurat_obj
  optimal_pk <- opt_result$optimal_pk
  min_pc <- opt_result$min_pc
  
  # Estimate homotypic doublet proportion
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # Adjust expected doublet number
  nExp.poi <- round(expected_doublet_rate * nrow(seurat_obj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  message("Expected doublets: ", nExp.poi, " (adjusted: ", nExp.poi.adj, ")")
  
  # Run DoubletFinder
  message("Running DoubletFinder...")
  seurat_obj <- doubletFinder_v3(
    seurat_obj,
    PCs = 1:min_pc,
    pN = 0.25,
    pK = optimal_pk,
    nExp = nExp.poi.adj,
    sct = TRUE
  )
  
  # Extract results
  df_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)[1]
  doublet_calls <- seurat_obj@meta.data[[df_col]]
  
  # Calculate statistics
  n_doublets <- sum(doublet_calls == "Doublet")
  doublet_rate <- n_doublets / cells_before * 100
  
  message("Detected doublets: ", n_doublets, " (", round(doublet_rate, 2), "%)")
  
  # Save results
  result <- list(
    sample_id = sample_id,
    cells_total = cells_before,
    doublets = n_doublets,
    doublet_rate = doublet_rate,
    optimal_pk = optimal_pk,
    metadata = seurat_obj@meta.data,
    seurat_obj = seurat_obj
  )
  
  output_file <- file.path(PATHS$qc_results, paste0('doublet_finder_', sample_id, '.rds'))
  saveRDS(result, output_file)
  
  # Generate report
  generate_doublet_report(result, sample_id)
  
  # Generate visualization
  generate_doublet_plots(seurat_obj, sample_id, df_col)
  
  return(result)
}

# ------------------------------------------------------------------------------
# Function: Generate Doublet Report
# ------------------------------------------------------------------------------
generate_doublet_report <- function(result, sample_id) {
  report_file <- file.path(PATHS$qc_results, paste0('doublet_report_', sample_id, '.txt'))
  
  report_lines <- c(
    "========== Doublet Detection Report ==========",
    paste("Sample ID:", sample_id),
    paste("Analysis time:", Sys.time()),
    "",
    "Detection Results:",
    paste("  Total cells:", result$cells_total),
    paste("  Detected doublets:", result$doublets),
    paste("  Doublet rate:", round(result$doublet_rate, 2), "%"),
    paste("  Optimal pK value:", result$optimal_pk),
    "",
    "Method Parameters:",
    "  Algorithm: DoubletFinder v3",
    "  Normalization: SCTransform",
    "  Expected doublet rate: 6.5%",
    "  pN: 0.25",
    "",
    "Quality Assessment:",
    if (result$doublet_rate < 2) "  WARNING: Doublet rate is unusually low" 
    else if (result$doublet_rate > 15) "  WARNING: Doublet rate is unusually high"
    else "  Doublet rate is within normal range (5-10%)"
  )
  
  writeLines(report_lines, report_file)
  message("  Report saved: ", report_file)
}

# ------------------------------------------------------------------------------
# Function: Generate Doublet Visualizations
# ------------------------------------------------------------------------------
generate_doublet_plots <- function(seurat_obj, sample_id, df_col) {
  output_dir <- file.path(PATHS$figures, "doublets")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. UMAP with doublet labels
  p1 <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = df_col,
    pt.size = 0.5
  ) +
    ggtitle(paste("Doublet Detection -", sample_id)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    file.path(output_dir, paste0(sample_id, "_doublet_umap.png")),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # 2. pANN score distribution
  pann_col <- grep("^pANN", colnames(seurat_obj@meta.data), value = TRUE)[1]
  
  if (!is.null(pann_col)) {
    p2 <- ggplot(seurat_obj@meta.data, aes(x = .data[[pann_col]])) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(
        xintercept = quantile(seurat_obj@meta.data[[pann_col]], 0.95),
        linetype = "dashed",
        color = "red"
      ) +
      labs(
        title = paste("pANN Score Distribution -", sample_id),
        x = "pANN Score",
        y = "Count"
      ) +
      theme_minimal()
    
    ggsave(
      file.path(output_dir, paste0(sample_id, "_pann_distribution.png")),
      plot = p2,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  message("  Visualizations saved to: ", output_dir)
}

# ------------------------------------------------------------------------------
# Function: Batch Process All Samples
# ------------------------------------------------------------------------------
batch_process_doublets <- function(sample_ids = NULL) {
  message("\n========== Starting Batch Doublet Detection ==========")
  
  # Get all sample IDs if not provided
  if (is.null(sample_ids)) {
    seurat_files <- list.files(
      PATHS$seurat_objects,
      pattern = "\\.rds$",
      full.names = FALSE
    )
    sample_ids <- gsub("\\.rds$", "", seurat_files)
  }
  
  message("Processing ", length(sample_ids), " samples")
  
  # Process each sample
  results <- list()
  for (sample_id in sample_ids) {
    result <- tryCatch({
      detect_doublets(sample_id)
    }, error = function(e) {
      warning("Failed to process ", sample_id, ": ", e$message)
      NULL
    })
    
    if (!is.null(result)) {
      results[[sample_id]] <- result
    }
  }
  
  # Generate summary statistics
  summary_stats <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      sample = x$sample_id,
      cells = x$cells_total,
      doublets = x$doublets,
      rate = x$doublet_rate
    )
  }))
  
  # Save summary
  write.csv(
    summary_stats,
    file.path(PATHS$qc_results, "doublet_summary.csv"),
    row.names = FALSE
  )
  
  message("\n========== Doublet Detection Complete ==========")
  message("Average doublet rate: ", round(mean(summary_stats$rate), 2), "%")
  message("Summary saved to: ", file.path(PATHS$qc_results, "doublet_summary.csv"))
  
  return(results)
}

# ------------------------------------------------------------------------------
# Function: Filter Doublets from Seurat Object
# ------------------------------------------------------------------------------
filter_doublets <- function(sample_id) {
  message("Filtering doublets from: ", sample_id)
  
  # Load doublet results
  doublet_file <- file.path(PATHS$qc_results, paste0('doublet_finder_', sample_id, '.rds'))
  if (!file.exists(doublet_file)) {
    stop("Doublet results not found for ", sample_id)
  }
  
  result <- readRDS(doublet_file)
  
  # Load original Seurat object
  seurat_file <- file.path(PATHS$seurat_objects, paste0(sample_id, '.rds'))
  seurat_obj <- readRDS(seurat_file)
  
  # Update metadata
  df_col <- grep("^DF.classifications", colnames(result$metadata), value = TRUE)[1]
  
  # Filter singlets only
  singlet_cells <- rownames(result$metadata)[result$metadata[[df_col]] == "Singlet"]
  seurat_obj <- subset(seurat_obj, cells = singlet_cells)
  
  message("  Cells retained: ", ncol(seurat_obj))
  
  # Save filtered object
  output_file <- file.path(PATHS$seurat_objects, paste0(sample_id, '_filtered.rds'))
  saveRDS(seurat_obj, output_file)
  
  message("  Filtered object saved: ", output_file)
  
  return(seurat_obj)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message(strrep("=", 70))
  message("Step 3: Doublet Detection and Removal")
  message(strrep("=", 70))
  
  # Run batch processing
  results <- batch_process_doublets()
  
  # Filter doublets from all samples
  message("\n========== Filtering Doublets ==========")
  for (sample_id in names(results)) {
    filter_doublets(sample_id)
  }
  
  message("\n", strrep("=", 70))
  message("✓ Step 3 completed successfully!")
  message(strrep("=", 70))
}

# Run if executed as script
if (!interactive()) {
  main()
}

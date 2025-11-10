#!/usr/bin/env Rscript
# ==============================================================================
# Step 7: Differential Expression Analysis with NEBULA
# ==============================================================================
# Description: Identify disease-associated genes using NEBULA mixed models
# Input: Annotated Seurat objects from Step 6
# Output: DEG lists per cell type
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(data.table)
  library(Matrix.utils)
})

# Install NEBULA if needed
if (!requireNamespace("nebula", quietly = TRUE)) {
  message("Installing NEBULA...")
  install.packages("nebula")
}
library(nebula)

# ------------------------------------------------------------------------------
# Function: Prepare NEBULA Data
# ------------------------------------------------------------------------------
prepare_nebula_data <- function(seurat_obj, cell_type, comparison = "AD.NORMAL") {
  message("\nPreparing NEBULA data for: ", cell_type)
  
  # Extract specific cell type
  cell_subset <- subset(seurat_obj, mayo.cell.type.number == cell_type)
  
  if (ncol(cell_subset) < 50) {
    warning("Insufficient cells (", ncol(cell_subset), "), may affect power")
    if (ncol(cell_subset) < 20) {
      return(NULL)
    }
  }
  
  DefaultAssay(cell_subset) <- "RNA"
  
  # Gene filtering
  message("  Filtering genes...")
  data <- cell_subset@assays$RNA@counts
  
  # Remove completely unexpressed genes
  genes_expressed <- rowSums(data > 0) > 0
  genes_keep <- names(genes_expressed)[genes_expressed]
  
  message("    Original genes: ", nrow(data))
  message("    Retained genes: ", length(genes_keep))
  
  # Apply filtering
  cell_subset <- subset(cell_subset, features = genes_keep)
  
  return(cell_subset)
}

# ------------------------------------------------------------------------------
# Function: Build NEBULA Input
# ------------------------------------------------------------------------------
build_nebula_input <- function(cell_subset, comparison = "AD.NORMAL") {
  message("  Building NEBULA input matrix...")
  
  # Extract count matrix
  count_matrix <- cell_subset@assays$RNA@counts
  
  # Transpose to cell x gene format
  n_cells <- ncol(count_matrix)
  n_genes <- nrow(count_matrix)
  
  message("    Matrix dimensions: ", n_genes, " genes × ", n_cells, " cells")
  
  # Handle large matrices with chunking
  if (n_cells > 10000) {
    message("    Using chunked transpose for large matrix...")
    chunk_size <- 5000
    n_chunks <- ceiling(n_cells / chunk_size)
    transposed_chunks <- list()
    
    for (i in 1:n_chunks) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_cells)
      chunk <- count_matrix[, start_idx:end_idx]
      transposed_chunks[[i]] <- as.data.frame(as.matrix(t(chunk)))
      
      if (i %% 5 == 0) message("      Processing chunk ", i, "/", n_chunks)
    }
    
    allgenes <- bind_rows(transposed_chunks)
    rm(transposed_chunks)
    gc()
  } else {
    allgenes <- as.data.frame(as.matrix(t(count_matrix)))
  }
  
  # Add metadata
  metadata <- cell_subset@meta.data
  allgenes$orig.ident <- metadata$orig.ident
  allgenes$sampID <- metadata$sampID
  allgenes$Sex <- metadata$Sex
  allgenes$age <- metadata$age
  allgenes$PMI <- metadata$PMI
  allgenes$no.nuclei.after.filter <- as.numeric(metadata$no.nuclei.after.filter)
  allgenes$wellKey <- rownames(metadata)
  
  # Set row names
  rownames(allgenes) <- allgenes$wellKey
  
  # Reorder columns (metadata first, then genes)
  n_genes_col <- ncol(allgenes) - 7
  gene_cols <- 1:n_genes_col
  meta_cols <- (n_genes_col + 1):(n_genes_col + 7)
  allgenes <- allgenes[, c(meta_cols, gene_cols)]
  
  return(allgenes)
}

# ------------------------------------------------------------------------------
# Function: Filter Genes by Expression
# ------------------------------------------------------------------------------
filter_genes_by_expression <- function(allgenes, groups, min_pct = 0.1) {
  message("  Filtering genes by expression percentage...")
  
  # Separate two groups
  group1_data <- allgenes[allgenes$orig.ident == groups[1], ]
  group2_data <- allgenes[allgenes$orig.ident == groups[2], ]
  
  if (nrow(group1_data) == 0 || nrow(group2_data) == 0) {
    stop("One group has no cells")
  }
  
  message("    ", groups[1], ": ", nrow(group1_data), " cells")
  message("    ", groups[2], ": ", nrow(group2_data), " cells")
  
  # Calculate expression percentage
  calculate_pct <- function(data, gene_start_col = 8) {
    gene_data <- as.matrix(data[, gene_start_col:ncol(data)])
    pct <- apply(gene_data, 2, function(x) sum(x > 0) / length(x))
    return(pct)
  }
  
  pct1 <- calculate_pct(group1_data)
  pct2 <- calculate_pct(group2_data)
  
  # Keep genes expressed in ≥min_pct in at least one group
  genes_keep <- names(which(pmax(pct1, pct2) >= min_pct))
  
  message("    Retained genes: ", length(genes_keep))
  
  return(genes_keep)
}

# ------------------------------------------------------------------------------
# Function: Run NEBULA Differential Expression
# ------------------------------------------------------------------------------
run_nebula_de <- function(allgenes, genes_keep, comparison = "AD.NORMAL") {
  message("  Running NEBULA differential expression...")
  
  # Prepare data
  gene_data <- as.matrix(allgenes[, genes_keep])
  
  # Prepare covariates
  metadata <- allgenes[, 1:7]
  metadata$orig.ident <- as.factor(metadata$orig.ident)
  metadata$sampID <- as.factor(metadata$sampID)
  
  # Standardize continuous variables
  for (col in c("age", "PMI", "no.nuclei.after.filter")) {
    if (col %in% colnames(metadata)) {
      metadata[[col]] <- scale(metadata[[col]])
    }
  }
  
  # Build design matrix
  design_formula <- ~ Sex + age + PMI + no.nuclei.after.filter + orig.ident
  design_matrix <- model.matrix(design_formula, data = metadata)
  
  # Ensure data alignment
  gene_data <- gene_data[rownames(design_matrix), ]
  metadata <- metadata[rownames(design_matrix), ]
  
  # Extract sample ID for random effect
  sample_id <- metadata$sampID
  
  message("    Genes: ", ncol(gene_data))
  message("    Cells: ", nrow(gene_data))
  message("    Samples: ", length(unique(sample_id)))
  
  # Run NEBULA-HL
  nebula_result <- tryCatch({
    nebula(
      count = t(gene_data),  # genes × cells
      id = sample_id,        # random effect
      pred = design_matrix,  # fixed effects
      method = 'HL'          # H-likelihood method
    )
  }, error = function(e) {
    warning("NEBULA failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(nebula_result)) {
    message("    NEBULA complete")
  }
  
  return(nebula_result)
}

# ------------------------------------------------------------------------------
# Function: Process NEBULA Results
# ------------------------------------------------------------------------------
process_nebula_results <- function(nebula_result, comparison = "AD.NORMAL") {
  if (is.null(nebula_result)) {
    return(NULL)
  }
  
  message("  Processing NEBULA results...")
  
  result <- nebula_result$summary
  
  # Get contrast column names
  contrast_col <- grep("logFC_orig.ident", colnames(result), value = TRUE)[1]
  pval_col <- grep("p_orig.ident", colnames(result), value = TRUE)[1]
  
  # Calculate log2FC
  result$log2FC <- result[[contrast_col]] / log(2)
  
  # FDR correction
  result$fdr <- p.adjust(result[[pval_col]], method = 'fdr')
  
  # Select key columns
  result <- result[, c('gene', 'log2FC', 'fdr')]
  
  # Sort by FDR
  result <- result[order(result$fdr), ]
  
  # Filter significant genes
  sig_result <- result[result$fdr < 0.05 & abs(result$log2FC) >= 0.2, ]
  
  message("    Total genes: ", nrow(result))
  message("    Significant genes: ", nrow(sig_result))
  message("      Upregulated: ", sum(sig_result$log2FC > 0))
  message("      Downregulated: ", sum(sig_result$log2FC < 0))
  
  return(list(
    all = result,
    significant = sig_result
  ))
}

# ------------------------------------------------------------------------------
# Function: Batch Process All Cell Types
# ------------------------------------------------------------------------------
run_nebula_batch <- function(comparison = "AD.NORMAL") {
  message("\n", strrep("=", 70))
  message("NEBULA Batch Analysis: ", comparison)
  message(strrep("=", 70))
  
  # Read annotated object
  input_file <- file.path(
    PATHS$annotation_results,
    paste0(comparison, ".integrated_PCA_UMAP_annotated.rds")
  )
  
  if (!file.exists(input_file)) {
    stop("Please run Step 6 (cell annotation) first")
  }
  
  seurat_obj <- readRDS(input_file)
  
  # Get all cell types
  cell_types <- unique(seurat_obj$mayo.cell.type.number)
  message("Cell types: ", length(cell_types))
  
  # Create output directory
  output_dir <- file.path(PATHS$deg_results, comparison, "NEBULA")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Store results
  all_results <- list()
  failed_types <- c()
  
  # Process each cell type
  for (cell_type in cell_types) {
    message("\n>>> Processing: ", cell_type)
    
    # Prepare data
    cell_subset <- prepare_nebula_data(seurat_obj, cell_type, comparison)
    
    if (is.null(cell_subset)) {
      failed_types <- c(failed_types, cell_type)
      next
    }
    
    # Build input
    allgenes <- build_nebula_input(cell_subset, comparison)
    
    # Get groups
    groups <- unique(allgenes$orig.ident)
    
    # Filter genes
    genes_keep <- filter_genes_by_expression(allgenes, groups)
    
    # Run NEBULA
    nebula_result <- run_nebula_de(allgenes, genes_keep, comparison)
    
    # Process results
    processed_result <- process_nebula_results(nebula_result, comparison)
    
    if (!is.null(processed_result)) {
      all_results[[cell_type]] <- processed_result
      
      # Save results
      write.csv(
        processed_result$all,
        file.path(output_dir, paste0(comparison, ".", cell_type, "_all_genes_NEBULA.csv")),
        row.names = FALSE
      )
      
      write.csv(
        processed_result$significant,
        file.path(output_dir, paste0(comparison, ".", cell_type, "_sig_genes_NEBULA.csv")),
        row.names = FALSE
      )
    } else {
      failed_types <- c(failed_types, cell_type)
    }
  }
  
  # Generate summary report
  generate_nebula_report(all_results, failed_types, comparison)
  
  return(all_results)
}

# ------------------------------------------------------------------------------
# Function: Generate NEBULA Report
# ------------------------------------------------------------------------------
generate_nebula_report <- function(all_results, failed_types, comparison) {
  report_file <- file.path(
    PATHS$deg_results,
    comparison,
    paste0(comparison, "_NEBULA_report.txt")
  )
  
  # Statistics
  total_types <- length(all_results) + length(failed_types)
  successful_types <- length(all_results)
  
  # Summarize DEGs
  deg_summary <- do.call(rbind, lapply(names(all_results), function(ct) {
    data.frame(
      cell_type = ct,
      total_genes = nrow(all_results[[ct]]$all),
      sig_genes = nrow(all_results[[ct]]$significant),
      up_genes = sum(all_results[[ct]]$significant$log2FC > 0),
      down_genes = sum(all_results[[ct]]$significant$log2FC < 0)
    )
  }))
  
  report_lines <- c(
    "========== NEBULA Differential Expression Report ==========",
    paste("Comparison:", comparison),
    paste("Analysis time:", Sys.time()),
    "",
    "Analysis Statistics:",
    paste("  Total cell types:", total_types),
    paste("  Successful analyses:", successful_types),
    paste("  Failed analyses:", length(failed_types)),
    "",
    "DEG Summary:"
  )
  
  # Add per cell type results
  for (i in 1:nrow(deg_summary)) {
    report_lines <- c(
      report_lines,
      paste("  ", deg_summary$cell_type[i], ":"),
      paste("    Significant genes:", deg_summary$sig_genes[i],
            "(Up:", deg_summary$up_genes[i],
            ", Down:", deg_summary$down_genes[i], ")")
    )
  }
  
  if (length(failed_types) > 0) {
    report_lines <- c(
      report_lines,
      "",
      "Failed Cell Types:",
      paste("  ", failed_types, collapse = "\n  ")
    )
  }
  
  report_lines <- c(
    report_lines,
    "",
    "Analysis Parameters:",
    "  Method: NEBULA-HL",
    "  FDR threshold: 0.05",
    "  |log2FC| threshold: 0.2",
    "  Min expression: 10%",
    "  Covariates: Sex, age, PMI, no.nuclei.after.filter",
    "  Random effect: sampID"
  )
  
  writeLines(report_lines, report_file)
  message("\n  NEBULA report saved: ", report_file)
  
  # Save summary table
  write.csv(
    deg_summary,
    file.path(PATHS$deg_results, comparison, paste0(comparison, "_DEG_summary.csv")),
    row.names = FALSE
  )
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message(strrep("=", 70))
  message("Step 7: Differential Expression Analysis with NEBULA")
  message(strrep("=", 70))
  
  # Process each comparison
  comparisons <- c("AD.NORMAL", "PD.NORMAL", "DLB.NORMAL")
  
  for (comparison in comparisons) {
    tryCatch({
      run_nebula_batch(comparison)
    }, error = function(e) {
      warning("Failed to process ", comparison, ": ", e$message)
    })
  }
  
  message("\n", strrep("=", 70))
  message("✓ Step 7 completed successfully!")
  message(strrep("=", 70))
}

# Run if executed as script
if (!interactive()) {
  main()
}

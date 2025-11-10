#!/usr/bin/env Rscript
# ==============================================================================
# Step 6: Cell Type Annotation with scMAYOmap
# ==============================================================================
# Description: Annotate cell clusters with brain-specific cell types
# Input: Clustered Seurat objects from Step 5
# Output: Annotated Seurat objects with cell type identities
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(presto)
  library(reshape2)
  library(tidyr)
  library(data.table)
})

# Install scMAYOmap if needed (custom package)
if (!requireNamespace("scMayoMap", quietly = TRUE)) {
  message("scMAYOmap not found. Please install from:")
  message("  devtools::install_github('Mayo-Clinic/scMAYOmap')")
  message("Or use alternative marker-based annotation")
}

# ------------------------------------------------------------------------------
# Function: Prepare Brain-Specific Database
# ------------------------------------------------------------------------------
prepare_mayo_database <- function(tissue_type = "brain") {
  message("\nPreparing ", tissue_type, "-specific scMAYOmap database...")
  
  # Load scMAYOmap database
  if (requireNamespace("scMayoMap", quietly = TRUE)) {
    library(scMayoMap)
    
    # Get tissue-specific cell types
    tissue_ct <- grep(tissue_type, 
                     colnames(scMayoMapDatabase[scMayoMapDatabase$tissue == tissue_type, ]),
                     value = TRUE)
    
    # Define 22 brain cell types
    brain_cell_types <- c(
      "brain:Astrocyte",
      "brain:Excitatory neuron", 
      "brain:Inhibitory neuron",
      "brain:Microglia",
      "brain:Oligodendrocyte",
      "brain:Oligodendrocyte precursor cell",
      "brain:Endothelial cell",
      "brain:Pericyte"
    )
    
    # Filter available types
    available_types <- brain_cell_types[brain_cell_types %in% tissue_ct]
    
    if (length(available_types) == 0) {
      warning("No brain cell types found in database")
      return(NULL)
    }
    
    # Build database subset
    tissue_db <- scMayoMapDatabase[, c('tissue', 'gene', available_types)]
    tissue_db <- tissue_db[tissue_db$tissue == tissue_type, ]
    
    message("  Database contains ", length(available_types), " cell types")
    message("  Gene count: ", nrow(tissue_db))
    
    return(list(database = tissue_db, cell_types = available_types))
    
  } else {
    # Alternative: Use known marker genes
    message("  Using marker-based annotation instead")
    return(create_marker_database())
  }
}

# ------------------------------------------------------------------------------
# Function: Alternative Marker-Based Database
# ------------------------------------------------------------------------------
create_marker_database <- function() {
  message("Creating marker gene database...")
  
  marker_genes <- list(
    Astrocyte = c("SLC1A2", "GFAP", "AQP4", "ALDH1L1", "SOX9"),
    Excitatory_neuron = c("SLC17A7", "CAMK2A", "RBFOX3", "SATB2", "TBR1"),
    Inhibitory_neuron = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP"),
    Microglia = c("AIF1", "CX3CR1", "P2RY12", "TMEM119", "CSF1R"),
    Oligodendrocyte = c("MOBP", "MBP", "PLP1", "MAG", "MOG"),
    OPC = c("PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2"),
    Endothelial = c("CLDN5", "FLT1", "PECAM1", "VWF"),
    Pericyte = c("PDGFRB", "RGS5", "ACTA2", "MYH11")
  )
  
  return(list(markers = marker_genes))
}

# ------------------------------------------------------------------------------
# Function: Find Cluster Markers
# ------------------------------------------------------------------------------
find_cluster_markers <- function(seurat_obj, method = "wilcox",
                                min_pct = 0.25, logfc_threshold = 0.25) {
  message("\nFinding cluster marker genes...")
  message("  Method: ", method)
  
  # Set parallelization if available
  if (requireNamespace("future", quietly = TRUE)) {
    library(future)
    plan("multicore", workers = 4)
  }
  
  # Find markers
  seurat_markers <- tryCatch({
    FindAllMarkers(
      seurat_obj,
      method = method,
      only.pos = TRUE,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      verbose = FALSE,
      return.thresh = 0.05
    )
  }, error = function(e) {
    warning("FindAllMarkers failed, using wilcox test")
    FindAllMarkers(
      seurat_obj,
      method = "wilcox",
      only.pos = TRUE,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      verbose = FALSE
    )
  })
  
  message("  Found ", nrow(seurat_markers), " marker genes")
  
  # Select top markers per cluster
  top_markers <- seurat_markers %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) %>%
    as.data.frame()
  
  return(top_markers)
}

# ------------------------------------------------------------------------------
# Function: Annotate with scMAYOmap
# ------------------------------------------------------------------------------
annotate_with_mayomap <- function(seurat_obj, mayo_db, comparison = "AD.NORMAL") {
  message("\nRunning scMAYOmap annotation...")
  
  DefaultAssay(seurat_obj) <- "integrated"
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  
  # Find marker genes
  seurat_markers <- find_cluster_markers(seurat_obj)
  
  # Run scMAYOmap
  if (requireNamespace("scMayoMap", quietly = TRUE)) {
    library(scMayoMap)
    
    scMayoMap_obj <- tryCatch({
      scMayoMap(
        data = seurat_markers,
        database = mayo_db$database,
        tissue = 'brain'
      )
    }, error = function(e) {
      stop("scMAYOmap annotation failed: ", e$message)
    })
    
    # Save visualization
    mayo_plot_file <- file.path(
      PATHS$figures,
      "annotation",
      paste0(comparison, "_scmayomap_heatmap.png")
    )
    
    dir.create(dirname(mayo_plot_file), recursive = TRUE, showWarnings = FALSE)
    
    png(mayo_plot_file, width = 2500, height = 1500, res = 300)
    print(scMayoMap.plot(scMayoMap.object = scMayoMap_obj))
    dev.off()
    
    message("  Annotation heatmap saved: ", mayo_plot_file)
    
    # Save scMAYOmap object
    mayo_obj_file <- file.path(
      PATHS$annotation_results,
      paste0(comparison, "_scMAYOmap.rds")
    )
    saveRDS(scMayoMap_obj, mayo_obj_file)
    
    return(scMayoMap_obj)
    
  } else {
    # Alternative annotation using markers
    return(annotate_with_markers(seurat_obj, mayo_db$markers, seurat_markers))
  }
}

# ------------------------------------------------------------------------------
# Function: Alternative Marker-Based Annotation
# ------------------------------------------------------------------------------
annotate_with_markers <- function(seurat_obj, marker_database, cluster_markers) {
  message("Using marker-based annotation...")
  
  clusters <- unique(cluster_markers$cluster)
  annotations <- list()
  
  for (clust in clusters) {
    cluster_genes <- cluster_markers %>%
      filter(cluster == clust) %>%
      pull(gene)
    
    # Calculate overlap with each cell type
    overlaps <- sapply(marker_database, function(markers) {
      length(intersect(cluster_genes, markers)) / length(markers)
    })
    
    # Assign cell type with highest overlap
    best_match <- names(which.max(overlaps))
    confidence <- max(overlaps)
    
    annotations[[as.character(clust)]] <- list(
      cluster = clust,
      celltype = best_match,
      score = confidence
    )
    
    message("  Cluster ", clust, " -> ", best_match, " (score: ", 
            round(confidence, 3), ")")
  }
  
  # Create result object
  result <- list(
    markers = data.frame(
      cluster = sapply(annotations, function(x) x$cluster),
      celltype = sapply(annotations, function(x) x$celltype),
      score = sapply(annotations, function(x) x$score)
    )
  )
  
  return(result)
}

# ------------------------------------------------------------------------------
# Function: Process Mayo Annotations
# ------------------------------------------------------------------------------
process_mayo_annotations <- function(mayo_obj, seurat_obj) {
  message("\nProcessing annotation results...")
  
  # Extract cluster numbers
  clusters <- as.numeric(unique(mayo_obj$markers$cluster))
  
  # For each cluster, select best match
  cell_type_annotations <- list()
  
  for (i in clusters) {
    cluster_matches <- mayo_obj$markers[mayo_obj$markers$cluster == i, ]
    
    # Select highest confidence match (last row)
    best_match <- cluster_matches$celltype[nrow(cluster_matches)]
    confidence_score <- cluster_matches$score[nrow(cluster_matches)]
    
    cell_type_annotations[[i]] <- list(
      cluster = i,
      cell_type = best_match,
      confidence = confidence_score
    )
    
    message("  Cluster ", i, " -> ", best_match,
            " (confidence: ", round(confidence_score, 3), ")")
  }
  
  # Convert to dataframe
  annotation_df <- do.call(rbind, lapply(cell_type_annotations, as.data.frame))
  
  # Check low confidence annotations
  low_confidence <- annotation_df[annotation_df$confidence < 0.5, ]
  if (nrow(low_confidence) > 0) {
    warning("Low confidence annotations (< 0.5):")
    print(low_confidence)
  }
  
  return(annotation_df)
}

# ------------------------------------------------------------------------------
# Function: Add Numeric Suffix to Cell Types
# ------------------------------------------------------------------------------
add_cell_type_suffix <- function(cell_types) {
  # Count occurrences
  type_counts <- table(cell_types)
  
  # Add numeric suffix for duplicates
  result <- character(length(cell_types))
  type_counter <- list()
  
  for (i in seq_along(cell_types)) {
    type <- cell_types[i]
    
    if (type_counts[type] > 1) {
      # Add suffix
      if (is.null(type_counter[[type]])) {
        type_counter[[type]] <- 1
      } else {
        type_counter[[type]] <- type_counter[[type]] + 1
      }
      result[i] <- paste0(type, type_counter[[type]])
    } else {
      # Still add "1" for consistency
      result[i] <- paste0(type, "1")
    }
  }
  
  return(result)
}

# ------------------------------------------------------------------------------
# Function: Integrate Annotations into Seurat Object
# ------------------------------------------------------------------------------
integrate_annotations <- function(seurat_obj, annotation_df) {
  message("\nIntegrating annotations into Seurat object...")
  
  # Create mapping dataframe
  cell_type_df <- data.frame(
    seurat_clusters = annotation_df$cluster,
    mayo.cell.type = annotation_df$cell_type,
    mayo.confidence = annotation_df$confidence
  )
  
  # Add numeric suffix
  cell_type_df$mayo.cell.type.number <- add_cell_type_suffix(cell_type_df$mayo.cell.type)
  
  # Merge with metadata
  metadata <- seurat_obj@meta.data
  metadata$row.names <- rownames(metadata)
  
  new_metadata <- merge(
    metadata,
    cell_type_df,
    by = "seurat_clusters",
    all.x = TRUE
  )
  
  # Restore row names
  rownames(new_metadata) <- new_metadata$row.names
  new_metadata$row.names <- NULL
  
  # Check for unannotated cells
  na_cells <- sum(is.na(new_metadata$mayo.cell.type))
  if (na_cells > 0) {
    warning(na_cells, " cells without annotation")
  }
  
  # Update Seurat object
  seurat_obj@meta.data <- new_metadata
  
  message("  Annotation complete")
  
  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# Function: Generate Annotation Visualizations
# ------------------------------------------------------------------------------
visualize_annotations <- function(seurat_obj, comparison = "AD.NORMAL") {
  message("\nGenerating annotation visualizations...")
  
  output_dir <- file.path(PATHS$figures, "annotation", comparison)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. UMAP by Mayo cell type
  ct <- unique(seurat_obj$mayo.cell.type)
  ct <- ct[order(ct)]
  
  p1 <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "mayo.cell.type",
    order = ct,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.1
  ) +
    labs(title = "Cell Type Annotation (scMAYOmap)") +
    theme(legend.position = "right")
  
  ggsave(
    file.path(output_dir, "mayo_cell_types.png"),
    p1, width = 12, height = 8, dpi = 300
  )
  
  # 2. UMAP with numeric labels
  p2 <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "mayo.cell.type.number",
    label = TRUE,
    label.size = 4,
    pt.size = 0.1,
    repel = TRUE
  ) +
    labs(title = "Cell Type Clusters")
  
  ggsave(
    file.path(output_dir, "mayo_clusters_numbered.png"),
    p2, width = 12, height = 8, dpi = 300
  )
  
  # 3. Confidence score distribution
  conf_data <- seurat_obj@meta.data %>%
    group_by(mayo.cell.type) %>%
    summarise(
      mean_confidence = mean(mayo.confidence, na.rm = TRUE),
      sd_confidence = sd(mayo.confidence, na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_confidence))
  
  p3 <- ggplot(conf_data, aes(x = reorder(mayo.cell.type, mean_confidence),
                               y = mean_confidence)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = mean_confidence - sd_confidence,
                     ymax = mean_confidence + sd_confidence),
                 width = 0.2) +
    coord_flip() +
    labs(
      title = "Annotation Confidence by Cell Type",
      x = "Cell Type",
      y = "Mean Confidence Score"
    ) +
    theme_minimal()
  
  ggsave(
    file.path(output_dir, "annotation_confidence.png"),
    p3, width = 10, height = 8, dpi = 300
  )
  
  # 4. Cell type proportions by condition
  prop_data <- seurat_obj@meta.data %>%
    group_by(orig.ident, mayo.cell.type) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(orig.ident) %>%
    mutate(proportion = n / sum(n) * 100)
  
  p4 <- ggplot(prop_data, aes(x = orig.ident, y = proportion, fill = mayo.cell.type)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Cell Type Proportions by Condition",
      x = "Condition",
      y = "Proportion (%)",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(
    file.path(output_dir, "cell_type_proportions.png"),
    p4, width = 12, height = 8, dpi = 300
  )
  
  message("  Visualizations saved to: ", output_dir)
}

# ------------------------------------------------------------------------------
# Function: Generate Annotation Report
# ------------------------------------------------------------------------------
generate_annotation_report <- function(seurat_obj, annotation_df, comparison) {
  report_file <- file.path(
    PATHS$annotation_results,
    paste0(comparison, "_annotation_report.txt")
  )
  
  # Collect statistics
  cell_type_counts <- table(seurat_obj$mayo.cell.type)
  mean_confidence <- mean(annotation_df$confidence, na.rm = TRUE)
  
  report_lines <- c(
    "========== scMAYOmap Annotation Report ==========",
    paste("Comparison:", comparison),
    paste("Analysis time:", Sys.time()),
    "",
    "Annotation Statistics:",
    paste("  Annotated clusters:", nrow(annotation_df)),
    paste("  Identified cell types:", length(unique(annotation_df$cell_type))),
    paste("  Mean confidence:", round(mean_confidence, 3)),
    "",
    "Cell Type Distribution:"
  )
  
  # Add cell counts per type
  for (ct in names(sort(cell_type_counts, decreasing = TRUE))) {
    report_lines <- c(
      report_lines,
      paste("  ", ct, ":", cell_type_counts[ct], "cells")
    )
  }
  
  report_lines <- c(
    report_lines,
    "",
    "Method Parameters:",
    "  Database: scMAYOmap brain database",
    "  Cell types: 22",
    "  Marker method: MAST",
    "  Min expression: 0.25",
    "  Min log2FC: 0.25"
  )
  
  writeLines(report_lines, report_file)
  message("  Annotation report saved: ", report_file)
}

# ------------------------------------------------------------------------------
# Function: Main Annotation Pipeline
# ------------------------------------------------------------------------------
run_mayo_annotation <- function(comparison = "AD.NORMAL") {
  message("\n", strrep("=", 70))
  message("Starting Cell Type Annotation: ", comparison)
  message(strrep("=", 70))
  
  # Read UMAP object
  input_file <- file.path(
    PATHS$integration_results,
    comparison,
    paste0(comparison, ".integrated_PCA_UMAP.rds")
  )
  
  if (!file.exists(input_file)) {
    stop("Please run Step 5 (UMAP and clustering) first")
  }
  
  seurat_obj <- readRDS(input_file)
  
  # Prepare Mayo database
  mayo_db <- prepare_mayo_database("brain")
  
  if (is.null(mayo_db)) {
    stop("Failed to prepare annotation database")
  }
  
  # Run scMAYOmap annotation
  mayo_obj <- annotate_with_mayomap(seurat_obj, mayo_db, comparison)
  
  # Process annotation results
  annotation_df <- process_mayo_annotations(mayo_obj, seurat_obj)
  
  # Integrate annotations
  seurat_obj <- integrate_annotations(seurat_obj, annotation_df)
  
  # Generate visualizations
  visualize_annotations(seurat_obj, comparison)
  
  # Generate report
  generate_annotation_report(seurat_obj, annotation_df, comparison)
  
  # Save annotated object
  output_file <- file.path(
    PATHS$annotation_results,
    paste0(comparison, ".integrated_PCA_UMAP_annotated.rds")
  )
  
  saveRDS(seurat_obj, output_file)
  
  message("\n", strrep("=", 70))
  message("✓ Annotation complete: ", output_file)
  message(strrep("=", 70))
  
  return(seurat_obj)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message(strrep("=", 70))
  message("Step 6: Cell Type Annotation with scMAYOmap")
  message(strrep("=", 70))
  
  # Process each comparison
  comparisons <- c("AD.NORMAL", "PD.NORMAL", "DLB.NORMAL")
  
  for (comparison in comparisons) {
    tryCatch({
      run_mayo_annotation(comparison)
    }, error = function(e) {
      warning("Failed to process ", comparison, ": ", e$message)
    })
  }
  
  message("\n", strrep("=", 70))
  message("✓ Step 6 completed successfully!")
  message(strrep("=", 70))
}

# Run if executed as script
if (!interactive()) {
  main()
}

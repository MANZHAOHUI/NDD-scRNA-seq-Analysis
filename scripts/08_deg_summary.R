#!/usr/bin/env Rscript
# ==============================================================================
# Step 8: DEG Summarization and Pathway Analysis Preparation
# ==============================================================================
# Description: Aggregate DEG results and prepare for pathway enrichment
# Input: NEBULA results from Step 7
# Output: Aggregated gene lists, Metascape input files
# ==============================================================================

# Load configuration
source("scripts/config.R")

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

# Count results by cluster
countResults <- function(datatable, variable, new_column_name) {
  tryCatch({
    if (variable == '.N') {
      counts <- datatable[, .N, by = cluster]
      colnames(counts)[2] <- new_column_name
    } else {
      variable <- parse(text = variable)
      counts <- datatable[, unique(eval(variable)), by = cluster]
      counts <- counts[, .N, by = cluster]
      colnames(counts)[2] <- new_column_name
    }
    
    if (exists("df_degs")) {
      counts <- full_join(x = df_degs, y = counts, by = 'cluster')
    }
    
    return(counts)
  }, error = function(e) {
    warning("countResults error: ", e$message)
    return(data.frame())
  })
}

# Generate gene name lists for Metascape
getResultNames <- function(datatable, variable) {
  tryCatch({
    variable <- parse(text = variable)
    res_names <- datatable[, paste(unique(eval(variable)), collapse = ','), by = ident]
    res_names$ident <- paste0(res_names$ident, '\t')
    res_names <- as.matrix(res_names) %>% t %>% as.character()
    return(res_names)
  }, error = function(e) {
    warning("getResultNames error: ", e$message)
    return(character())
  })
}

# ------------------------------------------------------------------------------
# Main Processing Loop
# ------------------------------------------------------------------------------
comparisons <- c("AD.NORMAL", "PD.NORMAL", "DLB.NORMAL")

for (comparison in comparisons) {
  message("\n", strrep("=", 70))
  message("Processing comparison: ", comparison)
  message(strrep("=", 70))
  
  # Set paths
  in.dir <- file.path(PATHS$deg_results, comparison, "NEBULA")
  out.dir <- file.path(PATHS$pathway_analysis, comparison)
  
  # Create output directory
  if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = TRUE)
    message("Created output directory: ", out.dir)
  }
  
  setwd(in.dir)
  
  # ============================================================================
  # Read NEBULA Results
  # ============================================================================
  degs_list <- list.files(pattern = "_all_genes_NEBULA.csv$")
  
  if (length(degs_list) == 0) {
    warning("No NEBULA results found, skipping ", comparison)
    next
  }
  
  message("Found ", length(degs_list), " cell type results")
  
  # Batch read and apply thresholds
  degs.all <- lapply(degs_list, function(x) {
    tryCatch({
      tmp <- read.csv(x)
      # Apply thresholds
      tmp <- tmp[abs(tmp$log2FC) >= 0.2 & tmp$fdr < 0.05, ]
      return(tmp)
    }, error = function(e) {
      warning("Failed to read: ", x, " - ", e$message)
      return(data.frame())
    })
  })
  
  names(degs.all) <- degs_list
  degs.all <- degs.all[sapply(degs.all, nrow) > 0]
  
  if (length(degs.all) == 0) {
    warning("No significant DEGs found, skipping ", comparison)
    next
  }
  
  message("Cell types with significant DEGs: ", length(degs.all))
  
  # ============================================================================
  # Data Organization
  # ============================================================================
  degs <- lapply(1:length(degs.all), function(x) {
    tmp <- degs.all[[x]]
    
    # Extract cell type name
    file_name <- names(degs.all)[x]
    cluster <- gsub(paste0(comparison, "\\.|_all_genes_NEBULA\\.csv"), "", file_name)
    
    tmp$cluster <- cluster
    tmp$disease <- comparison
    tmp$ident <- paste0(comparison, "_", cluster)
    
    # Separate upregulated and downregulated
    tmp.up <- tmp[tmp$log2FC > 0, ] %>% arrange(fdr, desc(abs(log2FC)))
    tmp.down <- tmp[tmp$log2FC < 0, ] %>% arrange(fdr, desc(abs(log2FC)))
    
    # Limit to top 3000 genes per direction
    if (nrow(tmp.up) > 3000) {
      message("  ", cluster, ": Limiting upregulated genes to 3000")
      tmp.up <- tmp.up[1:3000, ]
    }
    
    if (nrow(tmp.down) > 3000) {
      message("  ", cluster, ": Limiting downregulated genes to 3000")
      tmp.down <- tmp.down[1:3000, ]
    }
    
    return(bind_rows(tmp.up, tmp.down))
  }) %>% bind_rows() %>% as.data.table()
  
  # ============================================================================
  # Generate Output Files
  # ============================================================================
  setwd(out.dir)
  
  # Separate upregulated and downregulated genes
  degs.up <- degs[degs$log2FC > 0, c("gene", "cluster", "disease", "ident")]
  degs.down <- degs[degs$log2FC < 0, c("gene", "cluster", "disease", "ident")]
  
  # Save CSV files
  write.csv(degs.up, paste0(comparison, ".up_regulated_genes.csv"), row.names = FALSE)
  write.csv(degs.down, paste0(comparison, ".down_regulated_genes.csv"), row.names = FALSE)
  message("Gene lists saved")
  
  # Generate Metascape input files
  if (nrow(degs.up) > 0) {
    degs.up_names <- getResultNames(degs.up, 'gene')
    write(degs.up_names, paste0(comparison, ".up_regulated.metascape.txt"))
  }
  
  if (nrow(degs.down) > 0) {
    degs.down_names <- getResultNames(degs.down, 'gene')
    write(degs.down_names, paste0(comparison, ".down_regulated.metascape.txt"))
  }
  
  # ============================================================================
  # Generate Statistics
  # ============================================================================
  
  # Read complete results for statistics
  degs_complete <- lapply(degs_list, function(file_name) {
    tryCatch({
      cluster <- gsub(paste0(comparison, "\\.|_all_genes_NEBULA\\.csv"), "", file_name)
      tmp <- read.csv(file_name)
      tmp$cluster <- cluster
      return(tmp)
    }, error = function(e) {
      return(data.frame())
    })
  }) %>% bind_rows() %>% as.data.table()
  
  # Apply filters
  degs_filtered <- degs_complete[abs(degs_complete$log2FC) >= 0.2 & 
                                  degs_complete$fdr < 0.05, ]
  
  # Separate by direction
  degs_up_stats <- degs_filtered[degs_filtered$log2FC > 0, ]
  degs_down_stats <- degs_filtered[degs_filtered$log2FC < 0, ]
  
  # Count DEGs per cluster
  df_degs <- degs_filtered[, .(unique(gene)), by = cluster][, .N, by = cluster]
  colnames(df_degs) <- c('cluster', 'DEGs_total')
  
  df_degs <- countResults(degs_up_stats, 'gene', 'DEGs_upregulated')
  df_degs <- countResults(degs_down_stats, 'gene', 'DEGs_downregulated')
  
  df_degs[is.na(df_degs)] <- 0
  
  write.csv(df_degs, paste0(comparison, ".DEG_summary_statistics.csv"), row.names = FALSE)
  
  # ============================================================================
  # Generate Visualization
  # ============================================================================
  tryCatch({
    df_plot <- melt(df_degs, 
                   id.vars = "cluster",
                   measure.vars = c("DEGs_upregulated", "DEGs_downregulated"),
                   variable.name = "direction", 
                   value.name = "count")
    
    p <- ggplot(df_plot, aes(x = reorder(cluster, count), y = count, fill = direction)) +
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() +
      scale_fill_manual(
        values = c("DEGs_upregulated" = "red", "DEGs_downregulated" = "blue"),
        labels = c("Upregulated", "Downregulated")
      ) +
      labs(
        title = paste("Differential Expression Genes -", comparison),
        x = "Cell Type",
        y = "Number of Genes",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14))
    
    ggsave(
      paste0(comparison, ".DEG_barplot.png"),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    message("Visualization generated")
    
  }, error = function(e) {
    warning("Visualization failed: ", e$message)
  })
  
  # ============================================================================
  # Generate Analysis Report
  # ============================================================================
  report_lines <- c(
    paste("Differential Expression Summary Report -", comparison),
    paste("Generated:", Sys.Date()),
    "",
    "Data Overview:",
    paste("  Analyzed cell types:", length(unique(degs_complete$cluster))),
    paste("  Cell types with significant DEGs:", length(unique(degs_filtered$cluster))),
    paste("  Total significant genes:", nrow(degs_filtered)),
    paste("  Upregulated genes:", nrow(degs_up_stats)),
    paste("  Downregulated genes:", nrow(degs_down_stats)),
    "",
    "Filtering Criteria:",
    "  FDR < 0.05",
    "  |log2FC| >= 0.2",
    "  Metascape gene limit: 3000"
  )
  
  writeLines(report_lines, paste0(comparison, ".analysis_report.txt"))
  message("Analysis report generated")
}

message("\n", strrep("=", 70))
message("✓ All analyses complete!")
message(strrep("=", 70))

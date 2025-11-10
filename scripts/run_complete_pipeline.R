#!/usr/bin/env Rscript
# Master Pipeline - Run All Steps

source("scripts/config.R")

run_complete_pipeline <- function(steps = 1:13) {
  script_names <- c(
    "01_download_data.R",
    "02_cellranger_processing.sh",
    "03_doublet_detection.R",
    "04_integration.R",
    "05_umap_clustering.R",
    "06_cell_annotation.R",
    "07_differential_expression.R",
    "08_deg_summary.R",
    "09_depletion_analysis.R",
    "10_vulnerability_analysis.R",
    "11_cellchat_analysis.R",
    "12_ndd_comparison.R",
    "13_pairwise_comparison.R"
  )
  
  for (i in steps) {
    script <- script_names[i]
    message("\n", strrep("=", 70))
    message("Running Step ", i, ": ", script)
    message(strrep("=", 70))
    
    if (grepl("\\.sh$", script)) {
      system(paste("bash scripts/", script))
    } else {
      source(paste0("scripts/", script))
    }
  }
  
  message("\n", strrep("=", 70))
  message("✓ Pipeline complete!")
  message(strrep("=", 70))
}

# Run all steps
if (!interactive()) {
  run_complete_pipeline()
}

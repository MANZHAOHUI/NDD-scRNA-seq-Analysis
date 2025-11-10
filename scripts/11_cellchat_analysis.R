#!/usr/bin/env Rscript
# Step 11: Cell-Cell Communication Analysis

source("scripts/config.R")

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(ggplot2)
})

# Create CellChat object
create_cellchat <- function(seurat_obj, group_name) {
  message("Creating CellChat for: ", group_name)
  
  DefaultAssay(seurat_obj) <- "RNA"
  data.input <- seurat_obj@assays$RNA@data
  
  # Extract metadata
  meta <- seurat_obj@meta.data
  
  # Filter by group
  if (group_name == "AD") {
    cell.use <- rownames(meta)[meta$orig.ident == "LOAD"]
  } else {
    cell.use <- rownames(meta)[meta$orig.ident == "Normal"]
  }
  
  data.input.group <- data.input[, cell.use]
  meta.group <- meta[cell.use, ]
  
  message("  Cells: ", length(cell.use))
  
  # Prepare CellChat metadata
  meta.cellchat <- data.frame(
    labels = meta.group$cell_type,
    row.names = colnames(data.input.group)
  )
  
  # Create CellChat object
  cellchat <- createCellChat(
    object = data.input.group,
    meta = meta.cellchat,
    group.by = "labels"
  )
  
  cellchat <- addMeta(cellchat, meta = meta.cellchat)
  cellchat <- setIdent(cellchat, ident.use = "labels")
  
  return(cellchat)
}

# Run CellChat analysis
run_cellchat <- function(cellchat) {
  message("Running CellChat analysis...")
  
  # Set database
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  # Compute communication
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE, trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  df.net <- subsetCommunication(cellchat)
  message("  Detected ", nrow(df.net), " interactions")
  
  # Pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(list(cellchat = cellchat, df.net = df.net))
}

# Main
main <- function() {
  message("Step 11: Cell-Cell Communication")
  
  seurat_obj <- readRDS(file.path(PATHS$annotation_results,
                                  'AD.NORMAL.integrated_PCA_UMAP_annotated.rds'))
  
  # Create CellChat objects
  cellchat_ad <- create_cellchat(seurat_obj, "AD")
  cellchat_normal <- create_cellchat(seurat_obj, "Normal")
  
  # Run analyses
  results_ad <- run_cellchat(cellchat_ad)
  results_normal <- run_cellchat(cellchat_normal)
  
  # Save results
  output_dir <- file.path(PATHS$cellchat_results, "AD_vs_Normal")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(results_ad, file.path(output_dir, 'cellchat_AD.rds'))
  saveRDS(results_normal, file.path(output_dir, 'cellchat_Normal.rds'))
  
  message("✓ Step 11 complete")
}

if (!interactive()) main()

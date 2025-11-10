#!/usr/bin/env Rscript
# Step 9: Vulnerable Neuron Depletion Analysis

source("scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(betareg)
  library(lmtest)
  library(ggplot2)
})

# Extract and recluster neurons
extract_and_recluster_neurons <- function(input_file, output_dir) {
  message("Extracting neurons...")
  
  seurat_obj <- readRDS(input_file)
  DefaultAssay(seurat_obj) <- 'integrated'
  
  # Filter neurons (Excitatory and Inhibitory)
  Idents(seurat_obj) <- 'mayo.cell.type.id'
  neuro <- subset(seurat_obj, idents = c('Exc', 'Inh'))
  
  message("  Neuron count: ", ncol(neuro))
  
  # Recluster
  neuro <- FindVariableFeatures(neuro, nfeatures = 10000)
  neuro <- RunPCA(neuro, npcs = 30)
  neuro <- RunUMAP(neuro, dims = 1:30, seed.use = 42)
  neuro <- FindNeighbors(neuro, dims = 1:30)
  neuro <- FindClusters(neuro, resolution = 1)
  
  # Save
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(neuro, file.path(output_dir, 'AD.depletion.rds'))
  
  return(neuro)
}

# Beta regression for depletion
run_beta_regression <- function(neuro_obj, output_dir) {
  message("Running beta regression...")
  
  neuro.df <- neuro_obj@meta.data
  neuro.df$orig.ident <- as.factor(neuro.df$orig.ident)
  
  # Calculate proportions
  sc_data.cluster <- neuro.df %>%
    group_by(orig.ident, sampID, mayo.neuron.cell.type.number) %>%
    summarize(subcluster.n = n(), .groups = 'drop')
  
  sc_data.sample <- neuro.df %>%
    group_by(orig.ident, sampID) %>%
    summarize(sample.n = n(), .groups = 'drop')
  
  sc_data <- left_join(sc_data.cluster, sc_data.sample) %>%
    mutate(proportion = subcluster.n / sample.n)
  
  # Add covariates
  metadata_unique <- neuro.df %>%
    select(sampID, Sex, age, PMI, no.nuclei.after.filter) %>%
    distinct()
  
  sc_data <- left_join(sc_data, metadata_unique, by = "sampID")
  sc_data <- na.omit(sc_data)
  
  # Beta regression for each cluster
  results_list <- list()
  
  for (clust in unique(sc_data$mayo.neuron.cell.type.number)) {
    message("  Analyzing: ", clust)
    
    sc_data1 <- sc_data[sc_data$mayo.neuron.cell.type.number == clust, ]
    
    if (nrow(sc_data1) < 10) {
      warning("  Insufficient samples")
      next
    }
    
    # Standardize covariates
    sc_data1$age <- scale(sc_data1$age)
    sc_data1$PMI <- scale(sc_data1$PMI)
    sc_data1$no.nuclei.after.filter <- scale(sc_data1$no.nuclei.after.filter)
    
    tryCatch({
      model <- glmmTMB(
        proportion ~ orig.ident + (1|sampID) + Sex + age + PMI + no.nuclei.after.filter,
        data = sc_data1,
        family = beta_family(link = "logit")
      )
      
      model_summary <- summary(model)
      
      if (nrow(model_summary$coefficients$cond) >= 2) {
        p_value <- model_summary$coefficients$cond[2, 4]
        coefficient <- model_summary$coefficients$cond[2, 1]
        
        results_list[[clust]] <- data.frame(
          cluster = clust,
          coefficient = coefficient,
          p_value = p_value
        )
      }
    }, error = function(e) {
      warning("  Beta regression failed: ", e$message)
    })
  }
  
  # Combine and adjust
  if (length(results_list) > 0) {
    df <- do.call(rbind, results_list)
    df$adj_p <- p.adjust(df$p_value, method = 'holm')
    df <- df[order(df$adj_p), ]
    
    write.csv(df, file.path(output_dir, 'AD.beta.reg.depleted.adj.csv'), row.names = FALSE)
    
    sig_clusters <- df[df$adj_p < 0.05, ]
    if (nrow(sig_clusters) > 0) {
      message("\nSignificant depleted clusters:")
      print(sig_clusters)
    }
    
    return(df)
  }
  
  return(NULL)
}

# Main
main <- function() {
  message("Step 9: Depletion Analysis")
  
  input_file <- file.path(PATHS$annotation_results,
                         'AD.NORMAL.integrated_PCA_UMAP.mayo.annotated.rds')
  output_dir <- file.path(PATHS$depletion_results, 'AD')
  
  neuro_obj <- extract_and_recluster_neurons(input_file, output_dir)
  beta_results <- run_beta_regression(neuro_obj, output_dir)
  
  message("✓ Step 9 complete")
}

if (!interactive()) main()

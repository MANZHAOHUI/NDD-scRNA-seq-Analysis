#!/usr/bin/env Rscript
# ==============================================================================
# Step 1: Download Raw Data from Synapse
# ==============================================================================
# Description: Download FASTQ files from Synapse database
# Input: Synapse credentials
# Output: Raw FASTQ files organized by disease group
# ==============================================================================

# Load configuration
source("scripts/config.R")

# ------------------------------------------------------------------------------
# 1.1 Install and Load Synapse Client
# ------------------------------------------------------------------------------
install_synapse_client <- function() {
  message("Installing Synapse R client...")
  
  if (!requireNamespace("synapser", quietly = TRUE)) {
    install.packages("synapser", 
                    repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
  }
  
  if (!requireNamespace("synapserutils", quietly = TRUE)) {
    install.packages("synapserutils", 
                    repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
  }
  
  library(synapser)
  library(synapserutils)
  message("✓ Synapse client loaded successfully")
}

# ------------------------------------------------------------------------------
# 1.2 Secure Login to Synapse
# ------------------------------------------------------------------------------
synapse_login_safe <- function(auth_token = NULL) {
  message("\nLogging in to Synapse...")
  
  # Try environment variable first
  if (is.null(auth_token)) {
    auth_token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
  }
  
  if (auth_token == "") {
    stop("Please set SYNAPSE_AUTH_TOKEN environment variable\n",
         "Get token from: Synapse website -> Settings -> Personal Access Tokens")
  }
  
  tryCatch({
    synLogin(authToken = auth_token)
    message("✓ Successfully logged in to Synapse")
    return(TRUE)
  }, error = function(e) {
    stop("Synapse login failed: ", e$message)
  })
}

# ------------------------------------------------------------------------------
# 1.3 Download Function with Retry Mechanism
# ------------------------------------------------------------------------------
download_synapse_data <- function(synapse_id, target_dir, max_retries = 3) {
  for (attempt in 1:max_retries) {
    message(sprintf("Downloading %s (attempt %d/%d)", 
                   synapse_id, attempt, max_retries))
    
    tryCatch({
      dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      
      result <- synapserutils::syncFromSynapse(synapse_id, path = target_dir)
      
      message(sprintf("✓ Successfully downloaded %s", synapse_id))
      return(result)
      
    }, error = function(e) {
      if (attempt < max_retries) {
        message(sprintf("✗ Download failed, waiting 30 seconds before retry..."))
        Sys.sleep(30)
      } else {
        stop(sprintf("Download failed after %d attempts: %s", 
                    max_retries, e$message))
      }
    })
  }
}

# ------------------------------------------------------------------------------
# 1.4 Organize FASTQ Files by Disease Group
# ------------------------------------------------------------------------------
organize_fastq_files <- function(source_dir, sample_metadata) {
  message("\nOrganizing FASTQ files by disease group...")
  
  # Read sample metadata
  metadata <- read.csv(sample_metadata)
  
  # Find all FASTQ files
  fastq_files <- list.files(source_dir,
                            pattern = "\\.(fastq|fq)\\.gz$",
                            recursive = TRUE,
                            full.names = TRUE)
  
  message(sprintf("Found %d FASTQ files", length(fastq_files)))
  
  organized_count <- list(AD = 0, Normal = 0, PD = 0, DLB = 0)
  
  for (file in fastq_files) {
    basename_file <- basename(file)
    
    # Determine disease group from filename or metadata
    target_dir <- NULL
    for (i in 1:nrow(metadata)) {
      sample_id <- metadata$Sample.ID[i]
      group <- metadata$Group[i]
      
      if (grepl(sample_id, basename_file)) {
        target_dir <- switch(group,
                           "AD" = PATHS$raw_data_ad,
                           "Normal" = PATHS$raw_data_normal,
                           "PD" = PATHS$raw_data_pd,
                           "DLB" = PATHS$raw_data_dlb,
                           NULL)
        
        if (!is.null(target_dir)) {
          organized_count[[group]] <- organized_count[[group]] + 1
          break
        }
      }
    }
    
    # Copy file to target directory
    if (!is.null(target_dir)) {
      dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      target_file <- file.path(target_dir, basename_file)
      
      if (!file.exists(target_file)) {
        file.copy(file, target_file)
        message(sprintf("  Moved: %s -> %s", basename_file, basename(target_dir)))
      }
    }
  }
  
  # Print summary
  message("\n=== File Organization Summary ===")
  for (group in names(organized_count)) {
    message(sprintf("  %s: %d files", group, organized_count[[group]]))
  }
  
  return(organized_count)
}

# ------------------------------------------------------------------------------
# 1.5 Validate Downloads
# ------------------------------------------------------------------------------
validate_downloads <- function() {
  message("\n=== Validating Downloads ===")
  
  groups <- c("AD", "Normal", "PD", "DLB")
  group_dirs <- list(
    AD = PATHS$raw_data_ad,
    Normal = PATHS$raw_data_normal,
    PD = PATHS$raw_data_pd,
    DLB = PATHS$raw_data_dlb
  )
  
  validation_results <- list()
  
  for (group in groups) {
    dir_path <- group_dirs[[group]]
    
    if (!dir.exists(dir_path)) {
      validation_results[[group]] <- list(status = "Directory not found", count = 0)
      next
    }
    
    fastq_files <- list.files(dir_path, pattern = "\\.(fastq|fq)\\.gz$")
    r1_files <- grep("_R1", fastq_files, value = TRUE)
    r2_files <- grep("_R2", fastq_files, value = TRUE)
    
    if (length(r1_files) != length(r2_files)) {
      validation_results[[group]] <- list(
        status = "WARNING: R1/R2 mismatch",
        r1_count = length(r1_files),
        r2_count = length(r2_files)
      )
    } else {
      validation_results[[group]] <- list(
        status = "OK",
        pairs = length(r1_files)
      )
    }
  }
  
  # Print validation results
  for (group in names(validation_results)) {
    result <- validation_results[[group]]
    message(sprintf("  %s: %s", group, result$status))
    if (!is.null(result$pairs)) {
      message(sprintf("    Paired files: %d", result$pairs))
    }
  }
  
  return(validation_results)
}

# ------------------------------------------------------------------------------
# 1.6 Generate Download Report
# ------------------------------------------------------------------------------
generate_download_report <- function(validation_results) {
  report_file <- file.path(PATHS$logs,
                          paste0("download_report_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  report_lines <- c(
    "========== FASTQ Download Report ==========",
    paste("Generated:", Sys.time()),
    "",
    "Data Sources:",
    "  AD/Normal: syn51110813",
    "  PD/DLB: syn60245188",
    "",
    "Download Results:"
  )
  
  for (group in names(validation_results)) {
    result <- validation_results[[group]]
    report_lines <- c(report_lines,
                     paste("  ", group, ":", result$status))
  }
  
  writeLines(report_lines, report_file)
  message(sprintf("\n✓ Download report saved: %s", report_file))
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main <- function() {
  message("=" %R% 50)
  message("Step 1: Downloading Raw Data from Synapse")
  message("=" %R% 50)
  
  # Install and load Synapse client
  install_synapse_client()
  
  # Login to Synapse
  synapse_login_safe()
  
  # Download AD and Normal samples
  message("\n--- Downloading AD and Normal samples ---")
  ad_normal_files <- download_synapse_data(
    synapse_id = 'syn51110813',
    target_dir = PATHS$raw_data
  )
  
  # Download PD and DLB samples
  message("\n--- Downloading PD and DLB samples ---")
  pd_dlb_files <- download_synapse_data(
    synapse_id = 'syn60245188',
    target_dir = PATHS$raw_data
  )
  
  # Organize files by disease group
  sample_metadata <- file.path(PATHS$reference, "samples_metadata.csv")
  if (file.exists(sample_metadata)) {
    organize_fastq_files(PATHS$raw_data, sample_metadata)
  } else {
    warning("Sample metadata file not found, skipping organization")
  }
  
  # Validate downloads
  validation_results <- validate_downloads()
  
  # Generate report
  generate_download_report(validation_results)
  
  message("\n" %+% "=" %R% 50)
  message("✓ Step 1 completed successfully!")
  message("=" %R% 50)
}

# Run main function
if (!interactive()) {
  main()
}

#!/bin/bash
# ==============================================================================
# Step 2: Cell Ranger Processing
# ==============================================================================
# Description: Process FASTQ files with Cell Ranger to generate gene expression matrices
# Input: Raw FASTQ files
# Output: Cell Ranger outputs (filtered matrices, BAM files, QC metrics)
# ==============================================================================

set -euo pipefail
trap 'echo "Error on line $LINENO"' ERR

# Source configuration
source config.sh

# ------------------------------------------------------------------------------
# Function: Check Cell Ranger Installation
# ------------------------------------------------------------------------------
check_cellranger() {
    echo "Checking Cell Ranger installation..."
    
    if ! command -v cellranger &> /dev/null; then
        echo "ERROR: Cell Ranger not found in PATH"
        echo "Please install Cell Ranger from:"
        echo "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
        exit 1
    fi
    
    CELLRANGER_VERSION=$(cellranger --version | head -n1)
    echo "✓ Found: $CELLRANGER_VERSION"
    
    # Check reference genome
    if [ ! -d "${REFERENCE_GENOME}" ]; then
        echo "ERROR: Reference genome not found at ${REFERENCE_GENOME}"
        echo "Please download from:"
        echo "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
        exit 1
    fi
    
    echo "✓ Reference genome found: ${REFERENCE_GENOME}"
}

# ------------------------------------------------------------------------------
# Function: Run Cell Ranger Count
# ------------------------------------------------------------------------------
run_cellranger_count() {
    local SAMPLE=$1
    local FASTQ_PATH=$2
    local OUTPUT_DIR=$3
    
    echo ""
    echo "=========================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================="
    
    # Check if already processed
    if [ -f "${OUTPUT_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5" ]; then
        echo "✓ Sample ${SAMPLE} already processed, skipping..."
        return 0
    fi
    
    # Check input files exist
    if [ ! -d "${FASTQ_PATH}" ]; then
        echo "✗ ERROR: FASTQ path does not exist - ${FASTQ_PATH}"
        return 1
    fi
    
    # Create log directory
    mkdir -p "${LOG_DIR}"
    
    # Run Cell Ranger
    echo "Running Cell Ranger count..."
    cellranger count \
        --id="${SAMPLE}" \
        --transcriptome="${REFERENCE_GENOME}" \
        --fastqs="${FASTQ_PATH}" \
        --sample="${SAMPLE}" \
        --expect-cells=10000 \
        --localcores="${CORES:-16}" \
        --localmem="${MEMORY:-64}" \
        2>&1 | tee "${LOG_DIR}/cellranger_${SAMPLE}.log"
    
    # Check success
    if [ -f "${SAMPLE}/outs/web_summary.html" ]; then
        echo "✓ Successfully processed: ${SAMPLE}"
        
        # Move output to organized directory
        mv "${SAMPLE}" "${OUTPUT_DIR}/"
        
        # Generate summary
        generate_cellranger_summary "${OUTPUT_DIR}/${SAMPLE}"
        
        return 0
    else
        echo "✗ ERROR: Processing failed for ${SAMPLE}"
        return 1
    fi
}

# ------------------------------------------------------------------------------
# Function: Generate Cell Ranger Summary
# ------------------------------------------------------------------------------
generate_cellranger_summary() {
    local SAMPLE_DIR=$1
    local METRICS="${SAMPLE_DIR}/outs/metrics_summary.csv"
    
    if [ -f "${METRICS}" ]; then
        echo ""
        echo "=== QC Metrics Summary ==="
        
        # Extract key metrics
        ESTIMATED_CELLS=$(grep "Estimated Number of Cells" "${METRICS}" | cut -d',' -f2)
        MEAN_READS=$(grep "Mean Reads per Cell" "${METRICS}" | cut -d',' -f2)
        MEDIAN_GENES=$(grep "Median Genes per Cell" "${METRICS}" | cut -d',' -f2)
        SEQUENCING_SAT=$(grep "Sequencing Saturation" "${METRICS}" | cut -d',' -f2)
        
        echo "  Estimated Cells: ${ESTIMATED_CELLS}"
        echo "  Mean Reads/Cell: ${MEAN_READS}"
        echo "  Median Genes/Cell: ${MEDIAN_GENES}"
        echo "  Sequencing Saturation: ${SEQUENCING_SAT}"
        echo ""
    fi
}

# ------------------------------------------------------------------------------
# Function: Process All Samples
# ------------------------------------------------------------------------------
process_all_samples() {
    echo ""
    echo "================================================"
    echo "Starting Cell Ranger Batch Processing"
    echo "================================================"
    
    # Initialize counters
    local TOTAL=0
    local SUCCESS=0
    local FAILED=0
    
    # Process each disease group
    for GROUP_DIR in "${RAW_DATA_DIR}"/*; do
        if [ ! -d "${GROUP_DIR}" ]; then
            continue
        fi
        
        GROUP_NAME=$(basename "${GROUP_DIR}")
        echo ""
        echo "--- Processing ${GROUP_NAME} samples ---"
        
        # Process each sample in group
        for SAMPLE_DIR in "${GROUP_DIR}"/*; do
            if [ ! -d "${SAMPLE_DIR}" ]; then
                continue
            fi
            
            SAMPLE=$(basename "${SAMPLE_DIR}")
            OUTPUT_DIR="${CELLRANGER_OUTPUT}/${GROUP_NAME}"
            
            mkdir -p "${OUTPUT_DIR}"
            
            ((TOTAL++))
            
            if run_cellranger_count "${SAMPLE}" "${SAMPLE_DIR}" "${OUTPUT_DIR}"; then
                ((SUCCESS++))
            else
                ((FAILED++))
                echo "${SAMPLE}" >> "${LOG_DIR}/failed_samples.txt"
            fi
        done
    done
    
    # Print final summary
    echo ""
    echo "================================================"
    echo "Cell Ranger Processing Complete"
    echo "================================================"
    echo "Total samples: ${TOTAL}"
    echo "Successful: ${SUCCESS}"
    echo "Failed: ${FAILED}"
    
    if [ ${FAILED} -gt 0 ]; then
        echo ""
        echo "Failed samples listed in: ${LOG_DIR}/failed_samples.txt"
    fi
    
    echo "================================================"
}

# ------------------------------------------------------------------------------
# Function: Create Seurat Objects
# ------------------------------------------------------------------------------
create_seurat_objects() {
    echo ""
    echo "Creating Seurat objects from Cell Ranger outputs..."
    
    Rscript - <<'EOF'
source("scripts/config.R")
library(Seurat)
library(Matrix)

# Function to create Seurat object
create_seurat_from_cellranger <- function(sample_dir, sample_name) {
    message(sprintf("Processing: %s", sample_name))
    
    matrix_path <- file.path(sample_dir, "outs/filtered_feature_bc_matrix")
    
    if (!dir.exists(matrix_path)) {
        warning(sprintf("Matrix not found: %s", matrix_path))
        return(NULL)
    }
    
    tryCatch({
        # Read 10X data
        expression_matrix <- Read10X(matrix_path)
        
        # Create Seurat object
        seurat_obj <- CreateSeuratObject(
            counts = expression_matrix,
            project = sample_name,
            min.cells = 3,
            min.features = 200
        )
        
        # Add metadata
        seurat_obj$orig.ident <- sample_name
        seurat_obj$sampID <- sample_name
        
        # Calculate QC metrics
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
        seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
        
        # Basic filtering
        cells_before <- ncol(seurat_obj)
        seurat_obj <- subset(
            seurat_obj,
            subset = nFeature_RNA > 200 & 
                     nFeature_RNA < 8000 & 
                     percent.mt < 20
        )
        cells_after <- ncol(seurat_obj)
        
        message(sprintf("  Cells: %d -> %d (%.1f%% retained)",
                       cells_before, cells_after,
                       100 * cells_after / cells_before))
        
        # Save
        output_file <- file.path(PATHS$seurat_objects, paste0(sample_name, ".rds"))
        saveRDS(seurat_obj, output_file)
        
        return(seurat_obj)
    }, error = function(e) {
        warning(sprintf("Error processing %s: %s", sample_name, e$message))
        return(NULL)
    })
}

# Process all samples
cellranger_dirs <- list.dirs(PATHS$cellranger_output, recursive = FALSE)

for (group_dir in cellranger_dirs) {
    group_name <- basename(group_dir)
    message(sprintf("\n=== Processing %s samples ===", group_name))
    
    sample_dirs <- list.dirs(group_dir, recursive = FALSE)
    
    for (sample_dir in sample_dirs) {
        sample_name <- basename(sample_dir)
        create_seurat_from_cellranger(sample_dir, sample_name)
    }
}

message("\n✓ Seurat object creation complete")
EOF
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
main() {
    echo ""
    echo "================================================"
    echo "Step 2: Cell Ranger Processing Pipeline"
    echo "================================================"
    
    # Check prerequisites
    check_cellranger
    
    # Process all samples
    process_all_samples
    
    # Create Seurat objects
    create_seurat_objects
    
    echo ""
    echo "✓ Step 2 completed successfully!"
    echo "================================================"
}

# Run main function
main "$@"

#!/bin/bash
# Setup script for NDD single-cell analysis

echo "=========================================="
echo "NDD scRNA-seq Analysis Setup"
echo "=========================================="

# Check R installation
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not installed"
    echo "Please install R 4.2.0 or higher"
    exit 1
fi

echo "✓ R found: $(R --version | head -n1)"

# Check required system tools
for cmd in git curl wget; do
    if ! command -v $cmd &> /dev/null; then
        echo "WARNING: $cmd not found"
    else
        echo "✓ $cmd found"
    fi
done

# Create directory structure
echo ""
echo "Creating directory structure..."
mkdir -p data/{raw_data,reference,gwas,example}
mkdir -p results/{cellranger_output,seurat_objects,qc_results}
mkdir -p figures logs

# Install R packages
echo ""
echo "Installing R packages..."
Rscript -e '
packages <- c(
  "Seurat", "sctransform", "glmGamPoi",
  "ggplot2", "dplyr", "tidyr", "Matrix",
  "data.table", "reshape2", "cowplot"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    install.packages(pkg, repos = "http://cran.r-project.org")
  } else {
    message("✓ ", pkg, " already installed")
  }
}

message("\n========================================")
message("✓ Setup complete!")
message("========================================")
'

echo ""
echo "Next steps:"
echo "1. Set Synapse token: export SYNAPSE_AUTH_TOKEN='your_token'"
echo "2. Generate test data: Rscript scripts/generate_example_data.R"
echo "3. Start analysis: Rscript scripts/run_complete_pipeline.R"

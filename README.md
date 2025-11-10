# NDD Single-Cell Transcriptomics Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-4.2.0+-blue.svg)](https://www.r-project.org/)
[![GitHub Stars](https://img.shields.io/github/stars/MANZHAOHUI/NDD-scRNA-seq-Analysis?style=social)](https://github.com/MANZHAOHUI/NDD-scRNA-seq-Analysis)

> Comprehensive single-cell RNA-seq analysis pipeline comparing Alzheimer's Disease (AD), Parkinson's Disease (PD), and Dementia with Lewy Bodies (DLB)

## 🚀 Quick Start
```bash
# Clone repository
git clone https://github.com/MANZHAOHUI/NDD-scRNA-seq-Analysis.git
cd NDD-scRNA-seq-Analysis

# Run setup
bash setup.sh

# Set Synapse credentials
export SYNAPSE_AUTH_TOKEN="your_token_here"

# Test with example data
Rscript scripts/generate_example_data.R

# Run analysis
Rscript scripts/run_complete_pipeline.R
```

## 📊 Pipeline Overview

13-step automated analysis pipeline:

1. **Data Download** - Retrieve FASTQ files from Synapse
2. **Cell Ranger** - Process raw sequencing data
3. **Quality Control** - Doublet detection and removal
4. **Integration** - Batch correction with SCTransform
5. **UMAP & Clustering** - Dimensionality reduction
6. **Cell Annotation** - Identify cell types with scMAYOmap
7. **Differential Expression** - NEBULA mixed models
8. **DEG Summarization** - Prepare pathway analysis
9. **Depletion Analysis** - Identify vulnerable neurons
10. **Vulnerability Scoring** - AUCell risk assessment
11. **Cell Communication** - CellChat analysis
12. **Cross-Disease Comparison** - Compare AD, PD, DLB
13. **Pairwise Analysis** - Direct disease comparisons

## 📁 Repository Structure
```
NDD-scRNA-seq-Analysis/
├── scripts/              # All analysis scripts (Steps 1-13)
├── data/                 # Raw and processed data
├── results/              # Analysis outputs
├── figures/              # Generated visualizations
├── docs/                 # Documentation
└── tests/                # Unit tests
```

## 🔧 Requirements

### System
- Linux/macOS (Windows via WSL2)
- 128 GB RAM (minimum 64 GB)
- 32 CPU cores (minimum 16)
- 1 TB free disk space

### Software
- R 4.2.0+
- Cell Ranger 7.0.0
- Python 3.9+ (optional)

### R Packages
```r
# Core packages
install.packages(c(
  "Seurat", "sctransform", "glmGamPoi",
  "DoubletFinder", "nebula", "CellChat",
  "ggplot2", "dplyr", "Matrix"
))
```

## 📖 Documentation

- [Installation Guide](docs/installation.md)
- [Complete Tutorial](docs/tutorial.md)
- [API Reference](docs/api_reference.md)
- [Troubleshooting](docs/troubleshooting.md)

## 📝 Citation
```bibtex
@article{man2025ndd,
  title={Comparative mapping of single-cell transcriptomic landscapes in neurodegenerative diseases},
  author={Man, Zhaohui and colleagues},
  journal={Alzheimer's \& Dementia},
  year={2025}
}
```

## 📧 Contact

**Zhaohui Man**
- GitHub: [@MANZHAOHUI](https://github.com/MANZHAOHUI)
- Email: [your.email@institution.edu]

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Duke University School of Medicine
- Mayo Clinic (scMAYOmap database)
- 10X Genomics
- Synapse platform

---

**⭐ If you find this useful, please star the repository!**

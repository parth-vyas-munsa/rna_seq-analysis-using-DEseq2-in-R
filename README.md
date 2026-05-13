# rna_seq-analysis-using-DEseq2-in-R
RNA-seq differential gene expression analysis workflow using DESeq2 in R with normalization, DEG identification, PCA, volcano plots, and heatmap visualization.

# RNA-seq Differential Gene Expression Analysis using DESeq2

![R](https://img.shields.io/badge/R-Programming-blue)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green)
![RNAseq](https://img.shields.io/badge/Analysis-RNAseq-orange)
![License](https://img.shields.io/badge/License-MIT-yellow)

## Project Overview

This project presents a complete RNA-seq Differential Gene Expression (DGE) analysis workflow using DESeq2 in R. The pipeline processes raw gene count data, performs normalization, identifies significantly differentially expressed genes (DEGs), and generates publication-quality visualizations for transcriptomic analysis.

The workflow follows standard bulk RNA-seq analysis practices commonly used in bioinformatics and computational biology research.

---

## Objectives

- Perform RNA-seq differential expression analysis
- Normalize raw count data
- Identify upregulated and downregulated genes
- Generate statistical and biological insights
- Visualize transcriptomic patterns using multiple plots
- Build a reproducible bioinformatics workflow

---

## Tools & Technologies

| Tool | Purpose |
|------|----------|
| R | Statistical Computing |
| DESeq2 | Differential Expression Analysis |
| ggplot2 | Data Visualization |
| pheatmap | Heatmap Generation |
| EnhancedVolcano | Volcano Plot Visualization |
| Bioconductor | Bioinformatics Package Management |

---

## Workflow

```text
Raw Count Matrix
        ↓
Metadata Preparation
        ↓
Data Import in R
        ↓
Quality Checking
        ↓
Low Count Filtering
        ↓
DESeq2 Normalization
        ↓
Differential Expression Analysis
        ↓
Statistical Testing
        ↓
Visualization & Interpretation

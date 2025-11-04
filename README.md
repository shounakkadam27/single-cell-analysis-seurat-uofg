# Single-Cell RNA-seq Analysis (Seurat)
End-to-end single-cell RNA-seq and spatial transcriptomics analysis portfolio using Seurat (R). Includes QC, normalization, clustering, marker discovery, integration, doublet removal, CITE-sequencing and spatial workflows(Visium, Xenium, CosMx)

This repository contains hands-on **single-cell RNA-seq analysis workflows** in R using **Seurat**, demonstrating end-to-end analysis from raw data to biological insights.

### Goals
- Show practical understanding of the scRNA-seq pipeline  
- Build confidence in Seurat workflows used in research & industry  
- Practice reproducible analysis and Git version control


## 📂 Workflows

| Stage | Description | Link |
|------|------------|------|
| 1️⃣ Seurat Object + QC | Create Seurat object, calculate QC metrics, filter low-quality cells | [View](./01-seurat-object-scrna-qc/The_Seurat_object_-_scRNA_QC.md) |
| 2️⃣ Normalization + PCA + Clustering | Normalize counts (Log & SCT), find variable genes, dimensionality reduction (PCA), visualize cell structure (UMAP), and cluster cells | [View](./02-normalization-dim-reduction/Normalization_-_Dim_Reduction.md) |
| 3️⃣ Marker genes + Annotation | *(coming soon)* | – |

> Each folder contains `.Rmd`, rendered `.md`, & data used for that section.


## ✅ Skills Demonstrated

- Loading & preprocessing scRNA-seq data  
- QC metrics (genes, UMIs, mitochondrial content)
- Filtering & creating clean Seurat object
- Data handling and interpretation in R
- Version-controlled scientific workflows  
- Markdown-based reporting


## 💻 Requirements

- R (≥ 4.3)
- RStudio (recommended)
- Install core packages once:

```r
install.packages(c("Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork"))

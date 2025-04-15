# Cell-Cell Communication Analysis using Spatial Transcriptomics (GSE283269)

## Project Overview

This repository contains the complete source code, analysis pipeline, and output files for **Task 3** of the spatial transcriptomics project. The objective was to explore **cell-cell communication and interaction networks** using the **CellChat** package applied to the publicly available spatial gene expression dataset **GSE283269** from NCBI GEO.

> Candidate: **Tao Sun**  
> Task: **Cell-Cell Communication/Interaction Analysis**  
> Dataset: [GSE283269 – Spatial Transcriptomics](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283269)  
> Submission: 1 of 2 selected tasks  
> Software: R 4.3+, CellChat v1.6.1, Seurat v5

---

## Directory Structure

```bash
├── sample1a1_cellchat_outputs.tar.gz       # Output archive containing figures and tables
├── task3_cellchat.R                        # Full R script for this analysis
├── sample1a1_cellchat_all_communications.csv
├── sample1a1_cellchat_cluster_roles.csv
├── sample1a1_cellchat_[Pathway]_interactions.csv  # Pathway-specific tables (e.g., COLLAGEN)
├── sample1a1_umap_clusters.png             # Cluster visualization (UMAP)
├── sample1a1_cellchat_circle.png           # Global communication network
├── sample1a1_cellchat_dotplot_top10_pathways.png  # Top pathway interactions
├── sample1a1_cellchat_role_scatter.png     # Cluster role visualization
└── README.md                               # This file

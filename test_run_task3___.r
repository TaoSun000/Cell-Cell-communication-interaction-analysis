# =================================================================
# Task 3: Cell-Cell Communication Analysis Using Spatial Transcriptomics
# Author: Tao Sun | Date: April 14, 2025
# Dataset: GSE283269 (GSM8658911 - Sample1a1)
# =================================================================

# -------------------------------------------------------------
# Step 1: Install and Load Required Packages
# -------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install(c(
  "Seurat", "patchwork", "SingleCellExperiment", "SummarizedExperiment",
  "BiocNeighbors", "Biobase", "BiocGenerics", "GEOquery", "ComplexHeatmap"
))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("sqjin/CellChat", dependencies = TRUE)

library(Seurat)
library(CellChat)
library(patchwork)
library(GEOquery)
library(ComplexHeatmap)
library(R.utils)

# -------------------------------------------------------------
# Step 2: Download and Preprocess GSE283269 Data
# -------------------------------------------------------------

project_dir <- "E:/UBC_wang_qn3"
raw_tar <- file.path(project_dir, "GSE283269_RAW.tar")
untar_dir <- file.path(project_dir, "GSE283269")

# Download .tar file
options(timeout = 600)
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283269/suppl/GSE283269_RAW.tar",
              destfile = raw_tar, mode = "wb")
untar(raw_tar, exdir = untar_dir)

# Decompress .gz files
gz_files <- list.files(untar_dir, pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, gunzip, overwrite = TRUE)

# -------------------------------------------------------------
# Step 3: Load One Sample into Seurat and Perform Clustering
# -------------------------------------------------------------

sample1a1_files <- list(
  matrix = file.path(untar_dir, "GSM8658911_sample1a1_matrix.mtx"),
  features = file.path(untar_dir, "GSM8658911_sample1a1_features.tsv"),
  barcodes = file.path(untar_dir, "GSM8658911_sample1a1_barcodes.tsv")
)

expression_matrix <- ReadMtx(
  mtx = sample1a1_files$matrix,
  features = sample1a1_files$features,
  cells = sample1a1_files$barcodes,
  feature.column = 2
)

seurat_obj <- CreateSeuratObject(expression_matrix, project = "GSE283269_Sample1a1")

# Preprocessing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Save UMAP plot
png("E:/UBC_wang_qn3/sample1a1_umap_clusters.png", width = 1600, height = 1200, res = 200)
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Seurat UMAP Clusters")
dev.off()

# -------------------------------------------------------------
# Step 4: Create and Run CellChat Analysis
# -------------------------------------------------------------

seurat_obj$seurat_clusters_cellchat <- paste0("C", seurat_obj$seurat_clusters)
meta <- seurat_obj@meta.data
data.input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "seurat_clusters_cellchat")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "seurat_clusters_cellchat")
cellchat@DB <- CellChatDB.human

# Communication inference
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save overall communication network
png("E:/UBC_wang_qn3/sample1a1_cellchat_circle.png", width = 1600, height = 1200, res = 200)
netVisual_circle(cellchat@net$count,
                 vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE)
dev.off()

# -------------------------------------------------------------
# Step 5: Explore Pathways and Roles
# -------------------------------------------------------------

#  cellchat diagram for the CXCL network
png("E:/UBC_wang_qn3/sample1a1_cellchat_CXCL_pathway.png", width = 1600, height = 1200, res = 200)
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "circle")
dev.off()

# Save a chord diagram for the CXCL pathway
png("E:/UBC_wang_qn3/sample1a1_cellchat_CXCL_chord.png", width = 1600, height = 1200, res = 200)
netVisual_chord_cell(cellchat, signaling = "CXCL")
dev.off()

# Role analysis
cellchat <- netAnalysis_computeCentrality(cellchat)

# Heatmap (top 10 pathway)
pathway_counts <- sapply(cellchat@netP$pathways, function(path) {
  df <- subsetCommunication(cellchat, signaling = path)
  nrow(df)
})

pathway_counts <- sort(pathway_counts, decreasing = TRUE)
top10_pathways <- names(head(pathway_counts, 10))

available_probs <- names(cellchat@netP$prob)
valid_pathways <- intersect(top10_pathways, available_probs)

png("E:/UBC_wang_qn3/sample1a1_cellchat_pathway_heatmap.png", width = 1600, height = 1200, res = 200)
netVisual_heatmap(cellchat, signaling = valid_pathways)
dev.off()


# Scatter
png("E:/UBC_wang_qn3/sample1a1_cellchat_role_scatter.png", width = 1600, height = 1200, res = 200)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# -------------------------------------------------------------
# Step 6: Export Interactions and Summary
# -------------------------------------------------------------

df.net <- subsetCommunication(cellchat)
df.cxcl <- subsetCommunication(cellchat, signaling = "CXCL")

write.csv(df.net, "E:/UBC_wang_qn3/sample1a1_cellchat_all_interactions.csv", row.names = FALSE)
write.csv(df.cxcl, "E:/UBC_wang_qn3/sample1a1_cellchat_CXCL_interactions.csv", row.names = FALSE)

# Top pathways
pathway_counts <- sapply(cellchat@netP$pathways, function(path) {
  nrow(subsetCommunication(cellchat, signaling = path))
})
pathway_counts <- sort(pathway_counts, decreasing = TRUE)
top10_pathways <- names(head(pathway_counts, 10))

# Dot plot for top 10
png("E:/UBC_wang_qn3/sample1a1_cellchat_dotplot_top10_pathways.png", width = 1800, height = 1200, res = 200)
netVisual_bubble(cellchat,
                 signaling = top10_pathways,
                 sources.use = unique(cellchat@idents),
                 targets.use = unique(cellchat@idents),
                 remove.isolate = TRUE,
                 angle.x = 45)
dev.off()

# Save per-pathway tables
for (path in top10_pathways) {
  df <- subsetCommunication(cellchat, signaling = path)
  write.csv(df, file = paste0("E:/UBC_wang_qn3/sample1a1_cellchat_", path, "_interactions.csv"), row.names = FALSE)
}

# Save cluster role table
role_df <- netAnalysis_signalingRole_scatter(cellchat)$data
write.csv(role_df, "E:/UBC_wang_qn3/sample1a1_cellchat_cluster_roles.csv", row.names = FALSE)

message("Task 3 completed: Cell-Cell Communication results exported successfully.")

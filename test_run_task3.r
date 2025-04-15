# Cell-Cell Communication Analysis
# Run Test for Dr. Wang's Lab | Task 3 | Date: April 14
# -------------------------------------------------------------
#Step 1.Setup and Package Installation
# -------------------------------------------------------------
# Step 1.1: Set compatible Bioconductor version (for R 4.3.3)
# -------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.17")

# -------------------------------------------------------------
# Step 1.2: Install CRAN and Bioconductor packages
# -------------------------------------------------------------
BiocManager::install(c(
  "Seurat",
  "patchwork",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "BiocNeighbors",
  "Biobase",
  "BiocGenerics",
  "GEOquery",
  "ComplexHeatmap"
))

# -------------------------------------------------------------
# Step 1.3: Install CellChat from GitHub
# -------------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("sqjin/CellChat", dependencies = TRUE)

# -------------------------------------------------------------
# Step 1.4: Load all required libraries
# -------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(GEOquery)
library(ComplexHeatmap)

# -------------------------------------------------------------
# Step 1.5: Check for Rtools (Windows only)
# -------------------------------------------------------------
Sys.which("make")  # Should return path to Rtools 'make' if installed

# -------------------------------------------------------------
#Step 2: Download and Preprocess GSE283269 Data
# -------------------------------------------------------------
# Step 2.1 Download
# -------------------------------------------------------------
project_dir <- "E:/UBC_wang_qn3"
download_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283269/suppl/GSE283269_RAW.tar"
dest_file <- file.path(project_dir, "GSE283269_RAW.tar")

# Set longer timeout (10 minutes = 600 seconds)
options(timeout = 600)

# Redownload file (try again now)
download.file(url = download_url, destfile = dest_file, mode = "wb")

# Extract files
untar(dest_file, exdir = file.path(project_dir, "GSE283269"))

library(R.utils)
gz_files <- list.files(file.path(project_dir, "GSE283269"), pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, function(f) gunzip(f, overwrite = TRUE))

# -------------------------------------------------------------
# Step 2.2: Load One Sample into a Seurat Object to check if it works.
# -------------------------------------------------------------
library(Seurat)

# Define path to sample 1a1 (replace with other sample names for future runs)
sample_dir <- "E:/UBC_wang_qn3/GSE283269"

# Create a list object pointing to files
sample1a1_files <- list(
  "matrix" = file.path(sample_dir, "GSM8658911_sample1a1_matrix.mtx"),
  "features" = file.path(sample_dir, "GSM8658911_sample1a1_features.tsv"),
  "barcodes" = file.path(sample_dir, "GSM8658911_sample1a1_barcodes.tsv")
)

# Load expression matrix
expression_matrix <- ReadMtx(
  mtx = sample1a1_files$matrix,
  features = sample1a1_files$features,
  cells = sample1a1_files$barcodes,
  feature.column = 2
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expression_matrix, project = "GSE283269_Sample1a1")

# Check object structure
seurat_obj

# -------------------------------------------------------------
#Step 3: Normalize, Cluster, and Annotate the Data
# -------------------------------------------------------------
#Step 3.1: Standard Preprocessing in Seurat
# -------------------------------------------------------------
# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Find highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Find neighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Cluster the spots
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Define output path
output_file <- "E:/UBC_wang_qn3/sample1a1_umap_clusters.png"

# Save the UMAP cluster plot
png(filename = output_file, width = 1600, height = 1200, res = 200)
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP - Seurat Clusters (Sample1a1)")
dev.off()

# -------------------------------------------------------------
#Step 4: Run CellChat on the Seurat Object
# -------------------------------------------------------------
#Step 4.1: Create and Preprocess CellChat Object
# -------------------------------------------------------------
library(CellChat)
library(Seurat)

# Prefix seurat cluster names to avoid CellChat error with '0'
seurat_obj$seurat_clusters_cellchat <- paste0("C", seurat_obj$seurat_clusters)

# Re-extract metadata
meta.data <- seurat_obj@meta.data

# Use normalized log-expression (Seurat v5 layer)
data.input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Create CellChat object using modified cluster labels
cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "seurat_clusters_cellchat")

# Add metadata and set identity
cellchat <- addMeta(cellchat, meta = meta.data)
cellchat <- setIdent(cellchat, ident.use = "seurat_clusters_cellchat")

# Set the database (human)
cellchat@DB <- CellChatDB.human

# Continue with sub-setting and interaction inference
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# -------------------------------------------------------------
# Step 4.2: Inference of Communication Probability
# -------------------------------------------------------------
# Compute communication probability between clusters
cellchat <- computeCommunProb(cellchat)

# Filter out low-confidence interactions
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute probability at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network to get number and strength of interactions
cellchat <- aggregateNet(cellchat)

# -------------------------------------------------------------
# Step 4.3: Visualize the global communication network
# -------------------------------------------------------------
# Output file path
output_path <- "E:/UBC_wang_qn3/sample1a1_cellchat_circle.png"

# Save circular communication plot to file
png(filename = output_path, width = 1600, height = 1200, res = 200)

# Plot number of interactions between clusters
netVisual_circle(cellchat@net$count,
                 vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE,
                 label.edge = FALSE)

dev.off()

# -------------------------------------------------------------
#Step 5: Explore Specific Signaling Pathways and Communication Roles
# -------------------------------------------------------------
cellchat@netP$pathways
# -------------------------------------------------------------
#Step 5.1: Visualize a Specific Signaling Pathway (e.g., CXCL)
# -------------------------------------------------------------
# Plot cell-cell communication network for one signaling pathway
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "circle")

# Save to file (optional)
png("E:/UBC_wang_qn3/sample1a1_cellchat_CXCL_pathway.png", width = 1600, height = 1200, res = 200)
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "circle")
dev.off()

# -------------------------------------------------------------
#Step 5.2: Heatmap of Incoming and Outgoing Signals
# -------------------------------------------------------------
# Step 5.2.1: Compute centrality scores
# -------------------------------------------------------------
cellchat <- netAnalysis_computeCentrality(cellchat)

# -------------------------------------------------------------
# tep 5.2.2: Heatmap of incoming/outgoing signaling roles
# -------------------------------------------------------------
netAnalysis_signalingRole_heatmap(cellchat)

# Define output file path
output_path <- "E:/UBC_wang_qn3/sample1a1_cellchat_signalingRole_heatmap.png"

# Save heatmap to file
png(filename = output_path, width = 1600, height = 1200, res = 200)
netAnalysis_signalingRole_heatmap(cellchat)
dev.off()

# -------------------------------------------------------------
#Step 5.3: Compute and Visualize Functional Roles
# -------------------------------------------------------------
# Compute role scores for each cluster
cellchat <- netAnalysis_computeCentrality(cellchat)

# Visualize as scatter plot
netAnalysis_signalingRole_scatter(cellchat)

# Define output path
output_path <- "E:/UBC_wang_qn3/sample1a1_cellchat_role_scatter.png"

# Save the scatter plot to file
png(filename = output_path, width = 1600, height = 1200, res = 200)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# -------------------------------------------------------------
#Step 6: Extract and Interpret Cell-Cell Communication Results
# -------------------------------------------------------------
#Step 6.1: Extract All Significant Interactions
# -------------------------------------------------------------
# Extract the communication network data frame
df.net <- subsetCommunication(cellchat)

# View the top few interactions
head(df.net)

# -------------------------------------------------------------
#Step 6.2: Extract for a Specific Signaling Pathway (e.g., CXCL)
# -------------------------------------------------------------
# View only CXCL-related interactions
df.cxcl <- subsetCommunication(cellchat, signaling = "CXCL")

# Show top rows
head(df.cxcl)

# -------------------------------------------------------------
#Step 6.3: Save Interaction Table to CSV
# -------------------------------------------------------------
# Save full interaction table to CSV
write.csv(df.net, file = "E:/UBC_wang_qn3/sample1a1_cellchat_all_interactions.csv", row.names = FALSE)

# Save specific pathway interactions (e.g., CXCL)
write.csv(df.cxcl, file = "E:/UBC_wang_qn3/sample1a1_cellchat_CXCL_interactions.csv", row.names = FALSE)

# -------------------------------------------------------------
#Step 6.4: View Top Pathways or Ligand-Receptor Pairs
# -------------------------------------------------------------
#Step 6.4.1: Rank signaling pathways by overall strength
# -------------------------------------------------------------
# Visualize signaling network: which clusters send/receive most signals
netAnalysis_signalingRole_network(cellchat)
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "circle")

# Define PDF output file
pdf("E:/UBC_wang_qn3/sample1a1_cellchat_role_network_all_pathways.pdf", width = 8, height = 6)

# Loop through each signaling pathway and plot it
for (pathway in cellchat@netP$pathways) {
  netAnalysis_signalingRole_network(cellchat, signaling = pathway)
}

dev.off()

# -------------------------------------------------------------
#Step 6.4.2: Dot plot: top 2 signaling pathways and review the signaling pathways
# -------------------------------------------------------------
packageVersion("CellChat")

available_pathways <- cellchat@netP$pathways
length(available_pathways)  # Should return ~30+

top2_pathways <- names(sort(sapply(available_pathways, function(path) {
  nrow(subsetCommunication(cellchat, signaling = path))
}), decreasing = TRUE))[1:2]

# Generate dot plot for top 2 pathways only
netVisual_bubble(cellchat, 
                 signaling = top2_pathways,
                 sources.use = unique(cellchat@idents),
                 targets.use = unique(cellchat@idents),
                 remove.isolate = TRUE,
                 angle.x = 45)

png("E:/UBC_wang_qn3/sample1a1_cellchat_dotplot_top2_pathways.png", width = 900, height = 1200, res = 200)
netVisual_bubble(cellchat, 
                 signaling = top2_pathways,
                 sources.use = unique(cellchat@idents),
                 targets.use = unique(cellchat@idents),
                 remove.isolate = TRUE,
                 angle.x = 45)
dev.off()

#review the signaling pathways
head(cellchat@netP$pathways)

# Count number of LR pairs per pathway
pathway_counts <- sapply(cellchat@netP$pathways, function(path) {
  df <- subsetCommunication(cellchat, signaling = path)
  nrow(df)
})
pathway_counts <- sort(pathway_counts, decreasing = TRUE)

# View top 10
head(pathway_counts, 10)

# -------------------------------------------------------------
#Step 7: Export Summary Tables and Prepare for Reporting
# -------------------------------------------------------------
#Step 7.1: Export All Communication Events
# -------------------------------------------------------------
all_comm <- subsetCommunication(cellchat)

write.csv(all_comm, file = "E:/UBC_wang_qn3/sample1a1_cellchat_all_communications.csv", row.names = FALSE)

# -------------------------------------------------------------
#Step 7.2: Export Top 10 Pathways Interactions Separately
# -------------------------------------------------------------
top10_pathways <- names(head(pathway_counts, 10))

for (path in top10_pathways) {
  df <- subsetCommunication(cellchat, signaling = path)
  write.csv(df, file = paste0("E:/UBC_wang_qn3/sample1a1_cellchat_", path, "_interactions.csv"), row.names = FALSE)
}

# -------------------------------------------------------------
# Step 7.3: Export Cluster Role Scores
# -------------------------------------------------------------
# Recompute centrality scores (if not already done)
cellchat <- netAnalysis_computeCentrality(cellchat)

# Generate the scatter plot and extract the centrality scores from it
role_plot <- netAnalysis_signalingRole_scatter(cellchat)
role_df <- role_plot$data

# Save to CSV
write.csv(role_df, file = "E:/UBC_wang_qn3/sample1a1_cellchat_cluster_roles.csv", row.names = FALSE)



















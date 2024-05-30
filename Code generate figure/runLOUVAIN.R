# Load necessary libraries
library(dplyr)
library(tidyr)
library(Seurat)

# Define sources for methodological reference
# Reference links are not directly used in code but provide context and guidelines
# Reference: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# Reference: https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
# Reference: https://satijalab.org/seurat/articles/spatial_vignette

# Load and prepare data
data <- read.csv("path_to_data.csv")
orig_values <- as.matrix(data)  # Consider only marker data
rownames(orig_values) <- 1:nrow(data)
orig_values <- t(orig_values)
orig_values_metadata <- data.frame("CellID" = 1:nrow(data))
rownames(orig_values_metadata) <- orig_values_metadata$CellID

# Create a Seurat object from the dataset
scfp <- CreateSeuratObject(counts = orig_values, meta.data = orig_values_metadata)

# Data preprocessing
# Step 1: Normalize data
scfp <- NormalizeData(scfp, normalization.method = "CLR", margin = 2)

# Step 2: Identify highly variable features
scfp <- FindVariableFeatures(scfp, selection.method = "vst", nfeatures = 50)

# Step 3: Scale data
scfp <- ScaleData(scfp, features = rownames(scfp))

# Step 4: Perform linear dimensional reduction
scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp))
scfp <- RunUMAP(scfp, reduction = "pca", dims = 1:30)

# Step 5: Clustering
scfp <- FindNeighbors(scfp, dims = 1:30)
scfp <- FindClusters(scfp, resolution = 0.8)
clusters <- as.numeric(scfp@meta.data[["seurat_clusters"]])

# Step 6: Find differentially expressed genes
scfp.markers <- FindAllMarkers(scfp, only.pos = TRUE)
top5 <- scfp.markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  as.data.frame()

# Step 7: Assign cell type identity
DoHeatmap(scfp, features = top5$gene) + NoLegend()

# Save clustering results
write.csv(data.frame(Original_clusters = clusters, data), "Louvain_result.csv")






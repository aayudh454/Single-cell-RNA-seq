# --- Set working directory ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")

# --- Load required libraries ---
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)

# --- Define the gene list ---
genes <- c("CD34", "THY1", "PTPRC", "RUNX1", "HOXA9", 
           "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2",
           "PECAM1","CDH5","FLT1","TIE1","TEK","SPARC", "KCNK17", "IL33")

# --- Step 1: Load data and create Seurat objects ---
data_eng1   <- Read10X_h5("ENG-062824_filtered_feature_bc_matrix.h5")
data_gmp1   <- Read10X_h5("GDP-081924_filtered_feature_bc_matrix.h5")
data_pd160  <- Read10X_h5("PD160_EHT4_filtered_feature_bc_matrix.h5")

seu_eng1   <- CreateSeuratObject(counts = data_eng1, project = "ENG1", min.cells = 3, min.features = 200)
seu_gmp1   <- CreateSeuratObject(counts = data_gmp1, project = "GMP1", min.cells = 3, min.features = 200)
seu_pd160  <- CreateSeuratObject(counts = data_pd160, project = "PD160", min.cells = 3, min.features = 200)

# --- Step 2: Merge the Seurat objects with sample IDs ---
seuset <- merge(seu_eng1, y = list(seu_gmp1, seu_pd160),
                add.cell.ids = c("ENG1", "GMP1", "PD160"),
                project = "ST101_Merged")

# --- Step 3: Mitochondrial percentage and filtering ---
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^MT-")
seuset <- subset(
  seuset,
  subset = nFeature_RNA > 600 &
    nFeature_RNA < 5500 &
    nCount_RNA > 2000 &
    percent.mt < 10
)

# --- Step 4: Minimum processing to allow clustering ---
seuset <- NormalizeData(seuset)
seuset <- FindVariableFeatures(seuset)
seuset <- ScaleData(seuset)
seuset <- RunPCA(seuset)
seuset <- FindNeighbors(seuset, dims = 1:30)
seuset <- RunUMAP(seuset, dims = 1:30, verbose = FALSE)
seuset <- RunTSNE(seuset, dims = 1:30, verbose = FALSE)


# --- Step 5: Clustering at both resolutions ---
seuset <- FindClusters(seuset, resolution = 1.0)
cluster_res1 <- seuset$seurat_clusters

seuset <- FindClusters(seuset, resolution = 2.0)
cluster_res2 <- seuset$seurat_clusters

# --- Step 6: Store both clustering results ---
seuset$cluster_res1 <- cluster_res1
seuset$cluster_res2 <- cluster_res2

# --- Step 7: Extract raw counts for selected genes only ---
raw_counts <- LayerData(seuset[["RNA"]], layer = "counts")
genes_present <- intersect(genes, rownames(raw_counts))
subset_counts <- raw_counts[genes_present, ]

# --- Step 8: Transpose to cells x genes ---
cell_gene_matrix <- Matrix::t(subset_counts)

# --- Step 9: Create metadata dataframe with both cluster assignments ---
meta_df <- data.frame(
  cluster_res1 = seuset$cluster_res1,
  cluster_res2 = seuset$cluster_res2
)
rownames(meta_df) <- colnames(seuset)

# --- Step 10: Ensure rownames match and combine ---
meta_df <- meta_df[rownames(cell_gene_matrix), , drop = FALSE]
final_matrix <- cbind(meta_df, as.matrix(cell_gene_matrix))
head(final_matrix)

# --- Step 11: Write to CSV ---
write.csv(final_matrix, "ENG1cell_by_selected_genes_clusters_res1_res2_rawcounts.csv")

# UMAP - res 1.0
p_umap_res1 <- DimPlot(seuset, reduction = "umap", group.by = "cluster_res1", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Clusters (Resolution 1.0)") +
  theme_minimal()

# UMAP - res 2.0
p_umap_res2 <- DimPlot(seuset, reduction = "umap", group.by = "cluster_res2", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Clusters (Resolution 2.0)") +
  theme_minimal()

# t-SNE - res 1.0
p_tsne_res1 <- DimPlot(seuset, reduction = "tsne", group.by = "cluster_res1", label = TRUE, repel = TRUE) +
  ggtitle("t-SNE: Clusters (Resolution 1.0)") +
  theme_minimal()

# t-SNE - res 2.0
p_tsne_res2 <- DimPlot(seuset, reduction = "tsne", group.by = "cluster_res2", label = TRUE, repel = TRUE) +
  ggtitle("t-SNE: Clusters (Resolution 2.0)") +
  theme_minimal()

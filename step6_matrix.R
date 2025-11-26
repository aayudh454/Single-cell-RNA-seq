setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
#library(future)
#plan("multisession", workers = 8)

library(tictoc)
library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(tictoc)
library(svglite)
library(dplyr)
library(reshape2)


# Define files and sample names with updated names
files <- c(
  ENG1 = "ENG-062824_filtered_feature_bc_matrix.h5",
  GMP1 = "GDP-081924_filtered_feature_bc_matrix.h5",
  PD160 = "PD160_EHT4_filtered_feature_bc_matrix.h5"
)
sample_name <- names(files)

# Create Seurat objects
seurat_list <- list()
for (i in seq_along(files)) {
  tmp <- Read10X_h5(paste0("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1/", files[i]))
  seurat_list[[i]] <- CreateSeuratObject(counts = tmp, project = sample_name[i])
}
names(seurat_list) <- sample_name

# Merge objects
seurat <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.id = sample_name)
colnames(seurat@meta.data)[1] <- "sample"

# Set treatment column = sample name (optional but preserves info)
seurat$treatment <- seurat$sample
seurat$sample_norep <- seurat$treatment

metadata <- seurat@meta.data

# QC metrics
seurat$mitoRatio <- PercentageFeatureSet(seurat, pattern = "^MT-") / 100
seurat$riboRatio <- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)
seurat$cell_names <- colnames(seurat)

# Extract metadata
meta <- seurat@meta.data

nFeature_RNA_lowerLim <- 600 
nFeature_RNA_upperLim <- 5500 
nCount_RNA_lowerLim <- 2000
mitoPercent_upperLim <- 10


files <- c(
  ENG1 = "ENG-062824_filtered_feature_bc_matrix.h5",
  GMP1 = "GDP-081924_filtered_feature_bc_matrix.h5",
  PD160 = "PD160_EHT4_filtered_feature_bc_matrix.h5"
)

# Load samples
seuset_list <- lapply(names(files), function(name) {
  path <- files[name]
  message(paste("Reading", name, "from", path))
  data <- Read10X_h5(path)
  seu <- CreateSeuratObject(counts = data, project = name)
  seu$group <- name
  seu$sample <- name
  seu
})
names(seuset_list) <- names(files)

# Merge samples
seuset <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x$sample[1], y$sample[1])), seuset_list)
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^MT-")

# Filter cells
seuset <- subset(seuset, subset = nFeature_RNA > nFeature_RNA_lowerLim & 
                   nFeature_RNA < nFeature_RNA_upperLim &
                   nCount_RNA > nCount_RNA_lowerLim &
                   percent.mt < mitoPercent_upperLim)

# Access the "counts" layer (replaces slot = "counts")
counts_combined <- LayerData(seuset[["RNA"]], layer = "counts")

# Filter genes expressed in more than 5 cells
keep_genes <- rownames(counts_combined)[Matrix::rowSums(counts_combined > 0) > 5]

# Subset Seurat object to keep only those genes
seuset <- subset(seuset, features = keep_genes)

seuset$riboRatio <- PercentageFeatureSet(seuset, pattern = "^RPS|^RPL") / 100
seuset$mitoRatio <- PercentageFeatureSet(seuset, pattern = "^MT-") / 100
seuset$log10GenesPerUMI <- log10(seuset$nFeature_RNA) / log10(seuset$nCount_RNA)
metadata <- seuset@meta.data

# Normalize
library(future)
options(future.globals.maxSize = 10 * 1024^3)

tic()
seuset <- NormalizeData(seuset, verbose = FALSE)
toc()

tic()
seuset <- FindVariableFeatures(seuset, verbose = FALSE)
variable_features <- VariableFeatures(seuset)
print(paste("Found", length(variable_features), "variable genes:"))
print(head(variable_features))
seuset <- ScaleData(seuset, vars.to.regress = "percent.mt", verbose = FALSE)
toc()
#seuset <- SCTransform(seuset, vars.to.regress = "percent.mt", verbose = FALSE)
tic()
seuset <- RunPCA(seuset, verbose = FALSE)
toc()
tic()
seuset <- IntegrateLayers(seuset, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
toc()
tic()
seuset2 <- FindNeighbors(seuset, dims = 1:30, verbose = FALSE)
seuset2 <- FindClusters(seuset2, verbose = FALSE, resolution=1)
seuset2 <- RunUMAP(seuset2, dims = 1:30, verbose = FALSE)
seuset2 <- RunTSNE(seuset2, dims=1:30, verbose = FALSE) # check_duplicates=FALSE



seuset3 <- FindNeighbors(seuset, dims = 1:30, verbose = FALSE)
seuset3 <- FindClusters(seuset3, verbose = FALSE, resolution=2)
seuset3 <- RunUMAP(seuset3, dims = 1:30, verbose = FALSE)
seuset3 <- RunTSNE(seuset3, dims=1:30, verbose = FALSE) # check_duplicates=FALSE

umap_res1 <- DimPlot(seuset2, reduction = "umap", label = TRUE,label.size = 5) + theme_custom + NoLegend() + ggtitle("Res: 1.0; Clusters = 19")
umap_res1

umap_res1 <- DimPlot(object = seuset2, reduction = "umap",label = TRUE,label.size = 5,pt.size = 0.5,repel = TRUE) +
  ggtitle("Res: 1.0; Clusters = 19") + theme_void() +  # removes axis lines and ticks
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),legend.position = "none",plot.margin = margin(5, 5, 5, 5))
umap_res1

umap_res2 <- DimPlot(object = seuset3, reduction = "umap",label = TRUE,label.size = 5,pt.size = 0.5,repel = TRUE) +
  ggtitle("Res: 1.0; Clusters = 36") + theme_void() +  # removes axis lines and ticks
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),legend.position = "none",plot.margin = margin(5, 5, 5, 5))
umap_res2

png("Cluster_1.0_2.0_UMAP.png", width = 13, height = 6.2, units = "in", res = 300)
print(plot_grid(umap_res1,umap_res2,nrow = 1, align = "hv"))
dev.off()

# Get cluster identities from resolution 1.0 and 2.0
cluster_res1 <- seuset2$seurat_clusters
cluster_res2 <- seuset3$seurat_clusters

# Get gene expression values
gene_expr <- FetchData(seuset2, vars = c("CD34", "THY1", "PTPRC", "PECAM1"))

# Build the matrix (all cells in same order assumed)
meta_matrix <- data.frame(
  cell_id = rownames(seuset2@meta.data),
  cluster_res1 = cluster_res1,
  cluster_res2 = cluster_res2,
  CD34 = gene_expr$CD34,
  CD90 = gene_expr$THY1,
  CD45 = gene_expr$PTPRC,
  CD31 = gene_expr$PECAM1,
  row.names = NULL
)
meta_matrix_eng1 <- meta_matrix[grepl("^ENG1_ENG1", meta_matrix$cell_id), ]
# Optional: Check output
head(meta_matrix_eng1)
dim(meta_matrix_eng1)
write.csv(meta_matrix_eng1, "Res1.0_res2.0_ST-101_meta_matrix_ENG1.csv", row.names = FALSE)

#CD34+CD90+
subset_matrix <- meta_matrix_eng1[meta_matrix_eng1$CD34 > 0 & meta_matrix_eng1$CD90 > 0, ]
head(subset_matrix)
dim(subset_matrix)

write.csv(subset_matrix, "CD34+CD90+_Res1.0_res2.0_ST-101_subset_matrix_ENG1.csv", row.names = FALSE)

#CD34+CD90+CD45-
subset_matrix1 <- meta_matrix[
  meta_matrix$CD34 > 0 & 
    meta_matrix$CD90 > 0 & 
    meta_matrix$CD45 == 0, 
]
head(subset_matrix1)
dim(subset_matrix1)

# Keep only rows where cell_id starts with "ENG1_ENG1"
meta_matrix_eng1 <- subset_matrix1[grepl("^ENG1_ENG1", subset_matrix1$cell_id), ]


# --- Define gene list ---
genes <- c("PIEZO1", "CD34", "THY1", "PTPRC", "RUNX1", "HOXA9", 
           "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2", "DLL4", "CXCR4",
           "PECAM1", "CDH5", "FLT1", "TIE1", "TEK", "SPARC", "KCNK17", "IL33")

# --- Fetch gene expression values ---
gene_expr <- FetchData(seuset2, vars = genes)

# --- Build meta_matrix ---
meta_matrix <- data.frame(
  cell_id = rownames(seuset2@meta.data),
  cluster_res1 = seuset2$seurat_clusters,            # or use cluster_res1 if stored separately
  cluster_res2 = cluster_res2,                       # assuming this still exists in your environment
  gene_expr,                                         # include all gene expression values
  row.names = NULL
)

# --- Optional: Preview ---
head(meta_matrix)
# Keep only rows where cell_id starts with "ENG1_ENG1"
meta_matrix_eng1 <- meta_matrix[grepl("^ENG1_ENG1", meta_matrix$cell_id), ]

# Optional: Preview filtered matrix
head(meta_matrix_eng1)

write.csv(meta_matrix_eng1, "HSC_HE_Res1.0_res2.0_ST-101_meta_matrix.csv", row.names = FALSE)

#CD34+CD90+
# Filter rows where both CD34 > 0 and CD90 > 0
meta_matrix_filtered <- meta_matrix_eng1[
  meta_matrix_eng1$CD34 > 0 & meta_matrix_eng1$CD90 > 0,
]


write.csv(subset_matrix, "CD34+CD90+_Res1.0_res2.0_ST-101_subset_matrix_ENG1.csv", row.names = FALSE)

########-----------------------CLUSTER genes

# Step 1: Subset Seurat object to cluster 9 and 14 (from Res 1.0)
seuset2$cluster_res1 <- cluster_res1  # ensure cluster_res1 is stored
cells_to_keep <- rownames(seuset2@meta.data)[seuset2$cluster_res1 %in% c(9, 14)]
seu_subset <- subset(seuset2, cells = cells_to_keep)

# Step 2: Extract all gene expression values (RNA assay, data slot)
all_genes <- rownames(seu_subset)
expr_matrix <- FetchData(seu_subset, vars = all_genes)

# Step 3: Keep only genes with at least one non-zero value
expr_matrix_filtered <- expr_matrix[, colSums(expr_matrix > 0) > 0]

# Optional: Add metadata (cluster ID and cell_id)
expr_matrix_filtered$cell_id <- rownames(expr_matrix_filtered)
expr_matrix_filtered$cluster_res1 <- seu_subset$cluster_res1

# Reorder columns: cell_id, cluster, then gene values
expr_matrix_filtered <- expr_matrix_filtered[, c("cell_id", "cluster_res1", setdiff(colnames(expr_matrix_filtered), c("cell_id", "cluster_res1")))]

# View output
head(expr_matrix_filtered)
dim(expr_matrix_filtered)

expr_matrix_filtered_eng1 <- expr_matrix_filtered[grepl("^ENG1_ENG1", expr_matrix_filtered$cell_id), ]
dim(expr_matrix_filtered_eng1)
head(expr_matrix_filtered_eng1)

# --- Step 1: Subset matrix to cluster 9 ---
cluster9_matrix <- expr_matrix_filtered_eng1[expr_matrix_filtered_eng1$cluster_res1 == 9, ]

# --- Step 2: Remove metadata columns ---
gene_matrix <- cluster9_matrix[, !(colnames(cluster9_matrix) %in% c("cell_id", "cluster_res1"))]

# --- Step 3: Calculate max expression per gene ---
max_vals <- apply(gene_matrix, 2, max)

# --- Step 4: Sort and get top 20 genes ---
top100_genes <- sort(max_vals, decreasing = TRUE)[1:100]

# --- Step 5: Print result ---
top100_genes

#####-----------Cluster14
# --- Step 1: Subset matrix to cluster 14 ---
cluster14_matrix <- expr_matrix_filtered_eng1[expr_matrix_filtered_eng1$cluster_res1 == 14, ]

# --- Step 2: Remove metadata columns ---
gene_matrix_14 <- cluster14_matrix[, !(colnames(cluster14_matrix) %in% c("cell_id", "cluster_res1"))]

# --- Step 3: Calculate max expression per gene ---
max_vals_14 <- apply(gene_matrix_14, 2, max)

# --- Step 4: Sort and get top 100 genes ---
top100_genes_14 <- sort(max_vals_14, decreasing = TRUE)[1:100]

# --- Step 5: Print result ---
top100_genes_14


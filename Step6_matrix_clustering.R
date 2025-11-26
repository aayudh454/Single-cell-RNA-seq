# --- Load Required Library ---
library(pheatmap)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(devtools)
library(tictoc)
library(patchwork)

setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
meta_matrix_eng1 <- read.csv("Res1.0_res2.0_ST-101_meta_matrix_ENG1.csv")
head(meta_matrix_eng1)
list.files()

# --- Step 1: Extract and Order Marker Expression Matrix ---
expr_data <- meta_matrix_eng1[, c("CD34", "CD90", "CD45", "CD31")]
rownames(expr_data) <- meta_matrix_eng1$cell_id

# --- Step 2: Create Annotation Data (only Cluster_Res1) ---
annotation_row <- data.frame(
  Cluster_Res1 = as.factor(meta_matrix_eng1$cluster_res1)
)
rownames(annotation_row) <- meta_matrix_eng1$cell_id

# --- Step 3: Convert to numeric and clean ---
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- meta_matrix_eng1$cell_id

# Remove rows with any NA/NaN/Inf
keep_rows <- apply(expr_data, 1, function(x) all(is.finite(x)))
expr_data <- expr_data[keep_rows, ]
annotation_row <- annotation_row[rownames(expr_data), , drop = FALSE]

# Remove constant rows (zero standard deviation)
non_constant_rows <- apply(expr_data, 1, function(x) sd(x) > 0)
expr_data_clean <- expr_data[non_constant_rows, ]
annotation_row_clean <- annotation_row[rownames(expr_data_clean), , drop = FALSE]

# --- Step 4: Final heatmap without row names, ordered columns ---
# --- Step 4: Final heatmap without row names, ordered columns ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")

# Open PNG device
png("CD34_cd90_cd45_hierarch_res_1.png", width = 7, height = 6.2, units = "in", res = 300)

# Generate heatmap
print(pheatmap(
  mat = expr_data_clean[, c("CD34", "CD90", "CD45", "CD31")],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row_clean,
  scale = "row",
  show_rownames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize = 10,
  fontsize_col = 12,
  border_color = NA,
  main = "Hierarchical Clustering of Marker Expression"
))

# Close the PNG device to save the file
dev.off()

#####CD34CD90+--------------------------------------------------------------------
# --- Load Required Library ---
library(pheatmap)

# --- Step 1: Extract and Order Marker Expression Matrix ---
expr_data <- subset_matrix[, c("CD34", "CD90", "CD45", "CD31")]
rownames(expr_data) <- subset_matrix$cell_id

# --- Step 2: Create Annotation Data (only Cluster_Res1) ---
annotation_row <- data.frame(
  Cluster_Res1 = as.factor(subset_matrix$cluster_res1)
)
rownames(annotation_row) <- subset_matrix$cell_id

# --- Step 3: Convert to numeric and clean ---
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- subset_matrix$cell_id

# Remove rows with any NA/NaN/Inf
keep_rows <- apply(expr_data, 1, function(x) all(is.finite(x)))
expr_data <- expr_data[keep_rows, ]
annotation_row <- annotation_row[rownames(expr_data), , drop = FALSE]

# Remove constant rows (zero standard deviation)
non_constant_rows <- apply(expr_data, 1, function(x) sd(x) > 0)
expr_data_clean <- expr_data[non_constant_rows, ]
annotation_row_clean <- annotation_row[rownames(expr_data_clean), , drop = FALSE]

# --- Step 4: Final heatmap without row names, ordered columns ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("CD34_cd90_cd45_hierarch_res_1.png", width = 7, height = 6.2, units = "in", res = 300)
print(
  pheatmap(
    mat = expr_data_clean[, c("CD34", "CD90", "CD45", "CD31")],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_row = annotation_row_clean,
    scale = "row",
    show_rownames = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize = 10,
    fontsize_col = 12,
    border_color = NA,
    main = "Hierarchical Clustering CD34+CD90+ ST-101"
  )
)
dev.off()


####----------------all
# --- Load Required Library ---
library(pheatmap)

# --- Step 1: Extract and Order Marker Expression Matrix ---
expr_data <- meta_matrix_eng1[, c("PIEZO1", "CD34", "THY1", "PTPRC", "RUNX1", "HOXA9", 
                                  "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2", "DLL4", "CXCR4",
                                  "PECAM1", "CDH5", "FLT1", "TIE1", "TEK", "SPARC", "KCNK17", "IL33")]
rownames(expr_data) <- meta_matrix_eng1$cell_id

# --- Step 2: Create Annotation Data (only Cluster_Res1) ---
annotation_row <- data.frame(
  Cluster_Res1 = as.factor(meta_matrix_eng1$cluster_res1)
)
rownames(annotation_row) <- meta_matrix_eng1$cell_id

# --- Step 3: Convert to numeric and clean ---
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- meta_matrix_eng1$cell_id

# Remove rows with any NA/NaN/Inf
keep_rows <- apply(expr_data, 1, function(x) all(is.finite(x)))
expr_data <- expr_data[keep_rows, ]
annotation_row <- annotation_row[rownames(expr_data), , drop = FALSE]

# Remove constant rows (zero standard deviation)
non_constant_rows <- apply(expr_data, 1, function(x) sd(x) > 0)
expr_data_clean <- expr_data[non_constant_rows, ]
annotation_row_clean <- annotation_row[rownames(expr_data_clean), , drop = FALSE]

# --- Step 4: Final heatmap without row names, ordered columns ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("AllHSC_HE_hierarch_res_1.png", width = 7, height = 6.2, units = "in", res = 300)
print(
  pheatmap(
    mat = expr_data_clean[, c("PIEZO1", "CD34", "THY1", "PTPRC", "RUNX1", "HOXA9", 
                              "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2", "DLL4", "CXCR4",
                              "PECAM1", "CDH5", "FLT1", "TIE1", "TEK", "SPARC", "KCNK17", "IL33")],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_row = annotation_row_clean,
    scale = "row",
    show_rownames = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize = 10,
    fontsize_col = 12,
    border_color = NA,
    main = "Hierarchical Clustering of Marker Expression"
  )
)
dev.off()

library(pheatmap)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(devtools)
library(tictoc)
library(patchwork)

# --- Paths & Input ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
meta_matrix_all <- read.csv("Res1.0_res2.0_ST-101_meta_matrix.csv", stringsAsFactors = FALSE)
head(meta_matrix_all)
#CD34+CD90+
# Filter rows where both CD34 > 0 and CD90 > 0
meta_matrix_all <- meta_matrix_all[
  meta_matrix_all$CD34 > 0 & meta_matrix_all$CD90 > 0,
]

# --- Safety checks (required columns) ---
required_cols <- c("cell_id", "cluster_res1", "CD34", "CD90", "CD45", "CD31")
missing <- setdiff(required_cols, colnames(meta_matrix_all))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

# --- Step 1: Expression matrix (desired marker order) ---
markers   <- c("CD45","CD31","CD34","CD90")
expr_data <- meta_matrix_all[, markers]
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- meta_matrix_all$cell_id

# --- Step 2: Row annotation with explicit 0:19 cluster ordering ---
annotation_row <- data.frame(
  Cluster_Res1 = factor(meta_matrix_all$cluster_res1, levels = 0:19)
)
rownames(annotation_row) <- meta_matrix_all$cell_id

# --- Step 3: Clean rows: finite values only, remove zero-variance rows ---
keep_rows <- apply(expr_data, 1, function(x) all(is.finite(x)))
expr_data <- expr_data[keep_rows, , drop = FALSE]
annotation_row <- annotation_row[base::rownames(expr_data), , drop = FALSE]

non_constant_rows <- apply(expr_data, 1, function(x) stats::sd(x) > 0)
expr_data_clean   <- expr_data[non_constant_rows, , drop = FALSE]
annotation_row_clean <- annotation_row[base::rownames(expr_data_clean), , drop = FALSE]

# --- Step 4: Order rows by cluster_res1 (0 → 19) and compute gaps between clusters ---
ord_idx <- order(annotation_row_clean$Cluster_Res1, na.last = TRUE)
expr_ord <- expr_data_clean[ord_idx, , drop = FALSE]
ann_ord  <- annotation_row_clean[ord_idx, , drop = FALSE]

clust_counts <- table(droplevels(ann_ord$Cluster_Res1))
gaps_row <- if (length(clust_counts) > 0) {
  gr <- cumsum(clust_counts)
  if (length(gr) > 1) gr[-length(gr)] else NULL
} else NULL

# --- Step 5: Build row labels that show cluster name only at the start of each block ---
cluster_vec <- as.character(ann_ord$Cluster_Res1)
starts <- which(!duplicated(cluster_vec))       # first row of each cluster block
labels_row <- rep("", nrow(ann_ord))
labels_row[starts] <- paste0("Cluster ", cluster_vec[starts])

# --- Step 6: Plot (no dendrograms), landscape PNG, cluster names on left ---
 
png("3rep_CD_markers_res1_CD34_CD90.png", width = 7, height = 10, units = "in", res = 300)

pheatmap(
  mat = expr_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_row = labels_row,        # <-- cluster names on left at block starts
  show_rownames = TRUE,
  annotation_row = ann_ord,
  annotation_names_row = FALSE,   # hide "Cluster_Res1" text label
  annotation_legend = FALSE,      # optional: hide legend for the color bar
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize_row = 8,
  fontsize_col = 12,
  border_color = NA,
  main = "Hierarchial Clustering CD34+CD90+ ST-101",
  gaps_row = gaps_row,
  treeheight_row = 0,
  treeheight_col = 0
)

dev.off()
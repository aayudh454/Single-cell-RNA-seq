# --- Load library ---
library(ggplot2)
library(tidyr)
library(dplyr)

# --- Step 1: Set working directory ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")

# --- Step 2: Read datasets ---
HSC <- read.csv("clusters_9_14_vs_rest_significant.csv")
Progenitors1 <- read.csv("clusters_6_7_10_vs_rest_significant_CD34pos_progenitors.csv")
Progenitors2 <- read.csv("clusters_0,5,4_vs_rest_significant_CD34neg.csv")
Stroma <- read.csv("clusters_1_13_16_18_vs_rest_significant_stromal_cells.csv")

# --- Step 3: Define full gene list (explicit) ---
gene_list <- unique(c(
  "GJA5","IL33","ADGRF5","CLEC14A","DACT2","CLDN5","ROBO4","CDH5","SOX17",
  "KDR","HEY2","NOTCH4","SOX7","EFNB2","TEK","FOXC1","GJA4","HEY1","MECOM",
  "EPAS1","CD93","ESAM","ENG","FLT1","KITLG","NOTCH1","PECAM1","CD34","SPARC",
  "MLLT3","GATA2","RUNX1","SPINK2","PTPRC","CD38","MYCN","KIT","ITGA2B","HLA-DRA", "HLA-B", "HLA-C",
  "FN1","HSPG2","GPC6","ID3","PIEZO1","NT5E","MCAM","ALCAM","NFGR",
  "CXCL12","PDGFRA","SPP1","THY1","PROCR","HOXA9","DLL4","CXCR4","TIE1","KCNK17"
))

# --- Step 4: Define priority gene order ---
priority_order <- c(
  "PIEZO1", "CD34", "THY1", "PROCR", "PTPRC", "RUNX1", "HOXA9", "SPINK2",
  "MECOM", "MLLT3", "KIT", "GATA2", "DLL4", "CXCR4",
  "FN1", "GPC6", "HSPG2", "ID3",
  "PECAM1","CDH5","FLT1","TIE1","TEK","SPARC","KCNK17","IL33",
  "CD38","MYCN","ESAM","GATA2","KIT","ITGA2B","HLA-DRA", "HLA-B", "HLA-C",
  "ENG","NT5E","MCAM","ALCAM","NFGR","CXCL12","PDGFRA","SPP1"
)

# --- Step 5: Create a helper function to process each dataset ---
process_dataset <- function(df, gene_list, priority_order, name_col) {
  df_filtered <- df[, c("gene", "avg_log2FC")]
  gene_df <- data.frame(gene = gene_list, stringsAsFactors = FALSE)
  
  # Merge to include all genes
  m <- merge(gene_df, df_filtered, by = "gene", all.x = TRUE)
  m$avg_log2FC[is.na(m$avg_log2FC)] <- 0
  
  # Reorder genes: priority first, then remaining
  other_genes <- setdiff(m$gene, priority_order)
  final_order <- c(priority_order, other_genes)
  m <- m[match(final_order, m$gene), ]
  rownames(m) <- NULL
  
  # Rename the column to the dataset name
  colnames(m)[2] <- name_col
  return(m)
}

# --- Step 6: Process each dataset ---
HSC_std <- process_dataset(HSC, gene_list, priority_order, "HSC_log2FC")
P1_std  <- process_dataset(Progenitors1, gene_list, priority_order, "Progenitors1_log2FC")
P2_std  <- process_dataset(Progenitors2, gene_list, priority_order, "Progenitors2_log2FC")
STR_std <- process_dataset(Stroma, gene_list, priority_order, "Stroma_log2FC")

# --- Step 7: Merge all processed datasets into one final table ---
combined <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
                   list(HSC_std, P1_std, P2_std, STR_std))

# Replace remaining NAs with 0 (safety)
num_cols <- setdiff(colnames(combined), "gene")
for (cn in num_cols) combined[[cn]][is.na(combined[[cn]])] <- 0

# Keep custom order globally
other_genes_all <- setdiff(combined$gene, priority_order)
final_order_all <- c(priority_order, other_genes_all)
combined <- combined[match(final_order_all, combined$gene), ]
rownames(combined) <- NULL

# --- Step 8: Preview and save final dataset ---
head(combined, 20)

# --- Optional: Save the ordered dataset ---
#write.csv(combined, "All_clusters_combined_log2FC.csv", row.names = FALSE)

####----dotplot
gene_data <- read.csv("All_clusters_combined_log2FC_function.csv")
head(gene_data)

gene_data$Gene <- ifelse(
  is.na(gene_data$Function) | trimws(gene_data$Function) == "",
  gene_data$Gene,
  paste0(gene_data$Gene, " | ", trimws(gene_data$Function))
)

gene_data <- gene_data[, !(colnames(gene_data) %in% c("Function"))]
head(gene_data)

library(ggplot2)
library(tidyr)


gene_data_long <- gene_data %>%
  pivot_longer(
    cols = ends_with("_log2FC"),
    names_to = "Group",
    values_to = "log2FC"
  ) %>%
  mutate(Gene_label = sub(" \\|.*", "", Gene),
         Gene_label = factor(Gene_label, levels = rev(unique(Gene_label))))

pecam_pos <- which(levels(gene_data_long$Gene_label) == "ID3") - 0.5
il33_pos <- which(levels(gene_data_long$Gene_label) == "IL33") - 0.5
hlac_pos <- which(levels(gene_data_long$Gene_label) == "HLA-C") - 0.5
spp1_pos <- which(levels(gene_data_long$Gene_label) == "SPP1") - 0.5

png("gene_dotplot_single_alls.png", width = 10, height = 10, units = 'in', res = 300)
ggplot(gene_data_long, aes(x = Group, y = Gene_label, fill = log2FC, size = abs(log2FC))) +
  geom_point(shape = 21, color = "black", stroke = 0.7) +
  # Horizontal black line between ID3 and PECAM1
  geom_hline(yintercept = pecam_pos, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = il33_pos, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = hlac_pos, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = spp1_pos, color = "black", linewidth = 0.6) +
  scale_fill_gradient2(low = "#4BC9FF", mid = "#A4BDBC", high = "#D50032", midpoint = 0) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c(
    "HSC_log2FC" = "HSC+HE",
    "Progenitors1_log2FC" = "CD34+",
    "Progenitors2_log2FC" = "CD34-",
    "Stroma_log2FC" = "Stroma"
  )) +
  labs(x = NULL, y = NULL, fill = "log2FC", size = "|log2FC|") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(family = "Arial", size = 10, face = "italic", color = "black"),
    legend.position = "left",
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(color = "black", face = "bold", size = 14),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )
dev.off()
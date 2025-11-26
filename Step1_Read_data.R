setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
library(future)
plan("multisession", workers = 4)

list.files()
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(extrafont)

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
  tmp <- Read10X_h5(paste0("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/Only_ST_101/", files[i]))
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
metadata <- seurat@meta.data

# Define custom colors
custom_colors <- c(
  ENG1 = "#A4BDBC",
  GMP1 = "#D50032",
  PD160 = "#4BC9FF"
)

metadata$sample <- factor(metadata$sample, levels = c("PD160", "ENG1", "GMP1"))


# Combined plots
# Combined plots
combined_plot <- (
  (metadata %>% 
     ggplot(aes(x = sample, fill = sample)) + 
     geom_bar() +
     scale_fill_manual(values = custom_colors) +
     theme_classic() +
     theme(
       axis.ticks.x = element_blank(),
       axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),  # No angle specified
       axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
       axis.title.x = element_blank(),
       axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
       plot.title = element_text(hjust = 0.5, face = "bold", color = "black")
     ) +
     ggtitle("Cell count") +
     guides(fill = "none"))
  
  +
    (metadata %>% 
       ggplot(aes(x = nFeature_RNA, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       scale_x_log10() +
       theme_classic() +
       ylab("Cell density") +
       geom_vline(xintercept = 1200) +
       theme(
         axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black", angle = 0, hjust = 0.5),
         axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         plot.title = element_text(family = "Arial", size = 14, face = "bold", color = "black", hjust = 0.5)
       )
     +
       ggtitle("Number of detected genes/cell") +
       NoLegend())
  
  +
    (metadata %>% 
       ggplot(aes(x = nCount_RNA, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       scale_x_log10() +
       theme_classic() +
       ylab("Cell density") +
       geom_vline(xintercept = 2500) +
       theme(
         axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         plot.title = element_text(family = "Arial", size = 14, face = "bold", color = "black", hjust = 0.5)
       ) +
       ggtitle("RNA molecules detected/cell") +
       NoLegend())
  
  +
    (metadata %>%
       filter(mitoRatio > 0) %>%
       ggplot(aes(x = mitoRatio, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_vline(xintercept = 0.075) +
       theme(
         axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         plot.title = element_text(family = "Arial", size = 14, face = "bold", color = "black", hjust = 0.5)
       ) +
       ggtitle("Reads mapped to MT genes") +
       NoLegend())
  +
    (metadata %>% 
       filter(riboRatio > 0) %>%  
       ggplot(aes(x = riboRatio, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_vline(xintercept = 0.1) +
       ylab("Cell density") +
       xlab("riboRatio") +
       theme(
         axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         plot.title = element_text(family = "Arial", size = 14, face = "bold", color = "black", hjust = 0.5)
       ) +
       ggtitle("Reads mapped to ribosomal genes") +
       NoLegend())
  +
    (metadata %>% 
       ggplot(aes(x = sample, y = log10GenesPerUMI, fill = sample)) + 
       geom_violin() +
       geom_boxplot(width = 0.1, fill = alpha("white", 0.7)) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_hline(yintercept = 0.8) +
       ylab("log10GenesPerUMI") +
       xlab(NULL) +
       theme(
         axis.ticks.x = element_blank(),
         axis.text.x = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.text.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         axis.title.y = element_text(family = "Arial", size = 14, face = "bold", color = "black"),
         plot.title = element_text(family = "Arial", size = 14, face = "bold", color = "black", hjust = 0.5)
       ) +
       ggtitle("Library complexity") +
       NoLegend())
  
)

# Save plot
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("rawdata_ST101_only_all.png", width = 13, height = 6, units = 'in', res = 300)
print(combined_plot)
dev.off()

# Filter cells
keep_cell <- metadata$cell_names[
  metadata$nFeature_RNA > 1200 &
    metadata$nCount_RNA > 2500 &
    metadata$mitoRatio < 0.075
]
seurat <- seurat[, keep_cell]

# Combine counts from all layers
layer_names <- names(seurat[["RNA"]]@layers)
counts_list <- lapply(layer_names, function(ly) LayerData(seurat[["RNA"]], ly))
counts_combined <- do.call(cbind, counts_list)

# Filter genes expressed in >5 cells
keep_genes <- rownames(counts_combined)[Matrix::rowSums(counts_combined > 0) > 5]
seurat <- subset(seurat, features = keep_genes)

# Save final Seurat object
saveRDS(seurat, "step1_seurat_filtered_ST101_only_v5.RDS")

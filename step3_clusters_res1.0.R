setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
#library(future)
#options(future.globals.maxSize = 20 * 1024^3)


#load("step2_filtered_normalized_ST-101_res1.RData")
library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(tictoc)
library(svglite)
library(dplyr)
library(reshape2)

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
toc()
#ElbowPlot(seuset, reduction = "pca", ndims = 30)
#ElbowPlot

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.x = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.y = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", colour = "black", size = 14)
  )
custom_colors <- c(
  ENG1 = "#A4BDBC",
  GMP1 = "#D50032",
  PD160 = "#4BC9FF"
)

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

umap_res1 <- DimPlot(seuset2, reduction = "umap", label = TRUE,label.size = 6,repel = TRUE) + theme_custom + NoLegend() + ggtitle("Res: 1.0; Clusters = 19")
umap_res1

png("Cluster_1.0_UMAP.png", width = 5, height = 5, units = "in", res = 300)
print(plot_grid(umap_res1,nrow = 1, align = "hv"))
dev.off()

tsne_res1 <- DimPlot(seuset2, reduction = "tsne", label = TRUE,label.size = 6,repel = TRUE) + 
  theme_custom + NoLegend() + ggtitle("Res: 1.0; Clusters = 19")
tsne_res1 

png("Cluster_1.0_UMAP_tNSE.png", width = 13, height = 6.2, units = "in", res = 300)
print(plot_grid(umap_res1,tsne_res1, nrow = 1, align = "hv"))
dev.off()

umap_res2 <- DimPlot(seuset3, reduction = "umap", label = TRUE,label.size = 6,repel = TRUE) +  theme_custom +NoLegend() + ggtitle("Res: 2.0; Clusters = 36")
umap_res2

tsne_res2 <- DimPlot(seuset3, reduction = "tsne", label = TRUE,label.size = 6,repel = TRUE) + theme_custom + NoLegend() + ggtitle("Res: 2.0; Clusters = 36")
tsne_res2 

png("Cluster_1.0_2.0_UMAP_tNSE.png", width = 10, height = 6.2, units = "in", res = 300)
print(plot_grid(umap_res1,tsne_res1,umap_res2,tsne_res2, nrow = 2, align = "hv"))
dev.off()

#~~~~~~~~~PCA plots----------------------------------------------------------

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.x = element_text(family = "Arial", face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(family = "Arial", face = "bold", colour = "black", size = 12),
    axis.ticks = element_line(colour = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(colour = "black", size = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", colour = "black", size = 14)
  )

# Extract PCA variance explained
pca_variance <- Stdev(seuset2, reduction = "pca")^2
pca_percent_var <- round(100 * pca_variance / sum(pca_variance), 1)

# Custom labels with variance
x_label <- paste0("PC1 (", pca_percent_var[1], "%)")
y_label <- paste0("PC2 (", pca_percent_var[2], "%)")

# Generate DimPlot with labeled axes
umap_res1_pca <- DimPlot(seuset2, reduction = "pca", label = TRUE) +
  xlab(x_label) +
  ylab(y_label) +
  theme_custom +
  NoLegend()
umap_res1_pca 

# Extract standard deviations of PCs
sdev <- seuset@reductions$pca@stdev

# Calculate % variance explained
var_exp <- (sdev^2 / sum(sdev^2)) * 100

# Keep first 30 components
var_exp_30 <- var_exp[1:30]

# Replace NA values (if any) with 0
var_exp_30[is.na(var_exp_30)] <- 0

# Create data frame
pc_df <- data.frame(
  PC = factor(paste0("PC", 1:30), levels = paste0("PC", 1:30)),
  VarianceExplained = var_exp_30
)

# Plot
library(ggplot2)

bar_plot <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "#4BC9FF") +
  geom_text(aes(label = sprintf("%.1f", VarianceExplained)),
            vjust = -0.3, fontface = "bold", size = 3.5, family = "Arial") +
  theme_classic() +
  labs(
    title = NULL,
    x = "Principal Components",
    y = "Variance Explained (%)"
  ) +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


combined_plot <- plot_grid(umap_res1_pca, bar_plot, ncol = 2, rel_widths = c(1, 1.5))
combined_plot

setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("Combined_PCA_Plots.png", width = 13, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

#umap_p4 <- DimPlot(seuset, reduction = "umap", label = TRUE, group.by = "group", split.by = "sample", cols = custom_colors) + theme_custom
umap_p5 <- DimPlot(seuset2, reduction = "umap", label = TRUE,label.size = 3, group.by = "sample", split.by = "group", cols = custom_colors) + theme_custom
umap_p5
png("Sample_res1.0_UMAP.png", width = 6.5, height = 6, units = "in", res = 300)
print(umap_p5)
dev.off()

#save.image("step3_cluster_ST-101_res1.RData")

setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(devtools)
library(tictoc)
library(patchwork)
#plan("multisession", workers = 4)

load("step3_cluster_ST-101_res1.RData")

# join layers
seuset2 <- JoinLayers(seuset2)

tic()
#seuset <- PrepSCTFindMarkers(seuset)
##is identifying marker genes for each cluster in your Seurat single-cell object (seuset) 
##by performing differential expression analysis.
all.markers <- FindAllMarkers(seuset2, only.pos = TRUE)
toc()

#save.image("step4_ST-101_session.RData")

load("step4_ST-101_session.RData")


# Custom theme
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.x = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.y = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", colour = "black", size = 14)
  )

# Define custom colors
custom_colors <- c(
  ENG1 = "#A4BDBC",
  GMP1 = "#D50032",
  PD160 = "#4BC9FF"
)

gene <- "MYCN"

# Fetch expression and z-score normalize
expr <- FetchData(seuset2, vars = gene)
expr$zscore <- scale(expr[[gene]])[,1]
seuset2[[paste0(gene, "_zscore")]] <- expr$zscore

# Generate UMAP plot with custom styling
 FeaturePlot(object = seuset2,features = paste0(gene, "_zscore"),
  cols = c("grey90", "#D50032"),min.cutoff = -2,max.cutoff = 2,label = TRUE,order = TRUE,raster = TRUE,
  pt.size = 2) + theme_custom + theme(legend.position = "none") +ggtitle(gene)

 
# Individual plots
c <- FeaturePlot(seuset2, features = gene, label = TRUE, order = TRUE, raster = TRUE,
                 min.cutoff = "q10", max.cutoff = "q90") + theme_custom
c
a <- VlnPlot(seuset2, features = gene, split.by = "group") + 
  scale_fill_manual(values = custom_colors) + theme_custom
b <- RidgePlot(seuset2, features = gene, ncol = 1) + theme_custom
d <- FeaturePlot(seuset2, features = gene, label = TRUE, order = TRUE, raster = TRUE,
                 min.cutoff = "q10", max.cutoff = "q90", split.by = "group", ncol = 1) + theme_custom


# Save as PNG
png("scRNA2025_CD34_featurePlot.png", width = 36, height = 8, units = "in", res = 300)
print(e)
dev.off()

#####--------------------------------------------------------------------------------------------
# Load required libraries
library(Seurat)
library(ggplot2)


# Define custom clean theme
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# Define gene list
genes <- c("CD34", "THY1", "PTPRC", "RUNX1", "HOXA9", 
           "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2", "DLL4", "CXCR4")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, nrow = 2)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("HSC markers_res_1.png", width = 9, height = 5, units = "in", res = 300)
print(combined_plot)
dev.off()

####------------------------------------------------------------------
# Define gene list
genes <- c("CD34", "THY1", "PTPRC", "PECAM1")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("CD34_Cd90_cd45_cd31_res_1.png", width = 6.5, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

####------------------------------------------------------------------
# Define gene list
genes <- c("PECAM1","CDH5","FLT1","TIE1","TEK","SPARC", "KCNK17", "IL33")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, nrow = 2)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("HE_res_1.png", width = 10, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

####------------------------------------------------------------------
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
# Define gene list
genes <- c("PIEZO1","CD34", "THY1","PROCR", "PTPRC", "RUNX1", "HOXA9", 
           "SPINK2", "MECOM", "MLLT3", "KIT", "GATA2","DLL4", "CXCR4",
           "PECAM1","CDH5","FLT1","TIE1","TEK","SPARC", "KCNK17", "IL33")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    label.size = 3, 
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, nrow = 3)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("HSC_HE_res_1.png", width = 9, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

####------------------Cluster4/9 top genes------------------------------------------------
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
# Define gene list
genes <- c("CD34", "THY1", "PTPRC","PECAM1","FN1","GPC6","HSPG2","ID3")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    label.size = 3, 
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, nrow = 2)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("HSC_HE_res_1.png", width = 9, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

####------------------Progenitors------------------------------------------------
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
# Define gene list
genes <- c("CD34","PTPRC","CD38","MYCN", "ESAM","GATA2", "KIT","ITGA2B")

# Initialize plot list
plots <- list()

# Loop through genes
for (gene in genes) {
  # Fetch expression and z-score normalize
  expr <- FetchData(seuset2, vars = gene)
  expr$zscore <- scale(expr[[gene]])[,1]
  seuset2[[paste0(gene, "_zscore")]] <- expr$zscore
  
  # Generate UMAP plot with custom styling
  p <- FeaturePlot(
    object = seuset2,
    features = paste0(gene, "_zscore"),
    cols = c("grey90", "#D50032"),
    min.cutoff = -2,
    max.cutoff = 2,
    label = TRUE,  # show cluster numbers
    label.size = 3, 
    order = TRUE,
    raster = TRUE,
    pt.size = 2
  ) + 
    theme_custom +
    theme(legend.position = "none") +
    ggtitle(gene)
  
  # Store plot
  plots[[gene]] <- p
}

# Combine plots in a 2x5 grid
combined_plot <- wrap_plots(plots, nrow = 2)
combined_plot

# Save as PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("Progenitors_res_1.png", width = 9, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()
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

####------------------------------Calvanese suppl fig--------------------------------
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
# Define gene list
genes <- c(
  "RUNX1", "HOXA9", "MLLT3", "PTPRC", "SPN", "THY1", "SPINK2", "GFI1",
  "SELP", "STAT5A", "ITGA4", "SVOPL", "EMCN", "ACE", "PROCR", "HOXB9", "LIN28B",
  "CSF1R", "IL3RA", "HLA-DRA", "SELL", "MSI2", "PROM1", "CDH5", "IGFBP2",
  "MEIS2", "NRP2", "SOX17", "GJA5", "IL33", "DKK1", "ALDH1A1", "CD44",
  "KCNK17", "BCL11A", "LIN28A", "GAD1", "FGF23"
)

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

####------------MHC classI and II-------------------------------

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
genes <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M",
  "HLA-DMA", "HLA-DPB1", "HLA-DRA","HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "HLA-DQB1" 
)

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
png("MHCclassI_ClassII_res_1.png", width = 13, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()


##########-------stroma-----------------
theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

genes <- c(
  # MSC / stromal markers
  "THY1",
  "ENG",        # CD105 (Endoglin)
  "NT5E",       # CD73
  "ALCAM",      # CD166
  "ITGB1",      # CD29
  "MCAM",       # CD146
  "NGFR",       # CD271 (LNGFR)
  
  # HSC niche / stromal support markers
  "CXCL12",     # SDF1 chemokine
  "PDGFRA",     # Perivascular stromal marker
  "SPP1"       # Osteopontin
)

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

setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("stromal_markers_res_1.png", width = 10, height = 6.2, units = "in", res = 300)
print(combined_plot)
dev.off()

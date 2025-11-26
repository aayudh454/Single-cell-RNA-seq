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

qc_plot <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.4, color = "#D50032") +
  geom_hline(yintercept = 4200, color = "black", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 5400, color = "blue", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 600, color = "blue", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 1000, color = "black", linetype = "dashed", size = 0.6) +
  annotate("text", x = 18000, y = 4400, 
           label = "nFeature_RNA = 4200", hjust = 0, size = 3.5, fontface = "bold") +
  annotate("text", x = 18000, y = 1200, 
           label = "nFeature_RNA = 1000", hjust = 0, size = 3.5, fontface = "bold") +
  theme_classic() +
  xlab("Total RNA Counts (nCount_RNA)") +
  ylab("Number of Genes Detected (nFeature_RNA)") +
  ggtitle("nFeature_RNA vs nCount_RNA") +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5)
  )

qc_plot


png("QCplot_nFeature_RNA vs nCount_RNA.png", width = 6, height = 6, units = "in", res = 300)
print(qc_plot)
dev.off()

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    #axis.text.x = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    #axis.text.y = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", colour = "black", size = 14)
  )

p6 <- meta %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 300) +
  facet_wrap(~sample)+theme_custom
p6

png("QCplot_mitoRatio.png", width = 8, height = 6, units = "in", res = 300)
print(p6)
dev.off()

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

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.x = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    axis.text.y = element_text(family = "Arial", face = "bold", colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", colour = "black", size = 14)
  )

# Define colors if not already
custom_colors <- c(
  ENG1 = "#A4BDBC",
  GMP1 = "#D50032",
  PD160 = "#4BC9FF"
)

metadata$sample <- factor(metadata$sample, levels = c("PD160", "ENG1", "GMP1"))

# Create all QC plots
library(ggplot2)
library(cowplot)

combined_plot <- (
  (metadata %>% 
     ggplot(aes(x = sample, fill = sample)) + 
     geom_bar() +
     scale_fill_manual(values = custom_colors) +
     theme_classic() +
     theme_custom +
     ggtitle("Cell count") +
     guides(fill = "none")) +
    
    (metadata %>% 
       ggplot(aes(x = nFeature_RNA, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       scale_x_log10() +
       theme_classic() +
       ylab("Cell density") +
       geom_vline(xintercept = 1200) +
       theme_custom +
       ggtitle("Number of detected genes/cell") +
       guides(fill = "none")) +
    
    (metadata %>% 
       ggplot(aes(x = nCount_RNA, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       scale_x_log10() +
       theme_classic() +
       ylab("Cell density") +
       geom_vline(xintercept = 2500) +
       theme_custom +
       ggtitle("RNA molecules detected/cell") +
       guides(fill = "none")) +
    
    (metadata %>% 
       filter(mitoRatio > 0) %>%
       ggplot(aes(x = mitoRatio, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_vline(xintercept = 0.075) +
       theme_custom +
       ggtitle("Reads mapped to MT genes") +
       guides(fill = "none")) +
    
    (metadata %>% 
       filter(riboRatio > 0) %>%
       ggplot(aes(x = riboRatio, fill = sample)) + 
       geom_density(alpha = 0.2) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_vline(xintercept = 0.1) +
       ylab("Cell density") +
       xlab("riboRatio") +
       theme_custom +
       ggtitle("Reads mapped to ribosomal genes") +
       guides(fill = "none")) +
    
    (metadata %>% 
       ggplot(aes(x = sample, y = log10GenesPerUMI, fill = sample)) + 
       geom_violin() +
       geom_boxplot(width = 0.1, fill = alpha("white", 0.7)) +
       scale_fill_manual(values = custom_colors) +
       theme_classic() +
       geom_hline(yintercept = 0.8) +
       ylab("log10GenesPerUMI") +
       xlab(NULL) +
       theme_custom + 
       ggtitle("Library complexity") +
       guides(fill = "none"))
)
print(combined_plot)

# Save to PNG
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
png("filtereddata_ST101_only_all.png", width = 13, height = 6, units = 'in', res = 300)
print(combined_plot)
dev.off()


# --------------------------
# Correlation heatmap workflow
# --------------------------
#load("step3A_corrplot_ST101_analysis.RData")

# Normalize
library(future)
options(future.globals.maxSize = 10 * 1024^3)

tic()
seuset <- NormalizeData(seuset, verbose = FALSE)
toc()

library(svglite)
library(dplyr)
# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
avg_exp <- AggregateExpression(seuset, group.by="orig.ident", assays = "RNA") %>% as.data.frame()
colnames(avg_exp) <- gsub(pattern = "RNA.", "", colnames(avg_exp))
cormat <- round(cor(avg_exp, use="pairwise.complete.obs"), 2)
sample_order <- rownames(cormat)
lower_tri <- get_lower_tri(cormat)
melted_cormat <- reshape2::melt(lower_tri, na.rm = TRUE)


ggheatmap <- melted_cormat %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Var1 = factor(Var1, levels = sample_order)) %>%
  mutate(Var2 = factor(Var2, levels = sample_order)) %>%
  ggplot(aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4BC9FF",
    high = "#D50032",
    mid = "white",
    midpoint = 0.5, limit = c(0, 1), space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed() +
  geom_text(aes(label = value), color = "black", size = 4, family = "Arial", fontface = "bold") +
  guides(
    fill = guide_colorbar(
      barwidth = 7, barheight = 1,
      title.position = "top", title.hjust = 0.5
    )
  )
ggheatmap
png("step3_corrPlot.png", width = 6, height = 6, units = "in", res = 300)
print(ggheatmap)
dev.off()

save.image("step2_filtered_normalized_ST-101_res1.RData")

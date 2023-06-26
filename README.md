# Single-cell-RNA-seq

REF: https://github.com/SomenMistri/intro_to_scRNA-seq/tree/main

## Table of contents    
* [Page 1: 2023-25-06](#id-section1). Chapter 1: Data loading

* [Page 2: 2023-30-06](#id-section2). Chapter 2: Visualize the common QC metrics

* [Page 3: 2023-30-06](#id-section3). Chapter 3: Data integration

* [Page 4: 2023-30-06](#id-section3). Chapter 4: Clustering

------
<div id='id-section1'/>

## Chapter 1: Data loading

LOAD NECESSARY R libraries 

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("XVector",force = TRUE)
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")

install.packages('tidyverse')
install.packages('Matrix')
install.packages('RCurl')
install.packages('scales')
install.packages('metap')
install.packages('Seurat')
install.packages("ggplot2")
install.packages("sctransform")

library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```
### Load individual count matrices

The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).

```
data1 <- Read10X("PDAC_tissue_1_filtered_feature_bc_matrix")
data2 <- Read10X("PDAC_tissue_2_filtered_feature_bc_matrix")
data3 <- Read10X("PDAC_tissue_3_filtered_feature_bc_matrix")
data4 <- Read10X("PDAC_tissue_4_filtered_feature_bc_matrix")
```
### Create Seurat Objects
Let's use the individual count matrices to create separate Seurat objects. The seurat object serves as a container for both the data (like the count matrix) and analysis (e.g. PCA, metadata) for a single-cell dataset.

```
data_seurat1 <- CreateSeuratObject(counts = data1, project = "Human-1", min.cells = 3, min.features = 200)
data_seurat2 <- CreateSeuratObject(counts = data2, project = "Human-2", min.cells = 3, min.features = 200)
data_seurat3 <- CreateSeuratObject(counts = data3, project = "Human-3", min.cells = 3, min.features = 200)
data_seurat4 <- CreateSeuratObject(counts = data4, project = "Human-4", min.cells = 3, min.features = 200)
```

### Cell Cycle Scoring

In certain scenarios, it is necessary to address the impact of cell cycle heterogeneity in scRNA-seq data. One approach to mitigate this effect involves calculating cell cycle phase scores using known cell cycle markers and removing their influence during the pre-processing stage.

To accomplish this, the following steps are performed. First, individual Seurat objects are log normalized using the NormalizeData() function. Subsequently, the CellCycleScoring() function is utilized to assign each cell a cell cycle score based on its expression of markers associated with the G2/M and S phases. Seurat stores the S-phase genes and G2/M-phase genes in the "cc.genes.updated.2019" list.

```
#segregate the "cc.genes.updated.2019" list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Prior to running "CellCycleScoring" command, each seurat object needs to be Lognormalized using "NormalizeData" function
data_norm1 <- NormalizeData(data_seurat1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm2 <- NormalizeData(data_seurat2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm3 <- NormalizeData(data_seurat3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm4 <- NormalizeData(data_seurat4, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#Now perform CellCycleScoring for each seurat objects
data_norm1 <- CellCycleScoring(data_norm1, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm2 <- CellCycleScoring(data_norm2, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm3 <- CellCycleScoring(data_norm3, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm4 <- CellCycleScoring(data_norm4, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)

# view cell cycle scores and phase assignments
head(data_norm1[[]])
```
**nFeature_RNA** is the **number of genes detected in each cell**. Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet.
**nCount_RNA** is the **total number of molecules detected within a cell**. High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet). 

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/nCount_RNA_nFeature_RNA.png)

### Merge individual seurat objects into one
Merge all four seurat objects into one. The merge() function merges the raw count matrices of two or more Seurat objects creating a new Seurat object with a combined raw count matrix. Then, let's take a look at the metadata of the merged seurat object using the View() function.

```
#NOTE: By default, merge() function combines Seurat objects based on the raw count matrices, erasing any previous normalization
data_merged <- merge(data_norm1, y = c(data_norm2, data_norm3, data_norm4), add.cell.ids = c("H1", "H2", "H3","H4"), project = "Human_1234")
```

To make sure that cells from all the human samples were merged properly, you can use the table() function:

```
table(data_merged$orig.ident)
```

### Calculate additional quality control metrics
To determine the **mitochondrial and ribosomal transcript percentage per cell**, Seurat provides a convenient function called PercentageFeatureSet(). This function allows us to search for specific patterns within the dataset. In this case, we can utilize it to identify mitochondrial and ribosomal genes. To calculate the mitochondrial transcript percentage, we can use the pattern "MT-" to search for genes associated with mitochondria. Cells with a high proportion of mitochondrial genes are typically considered low quality. For ribosomal transcript percentages, we can search for genes using the pattern "^RP[SL]". It's important to note that the percentage of ribosomal transcript varies significantly across different cell types. Therefore, caution should be exercised when using percent.RIBO values to filter out low-quality cells, as they can exhibit substantial variation depending on the cell type.

```
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#First add column with mitochondiral gene expression
data_merged[["percent.MT"]] <- PercentageFeatureSet(data_merged, pattern = "^MT-")
#Add column with ribosomal gene expression
data_merged[["percent.RIBO"]] <- PercentageFeatureSet(data_merged, pattern = "^RP[SL]")
#NOTE: this calculation is performed per cell. That is why this step can be performed on merged data

head (data_merged)

                      orig.ident nCount_RNA nFeature_RNA      S.Score   G2M.Score Phase percent.MT percent.RIBO
H1_AAACGAAAGTGGAAAG-1    Human-1        501          269 -0.021500501  0.02518923   G2M  12.574850    0.3992016
H1_AAACGAAGTAGGGTAC-1    Human-1      26196         4641 -0.043996838 -0.08059179    G1  45.132845    3.8402810
H1_AAACGAAGTCATAGTC-1    Human-1      72725         7990 -0.055117973 -0.09913018    G1  14.349948   14.3045720
H1_AAAGAACCATTAAAGG-1    Human-1      12734         4046 -0.027648932 -0.08673782    G1   3.683053    8.0021988
```
-----
<div id='id-section2'/>

## Chapter 2: Visualize the common QC metrics

### Before filering
Now plot the common QC metrics to see-
* #### UMI counts per cell (nCount_RNA)
* #### Genes detected per cell (nFeature_RNA_)
* #### Mitochondrial counts ratio
* #### Ribosomal counts ratio
```
VlnPlot(data_merged, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(data_merged, features = c("percent.MT","percent.RIBO"), ncol = 2)
FeatureScatter(data_merged, feature1 = "percent.RIBO", feature2 = "percent.MT")
FeatureScatter(data_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/before_filtering.png)

#### After filtering

##### Doublets
Doublets pose a challenge in single-cell RNA sequencing experiments as they result from the unintended combination of two cells. These occurrences commonly happen during cell sorting or capture, particularly in droplet-based protocols that handle a large number of cells. The presence of doublets is undesirable when the goal is to accurately characterize populations at the single-cell level. It can lead to erroneous indications of intermediate populations or transient states that do not genuinely exist. Consequently, it is crucial to eliminate doublet libraries to ensure the integrity of the interpretation derived from the results.

It is worth noting that many existing workflows rely on setting maximum thresholds for unique molecular identifiers (UMIs) or genes as an indicator of multiple cells. Although this approach may seem intuitive, it is not always reliable. Moreover, several tools used to detect doublets have a tendency to discard cells with intermediate or continuous phenotypes, which can be problematic when dealing with datasets containing cell types that exhibit more nuanced characteristics.

Good cells will generally exhibit both higher number of genes per cell (nFeature_RNA) and higher numbers of UMIs (nCount_RNA) per cell. Cells that are poor quality are likely to have low nFeature_RNA and nCount_RNA. Also Mitochondrial read fractions are only high in particularly low count cells with few detected genes.

Now that we have visualized the various metrics, we can decide on the thresholds to apply which will result in the removal of low quality cells. Often the recommendations mentioned earlier are a rough guideline, and the specific experiment needs to inform the exact thresholds chosen. We will use the following thresholds:

* nFeature_RNA > 500
* nCount_RNA < 700
* percent.MT < 25
* percent.RIBO > 3

```
data_filtered <- subset(data_merged, subset = nFeature_RNA > 500 & nCount_RNA < 200000 & percent.MT < 25 & percent.RIBO > 3)
VlnPlot(data_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(data_filtered, features = c("percent.MT","percent.RIBO"), ncol = 2)
FeatureScatter(data_filtered, feature1 = "percent.RIBO", feature2 = "percent.MT")
FeatureScatter(data_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/after_filtering.png)

REF: Jihe Liu, William J. Gammerdinger, Meeta Mistry, Mary E. Piper, & Radhika S. Khetani. (2022, January 6). hbctraining/Intro-to-shell-flipped: Shell and HPC Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.5826091

#### Save as .rds

```
saveRDS(data_filtered, file = "data_filtered.rds")
```

-----
<div id='id-section3'/>

## Chapter 3: Data integration

### Normalization
Normalization is an important initial step in analyzing mRNA expression data. It involves adjusting the expression counts of genes to account for systematic variations, making them comparable across different genes and samples. When working with single-cell RNA sequencing (scRNA-seq) data, specific normalization methods are employed.

### Methods
Some scRNA-seq normalization methods resemble those used in bulk RNA-seq analysis, where global scale factors are applied to account for a common count-depth relationship across genes. However, these simple methods can lead to over-correction for low and moderately expressed genes and under-normalization for highly expressed genes if the assumptions underlying them are not valid. More complex normalization methods take into consideration the characteristics of each individual gene and apply corrections accordingly.

Simple transformations involve applying the same mathematical function to each measurement. Common examples include **log transforms (used in the original Seurat workflow) or square root transforms (less frequently used)**. However, a study by Hafemeister and Satija in 2019 highlighted issues with simple transformations. They found that the standard log normalization approach affects genes with different abundances in varying ways, and effective normalization using log transforms is only observed with low to medium abundance genes. Additionally, significant imbalances in variance were observed when working with log-normalized data, indicating that all genes cannot be treated equally during normalization.

#### Pearson residuals for transformation
The suggested solution to address the issues with simple transformations is the utilization of **Pearson residuals for transformation**, which is implemented in Seurat's SCTransform function. This approach involves the following steps:

* Multiplication by gene-specific weights: The measurements of each gene are multiplied by a weight that is specific to that gene.

* Weighting based on non-uniform expression: Each gene is assigned a weight based on the amount of evidence indicating that it is expressed non-uniformly across cells. Genes that exhibit expression in only a small fraction of cells are given higher weights. This is particularly useful for identifying rare cell populations.

* Consideration of both expression level and distribution: The approach not only takes into account the expression level of genes but also considers the distribution of their expression across cells.

By incorporating these steps, the Pearson residuals transformation provides a more nuanced and effective normalization method for scRNA-seq data, addressing the limitations observed with simple transformations.

```
data.filtered <- readRDS ("data_filtered.rds")
```
During this lesson, our main approach for data normalization will be to utilize the **"SCTransform"** function within Seurat. It's important to note that this single command, SCTransform(), replaces the previous steps of NormalizeData(), ScaleData(), and FindVariableFeatures() in the original Seurat workflow, which involved log-normalization. By using SCTransform(), we can perform normalization in a more efficient and streamlined manner.

```
data_SCT <- SCTransform(data.filtered, verbose = TRUE)
```

### Clustering by Principal Component Analysis (PCA)

In order to handle the technical noise in the gene expression data obtained from scRNA-seq, Seurat uses PCA scores derived from the expression of the most variable genes to assign cells to clusters. Each principal component (PC) represents a combination of genes that are correlated, acting as a "metagene." Determining the appropriate number of PCs to include in the clustering process is crucial to ensure that we capture most of the variation and cell types present in the dataset.

Before deciding which PCs to include for clustering, it is beneficial to explore them. One approach is to create a heatmap that visualizes the most variant genes for selected PCs, arranging the genes and cells based on their PCA scores. The purpose is to examine the PCs and assess whether the genes driving them make sense for distinguishing different cell types.

Another useful method for determining the number of PCs to use in clustering is the elbow plot. This plot illustrates the standard deviation of each PC, and we look for the point where the standard deviations start to level off. Typically, the elbow indicates the threshold for capturing the majority of the variation. However, it's important to note that this approach can be subjective.

```
data_PCA<- RunPCA(data_SCT, npcs = 40, verbose = FALSE)

# Explore heatmap of PCs
DimHeatmap(data_PCA, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)

# Let's draw the elbow plot using the top 40 PCs
ElbowPlot(object = data_PCA, 
          ndims = 40)
```

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/joined_pca_plots.png)



-----
<div id='id-section4'/>

## Chapter 4: Clustering

Clusters of cells are formed by grouping cells together based on the **similarity of their gene expression profiles**. To determine expression profile similarity, distance metrics are used, often relying on dimensionality-reduced representations as input. Seurat assigns cells to clusters by considering their PCA scores derived from the expression of the most variable genes.

Although PCA determines all principal components (PCs), we can only visualize two at a time. On the other hand, dimensionality reduction techniques like **Uniform Manifold Approximation and Projection (UMAP)** utilize **information from any number of top PCs to position cells in a multidimensional space**. UMAP calculates distances in this multidimensional space and then plots them in two dimensions, aiming to **preserve both local and global structure**. Thus, the distances between cells in the plot represent their similarity in gene expression.

For cell clustering, Seurat employs a graph-based approach using a K-nearest neighbor method. It attempts to partition the graph into densely interconnected "quasi-cliques" or "communities." The initial step involves constructing a **K-nearest neighbor (KNN) graph** based on the Euclidean distance in PCA space. Seurat accomplishes this through the FindNeighbors() function. Subsequently, cells are grouped iteratively in order to optimize the standard modularity function. The FindClusters() function in Seurat handles the graph-based clustering. The "resolution" parameter plays a crucial role at this stage, as it determines the level of detail in the resulting clusters and needs to be optimized for each experiment. In general, for **datasets containing 3,000 to 5,000 cells, a resolution between 0.4 and 1.4 tends to yield good clustering**. Larger datasets often require higher resolution values to generate a greater number of clusters.

Note: When running the RunUMAP(), FindNeighbors(), and FindClusters() functions sequentially (as shown in chunk 6), it is advisable to use the same number of PCA dimensions as input for both RunUMAP() and FindNeighbors().

```
# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
data_clust <- RunUMAP(data_PCA, reduction = "pca", dims = 1:30)
# Determine the K-nearest neighbor graph
data_clust <- FindNeighbors(object = data_clust, 
                                dims = 1:30)
# Perform graph based clustering
data_clust <- FindClusters(object = data_clust,
                               resolution = 0.6)

# Visualize clustered cells
DimPlot(data_clust, reduction = "umap", label = TRUE) + NoLegend()

# Save the clustered plot
ggsave(path = "Figs", filename = "Clusters.png",  height=5, width=6, units='in', dpi = 300, bg = "transparent", device='png')
```
![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/Clusters.png)

### Undesired variations 
To begin the workflow, it is important to assess whether our data exhibits any undesired variations. In single-cell RNA-seq data, the most frequently examined biological effect is the **influence of the cell cycle on the transcriptome**. Additionally, there may be technical sources of variation, such as **batch effects**. This step of the workflow entails thoroughly examining the data to determine the specific covariates that should be accounted for and corrected.

#### Explore the effects of cell cycle genes:
```
head(data_clust@meta.data)

# Evaluating effects of cell cycle (Phase)
DimPlot(data_clust, group.by = "Phase", label = FALSE)

# Save the plot
ggsave(path = "Figs", filename = "CellCycle_Phase.png",  height=5, width=7, units='in', dpi = 300, bg = "transparent", device='png')
```

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/CellCycle_Phase.png)

#### Now, explore technical sources of variation such as the Batch Effect:
# Set identity classes to seurat_clusters
Idents(object = data_clust) <- "seurat_clusters"

# Explore the significance of back effect on clustering
DimPlot(data_clust, split.by = "orig.ident", label = TRUE, ncol = 2)+ NoLegend()

# Save the plot
ggsave(path = "Figs", filename = "Batch_effect.png",  height=6, width=8, units='in', dpi = 300, bg = "transparent", device='png')

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/Batch_effect.png)

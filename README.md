# Single-cell-RNA-seq

REF: https://github.com/SomenMistri/intro_to_scRNA-seq/tree/main

## Table of contents    
* [Page 1: 2023-25-06](#id-section1). Chapter 1: Data loading

* [Page 2: 2023-30-06](#id-section2). Chapter 2: Visualize the common QC metrics

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

In some cases, there is a need for mitigating the effects of cell cycle heterogeneity in scRNA-seq data.This can be done by calculating cell cycle phase scores based on known cell cycle markers , and regressing these out of the data during pre-processing.

Here we first Log Normalizing individual seurat objects using the NormalizeData() function. Then, we are using the CellCycleScoring() function to assign each cell a cell cycle score, based on its expression of G2/M and S phase markers. Seurat stores the s.genes and g2m.genes in the "cc.genes.updated.2019" list.


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

Calculate the **mitochondrial and ribosomal transcript percentage per cell**. Seurat has a function that enables us to do this. The PercentageFeatureSet() function can take a specific pattern and search through the dataset for that pattern. We can search for mitochondrial genes by looking for the pattern "MT-". Similarly, for the ribosomal genes, we can look for the pattern "^RP[SL]". Usually, cells with high proportions of mitochondrial genes are considered as poor-quality cells. On the other hand, percentage of ribosomal transcript per cell varies greatly from cell type to cell type. Therefore, caution should be taken to use percent.RIBO values to filter out low quality cells.

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

**nFeature_RNA > 500
nCount_RNA < 700
percent.MT < 25
percent.RIBO > 3**

```
data_filtered <- subset(data_merged, subset = nFeature_RNA > 500 & nCount_RNA < 200000 & percent.MT < 25 & percent.RIBO > 3)
VlnPlot(data_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(data_filtered, features = c("percent.MT","percent.RIBO"), ncol = 2)
FeatureScatter(data_filtered, feature1 = "percent.RIBO", feature2 = "percent.MT")
FeatureScatter(data_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

![alt text](https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/after_filtering.png)

REF: Jihe Liu, William J. Gammerdinger, Meeta Mistry, Mary E. Piper, & Radhika S. Khetani. (2022, January 6). hbctraining/Intro-to-shell-flipped: Shell and HPC Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.5826091






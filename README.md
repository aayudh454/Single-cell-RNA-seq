# Single-cell-RNA-seq

REF: https://github.com/SomenMistri/intro_to_scRNA-seq/tree/main

## Table of contents    
* [Page 1: 2023-25-06](#id-section1). Chapter 1: Data loading

* [Page 2: 2023-30-06](#id-section2). Chapter 2: Processing sc-RNA data



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

![alt text]([(https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/nCount_RNA_nFeature_RNA.png)]

https://github.com/aayudh454/Single-cell-RNA-seq/blob/main/nCount_RNA_nFeature_RNA.png

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

## Chapter 2: Processing sc-RNA data

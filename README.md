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


-----
<div id='id-section2'/>

## Chapter 2: Processing sc-RNA data

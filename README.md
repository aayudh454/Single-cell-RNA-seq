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


-----
<div id='id-section2'/>

## Chapter 2: Processing sc-RNA data

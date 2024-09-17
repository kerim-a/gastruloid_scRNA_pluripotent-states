# load libraries (not all are required)

library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scCustomize)
library(dittoSeq)
library(SeuratWrappers)
library(RColorBrewer)
library(rcartocolor)

#This script describes how scRNA datasets (CellRanger output) were loaded and processed, i.e. initial QC, prior to merging etc.


### 0h / mESC datasets


#0h serum+LIF (base) 1
#Load data
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022\\ESL_0h\\filtered_feature_bc_matrix')
ESL_0h_old <- CreateSeuratObject(counts = matrix_0h.data, project = "ESL0h_old", min.cells = 3, min.features = 200)
ESL_0h_old[["percent.mt"]] <- PercentageFeatureSet(ESL_0h_old, pattern = "^mt-")
VlnPlot(ESL_0h_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC filtering
ESL_0h_filtered_old <- subset(ESL_0h_old, subset = nFeature_RNA > 2750 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(ESL_0h_filtered_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#For later identification
ESL_0h_filtered_old$merge.ident <- "ESL0h"
ESL_0h_filtered_old$fuse.ident <- "ESL0h_old"
#Data Normalization
ESL_0h_filtered_old <- NormalizeData(ESL_0h_filtered_old, normalization.method = "LogNormalize", scale.factor = 10000)


#0h serum+LIF (base) 2
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_10\\ESL_0h\\filtered_feature_bc_matrix')
ESL_0h <- CreateSeuratObject(counts = matrix_0h.data, project = "ESL_0h_new", min.cells = 3, min.features = 200)
ESL_0h[["percent.mt"]] <- PercentageFeatureSet(ESL_0h, pattern = "^mt-")
VlnPlot(ESL_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_new <- subset(ESL_0h, subset = nFeature_RNA > 2600 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(ESL_0h_filtered_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_new$merge.ident <- "ESL0h"
ESL_0h_filtered_new$fuse.ident <- "ESL0h_new"

ESL_0h_filtered_new <- NormalizeData(ESL_0h_filtered_new, normalization.method = "LogNormalize", scale.factor = 10000)


#0h serum+LIF+PD03 (base+1i) 
matrix_PDnew0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\09_S02_PD0h\\filtered_feature_bc_matrix')
PDnew_0h <- CreateSeuratObject(counts = matrix_PDnew0h.data, project = "PD0h_3", min.cells = 3, min.features = 200)
PDnew_0h[["percent.mt"]] <- PercentageFeatureSet(PDnew_0h, pattern = "^mt-")
VlnPlot(PDnew_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_0h_filtered_newest <- subset(PDnew_0h, subset = nFeature_RNA > 2750 & percent.mt < 15 & percent.mt > 0.5 & nCount_RNA < 130000)
VlnPlot(PD_0h_filtered_newest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_0h_filtered_newest$merge.ident <- "PD0h"
PD_0h_filtered_newest$fuse.ident <- "PD0h_new"

PD_0h_filtered_newest <- NormalizeData(PD_0h_filtered_newest, normalization.method = "LogNormalize", scale.factor = 10000)


#0h serum+LIF+CGP77+CHIR99 (base+2i)
matrix_a2i0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\08_S01_a2i0h\\filtered_feature_bc_matrix')
a2i_0h <- CreateSeuratObject(counts = matrix_a2i0h.data, project = "a2i0h", min.cells = 3, min.features = 200)
a2i_0h[["percent.mt"]] <- PercentageFeatureSet(a2i_0h, pattern = "^mt-")
VlnPlot(a2i_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_0h_filtered <- subset(a2i_0h, subset = nFeature_RNA > 2500 & percent.mt < 15 & percent.mt > 0.5 & nCount_RNA < 130000)
VlnPlot(a2i_0h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_0h_filtered$merge.ident <- "a2i0h"
a2i_0h_filtered$fuse.ident <- "a2i0h_1"

a2i_0h_filtered <- NormalizeData(a2i_0h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


### Pre-CHIR99 datasets


#48hpa serum+LIF (base) 

matrix_48h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2021_newp_esl\\just_matrices\\48h\\filtered_feature_bc_matrix')
ESL_48h <- CreateSeuratObject(counts = matrix_48h.data, project = "ESL48h", min.cells = 3, min.features = 200)
ESL_48h[["percent.mt"]] <- PercentageFeatureSet(ESL_48h, pattern = "^mt-")
VlnPlot(ESL_48h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_48h_filtered <- subset(ESL_48h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 100000)
VlnPlot(ESL_48h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_48h_filtered$merge.ident <- "ESL48h"
ESL_48h_filtered$fuse.ident <- "ESL48h_1"

ESL_48h_filtered <- NormalizeData(ESL_48h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#48hpa serum+LIF+PD03 (base+1i) 1

matrix_48h.data2 <- Read10X(data.dir = 'Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\2021_data\\matrices\\20210911_CellRangerCount_VT10x03_S02_PD48H\\filtered_feature_bc_matrix')
g_48h <- CreateSeuratObject(counts = matrix_48h.data2, project = "pd48h_3", min.cells = 3, min.features = 200)
g_48h[["percent.mt"]] <- PercentageFeatureSet(g_48h, pattern = "^mt-")
VlnPlot(g_48h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_48h_filtered3 <- subset(g_48h, subset = nFeature_RNA > 3000 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 250000)
VlnPlot(PD_48h_filtered3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_48h_filtered3$fuse.ident <- "PD48h_2"
PD_48h_filtered3$merge.ident <- "PD48h"

PD_48h_filtered3 <- NormalizeData(PD_48h_filtered3, normalization.method = "LogNormalize", scale.factor = 10000)


#48hpa serum+LIF+PD03 (base+1i) 2

matrix_PD48h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_10\\PD_48h\\filtered_feature_bc_matrix')
PD_48h <- CreateSeuratObject(counts = matrix_PD48h.data, project = "PD_48h", min.cells = 3, min.features = 200)
PD_48h[["percent.mt"]] <- PercentageFeatureSet(PD_48h, pattern = "^mt-")
VlnPlot(PD_48h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_48h_filtered <- subset(PD_48h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 130000)
VlnPlot(PD_48h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_48h_filtered$merge.ident <- "PD48h"
PD_48h_filtered$fuse.ident <- "PD48h_3"

PD_48h_filtered <- NormalizeData(PD_48h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#72hpa serum+LIF+CGP77+CHIR99 (base+2i)

matrix_72h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\08_S02_a2i72h\\filtered_feature_bc_matrix')
a2i_72h <- CreateSeuratObject(counts = matrix_72h.data, project = "a2i_72h", min.cells = 3, min.features = 200)
a2i_72h[["percent.mt"]] <- PercentageFeatureSet(a2i_72h, pattern = "^mt-")
VlnPlot(a2i_72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_72h_filtered <- subset(a2i_72h, subset = nFeature_RNA > 2000 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(a2i_72h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_72h_filtered <- NormalizeData(a2i_72h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

a2i_72h_filtered$merge.ident <- "a2i_72h"
a2i_72h_filtered$fuse.ident <- "a2i_72h_1"


### Post-CHIR99 datasets


#72hpa serum+LIF (base) 1

matrix_72h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022\\ESL_72h\\filtered_feature_bc_matrix')
ESL_72h <- CreateSeuratObject(counts = matrix_72h.data, project = "ESL_72h", min.cells = 3, min.features = 200)
ESL_72h[["percent.mt"]] <- PercentageFeatureSet(ESL_72h, pattern = "^mt-")
VlnPlot(ESL_72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_72h_filtered <- subset(ESL_72h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 125000)
VlnPlot(ESL_72h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_72h_filtered$merge.ident <- "ESL_72h_new"
ESL_72h_filtered$fuse.ident <- "ESL_72h_new_1"

ESL_72h_filtered <- NormalizeData(ESL_72h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#72hpa serum+LIF (base) 2

matrix_72h.data_2 <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2021_newp_esl\\just_matrices\\72h\\filtered_feature_bc_matrix')
ESL_72h_2 <- CreateSeuratObject(counts = matrix_72h.data_2, project = "ESL72h", min.cells = 3, min.features = 200)
ESL_72h_2[["percent.mt"]] <- PercentageFeatureSet(ESL_72h_2, pattern = "^mt-")
VlnPlot(ESL_72h_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_72h_filtered_2 <- subset(ESL_72h_2, subset = nFeature_RNA > 3500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 300000)
VlnPlot(ESL_72h_filtered_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_72h_filtered_2$merge.ident <- "ESL72h"
ESL_72h_filtered_2$fuse.ident <- "ESL72h_2"

ESL_72h_filtered_2 <- NormalizeData(ESL_72h_filtered_2, normalization.method = "LogNormalize", scale.factor = 10000)


#72hpa serum+LIF+PD03 (base+1i) 

matrix_PD72h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\09_S04_PD72h\\filtered_feature_bc_matrix')
PD_72h <- CreateSeuratObject(counts = matrix_PD72h.data, project = "PD_72h", min.cells = 3, min.features = 200)
PD_72h[["percent.mt"]] <- PercentageFeatureSet(PD_72h, pattern = "^mt-")
VlnPlot(PD_72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_72h_filtered <- subset(PD_72h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 170000)
VlnPlot(PD_72h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_72h_filtered$merge.ident <- "PD72h"
PD_72h_filtered$fuse.ident <- "PD72h_1"

PD_72h_filtered <- NormalizeData(PD_72h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#96hpa serum+LIF+CGP77+CHIR99 (base+2i)

matrix_96h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\08_S03_a2i96h\\filtered_feature_bc_matrix')
a2i_96h <- CreateSeuratObject(counts = matrix_96h.data, project = "a2i_96h", min.cells = 3, min.features = 200)
a2i_96h[["percent.mt"]] <- PercentageFeatureSet(a2i_96h, pattern = "^mt-")
VlnPlot(a2i_96h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_96h_filtered <- subset(a2i_96h, subset = nFeature_RNA > 2200 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(a2i_96h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_96h_filtered$merge.ident <- "a2i_96h"
a2i_96h_filtered$fuse.ident <- "a2i_96h_1"

a2i_96h_filtered <- NormalizeData(a2i_96h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


### Max. Elongation Datasets


#120hpa serum+LIF (base)

matrix_120h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2021_newp_esl\\just_matrices\\120h\\filtered_feature_bc_matrix')
ESL_120h <- CreateSeuratObject(counts = matrix_120h.data, project = "ESL120h", min.cells = 3, min.features = 200)
ESL_120h[["percent.mt"]] <- PercentageFeatureSet(ESL_120h, pattern = "^mt-")
VlnPlot(ESL_120h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_120h_filtered <- subset(ESL_120h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 100000)
VlnPlot(ESL_120h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_120h_filtered$merge.ident <- "ESL120h"
ESL_120h_filtered$fuse.ident <- "ESL120h_1"

ESL_120h_filtered <- NormalizeData(ESL_120h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#120hpa serum+LIF+PD03 (base+1i) 

matrix_120h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\07_S04_PD120h\\filtered_feature_bc_matrix')
PD_120h <- CreateSeuratObject(counts = matrix_120h.data, project = "PD_120h", min.cells = 3, min.features = 200)
PD_120h[["percent.mt"]] <- PercentageFeatureSet(PD_120h, pattern = "^mt-")
VlnPlot(PD_120h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_120h_filtered <- subset(PD_120h, subset = nFeature_RNA > 2400 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(PD_120h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_120h_filtered$merge.ident <- "PD120h"
PD_120h_filtered$fuse.ident <- "PD120h_1"

PD_120h_filtered <- NormalizeData(PD_120h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)


#144hpa serum+LIF+CGP77+CHIR99 (base+2i)

matrix_144h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\08_S04_a2i144h\\filtered_feature_bc_matrix')
a2i_144h <- CreateSeuratObject(counts = matrix_144h.data, project = "a2i_144h", min.cells = 3, min.features = 200)
a2i_144h[["percent.mt"]] <- PercentageFeatureSet(a2i_144h, pattern = "^mt-")
VlnPlot(a2i_144h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_144h_filtered <- subset(a2i_144h, subset = nFeature_RNA > 2250 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(a2i_144h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_144h_filtered$merge.ident <- "a2i_144h"
a2i_144h_filtered$fuse.ident <- "a2i_144h_1"

a2i_144h_filtered <- NormalizeData(a2i_144h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
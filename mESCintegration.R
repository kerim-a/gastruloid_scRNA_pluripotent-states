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

#Load and prepare individual datasets to be integrated

#0h serum+LIF (base) 1
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022\\ESL_0h\\filtered_feature_bc_matrix')
ESL_0h_old <- CreateSeuratObject(counts = matrix_0h.data, project = "ESL0h_old", min.cells = 3, min.features = 200)
ESL_0h_old[["percent.mt"]] <- PercentageFeatureSet(ESL_0h_old, pattern = "^mt-")
VlnPlot(ESL_0h_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_old <- subset(ESL_0h_old, subset = nFeature_RNA > 2750 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(ESL_0h_filtered_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_old$merge.ident <- "ESL0h"
ESL_0h_filtered_old$fuse.ident <- "ESL0h_old"

#0h serum+LIF (base) 2
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_10\\ESL_0h\\filtered_feature_bc_matrix')
ESL_0h <- CreateSeuratObject(counts = matrix_0h.data, project = "ESL_0h_new", min.cells = 3, min.features = 200)
ESL_0h[["percent.mt"]] <- PercentageFeatureSet(ESL_0h, pattern = "^mt-")
VlnPlot(ESL_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_new <- subset(ESL_0h, subset = nFeature_RNA > 2600 & percent.mt < 20 & percent.mt > 0.5 & nCount_RNA < 150000)
VlnPlot(ESL_0h_filtered_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ESL_0h_filtered_new$merge.ident <- "ESL0h"
ESL_0h_filtered_new$fuse.ident <- "ESL0h_new"

#0h serum+LIF+PD03 (base+1i) 
matrix_PDnew0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\09_S02_PD0h\\filtered_feature_bc_matrix')
PDnew_0h <- CreateSeuratObject(counts = matrix_PDnew0h.data, project = "PD0h_3", min.cells = 3, min.features = 200)
PDnew_0h[["percent.mt"]] <- PercentageFeatureSet(PDnew_0h, pattern = "^mt-")
VlnPlot(PDnew_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_0h_filtered_newest <- subset(PDnew_0h, subset = nFeature_RNA > 2750 & percent.mt < 15 & percent.mt > 0.5 & nCount_RNA < 130000)
VlnPlot(PD_0h_filtered_newest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PD_0h_filtered_newest$merge.ident <- "PD0h"
PD_0h_filtered_newest$fuse.ident <- "PD0h_new"

#0h serum+LIF+CGP77+CHIR99 (base+2i)
matrix_a2i0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2022_11_morereads\\08_S01_a2i0h\\filtered_feature_bc_matrix')
a2i_0h <- CreateSeuratObject(counts = matrix_a2i0h.data, project = "a2i0h", min.cells = 3, min.features = 200)
a2i_0h[["percent.mt"]] <- PercentageFeatureSet(a2i_0h, pattern = "^mt-")
VlnPlot(a2i_0h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_0h_filtered <- subset(a2i_0h, subset = nFeature_RNA > 2500 & percent.mt < 15 & percent.mt > 0.5 & nCount_RNA < 130000)
VlnPlot(a2i_0h_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

a2i_0h_filtered$merge.ident <- "a2i0h"
a2i_0h_filtered$fuse.ident <- "a2i0h_1"

#Set the assay to RNA (should be default)
DefaultAssay(ESL_0h_filtered_new) <- "RNA"
DefaultAssay(ESL_0h_filtered_old) <- "RNA"
DefaultAssay(PD_0h_filtered_newest) <- "RNA"
DefaultAssay(a2i_0h_filtered) <- "RNA"

#Merge into single object, remove the scale.data slot, normalize data and compute MVGs altogether
mESC_merged <- merge(ESL_0h_filtered_new, y = c(ESL_0h_filtered_old,PD_0h_filtered_newest,a2i_0h_filtered), add.cell.ids = c('ESL0h_new','ESL0h_old','PD0h','a2i0h'), project = "mESCmerged")
mESC_merged@assays$RNA@scale.data <- as.matrix(0)
mESC_merged <- NormalizeData(mESC_merged)
mESC_merged <- FindVariableFeatures(mESC_merged, selection.method = "vst", nfeatures = 3000)

#Select MVGs for integration and remove cell cycle genes from Seurats built-in list (mouse orthologues)
features <- SelectIntegrationFeatures(object.list = SplitObject(mESC_merged, split.by = "orig.ident"), nfeatures = 3000)
features_filtered1 <- setdiff(features, g2m_genes)
features_filtered2 <- setdiff(features_filtered1, s_genes)

#Run the integration
mESC_merged <- RunFastMNN(object.list = SplitObject(mESC_merged, split.by = "orig.ident"), features = features_filtered2)

#UMAP and clustering
mESC_merged <- RunUMAP(mESC_merged, reduction = "mnn", dims = 1:25)
mESC_merged <- FindNeighbors(mESC_merged, reduction = "mnn", dims = 1:25)
mESC_merged <- FindClusters(mESC_merged, resolution = 0.3)

#Visualization of UMAP
p <- DimPlot(mESC_merged, reduction = "umap")
AugmentPlot(plot = p)
plot(p)

#Evaluate marker genes and save
DefaultAssay(mESC_merged) <- "RNA"
mESC_merged_markers <- FindAllMarkers(mESC_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.35)

#Clusters 0 and 1 appear very similar, so they can reasonably be merged
new.cluster.ids <- c("0", "0", "2", "3", "4", "5")
names(new.cluster.ids) <- levels(mESC_merged)
mESC_merged <- RenameIdents(mESC_merged, new.cluster.ids)
#Make a combined ident slot
mESC_merged@meta.data$combined.ident <- mESC_merged@active.ident

#Check UMAP again
p <- DimPlot(mESC_merged, reduction = "umap")
AugmentPlot(plot = p)
plot(p)

#Arrange cluster order and color clusters as desired
levels(x = mESC_merged) <- c("5","3","4","0","2")

p <- DimPlot(mESC_merged, reduction = "umap", cols = c("#5E4FA2","#3288BD","#2E8B57","#66C2A5","#ABDDA4"))
AugmentPlot(plot = p)
plot(p)

#UMAP split by medium condition
p <- DimPlot(mESC_merged, reduction = "umap", group.by = "merge.ident", cols = c("azure4","deepskyblue4","darkseagreen")) #merge.ident #fuse.ident?
AugmentPlot(plot = p)
plot(p)

#Re-calculate marker genes for clusters and for the medium conditions
mESC_merged_markers <- FindAllMarkers(mESC_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.35)
write.csv(mESC_merged_markers, file = "F:\\paper_2\\data_scRNAseq\\0h_fastMNN\\allmarkers_onlypos_03res_01merged_025_035.csv")

mESC_merged <- SetIdent(mESC_merged, value = mESC_merged@meta.data$merge.ident)
levels(x = mESC_merged) <- c("ESL0h","PD0h","a2i0h")
mESC_merged$merge.ident <- Idents(mESC_merged)

mESC_merged_markers <- FindAllMarkers(mESC_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
write.csv(mESC_merged_markers, file = "F:\\paper_2\\data_scRNAseq\\0h_fastMNN\\mergeidentmarkers_onlypos_03res_01merged_01_05.csv")

#Gene expression dotplot generation, depending on the active Ident this will be done either between the clusters or the medium conditions
cd_genes2 <- c("Zfp42","Dppa5a","Klf4","Tdh","Trh","Sox2","Nanog","Pou5f1","Utf1","Pim2","Dnmt3b","Pou3f1","Lefty1","Cdh1","Cdh2","Tagln","Krt8","Krt18","Mylpf","Mymx","Nodal","T","Fgf5","Fgf8","Wnt3","Eomes","Sox1","Sox3","Sox11","Calcoco2","Rpl39l","Ckb","Gpi1","Pgk1","Aldoa","mt-Nd3","Pdia3","Hspa5")
result = rev(cd_genes2)
dplot <- DotPlot(mESC_merged, features = result, cols = "Blues", scale.by = "size", scale = TRUE) + coord_flip()
dplot + scale_color_distiller(direction = 0)

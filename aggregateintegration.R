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


#Load Pre-CHIR99 datasets, set idents for later identification and correct data slot
ESL_48h_filtered <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\ESL_48h_filtered_allcc.rds")
a2i_72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_72h_filtered_nocc_morereads.rds")
PD_48h_filtered <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\PD_48h_forglobalint.rds")
PD_48h_filtered3 <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\PD_48h_forglobalint3.rds")

ESL_48h_filtered$fuse.ident <- "ESL_48h_1"
ESL_48h_filtered$merge.ident <- "ESL_48h"

a2i_72h_filtered$fuse.ident <- "a2i_72h_1"
a2i_72h_filtered$merge.ident <- "a2i_72h"

PD_48h_filtered$fuse.ident <- "PD_48h_1"
PD_48h_filtered$merge.ident <- "PD_48h"
PD_48h_filtered3$fuse.ident <- "PD_48h_2"
PD_48h_filtered3$merge.ident <- "PD_48h"

DefaultAssay(PD_48h_filtered) <- "RNA"
DefaultAssay(PD_48h_filtered3) <- "RNA"
DefaultAssay(ESL_48h_filtered) <- "RNA"
DefaultAssay(a2i_72h_filtered) <- "RNA" 


#Load Post-CHIR99 datasets, set idents for later identification and correct data slot
PD_72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_72h_filtered2022_nocc_morereads.rds")
a2i_96h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_96h_filtered_nocc_morereads.rds")
ESL_72h_filtered_new <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\ESL_72h_filtered_new_nocc.rds")
ESL_72h_filtered_old <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\ESL_72h_filtered_old_allcc.rds")

DefaultAssay(PD_72h_filtered) <- "RNA"
DefaultAssay(a2i_96h_filtered) <- "RNA"
DefaultAssay(ESL_72h_filtered_new) <- "RNA"
DefaultAssay(ESL_72h_filtered_old) <- "RNA" 

PD_72h_filtered$fuse.ident <- "PD_72h_1"
PD_72h_filtered$merge.ident <- "PD_72h"

a2i_96h_filtered$fuse.ident <- "a2i_96h_1"
a2i_96h_filtered$merge.ident <- "a2i_96h"

ESL_72h_filtered_new$fuse.ident <- "ESL_72h_1"
ESL_72h_filtered_new$merge.ident <- "ESL_72h"
ESL_72h_filtered_old$fuse.ident <- "ESL_72h_2"
ESL_72h_filtered_old$merge.ident <- "ESL_72h"


#Load Max. Elongation datasets
PD_120h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_120h_filtered_nocc_morereads.rds")
ESL_120h_filtered <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\ESL_120h_filtered_allcc.rds")
a2i_144h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_144h_filtered_nocc_morereads.rds")

PD_120h_filtered$fuse.ident <- "PD_120h_1"
ESL_120h_filtered$fuse.ident <- "ESL_120h_1"
a2i_144h_filtered$fuse.ident <- "a2i_144h_1"

#Add more idents
PD_120h_filtered$time.ident <- "elongated"
ESL_120h_filtered$time.ident <- "elongated"
a2i_144h_filtered$time.ident <- "elongated"

PD_72h_filtered$time.ident <- "postCHIR"
a2i_96h_filtered$time.ident <- "postCHIR"
ESL_72h_filtered_new$time.ident <- "postCHIR"
ESL_72h_filtered_old$time.ident <- "postCHIR"

ESL_48h_filtered$time.ident <- "preCHIR"
a2i_72h_filtered$time.ident <- "preCHIR"
PD_48h_filtered$time.ident <- "preCHIR"
PD_48h_filtered3$time.ident <- "preCHIR"

PD_120h_filtered$cond.ident <- "PD"
ESL_120h_filtered$cond.ident <- "SL"
a2i_144h_filtered$cond.ident <- "a2i"

PD_72h_filtered$cond.ident <- "PD"
a2i_96h_filtered$cond.ident <- "a2i"
ESL_72h_filtered_new$cond.ident <- "SL"
ESL_72h_filtered_old$cond.ident <- "SL"

ESL_48h_filtered$cond.ident <- "SL"
a2i_72h_filtered$cond.ident <- "a2i"
PD_48h_filtered$cond.ident <- "PD"
PD_48h_filtered3$cond.ident <- "PD"

#Merge all datasets to be integrated into a single object
all_objects <- merge(ESL_48h_filtered, y = c(PD_48h_filtered,PD_48h_filtered3,a2i_72h_filtered,ESL_72h_filtered_new,ESL_72h_filtered_old,PD_72h_filtered,a2i_96h_filtered,ESL_120h_filtered,PD_120h_filtered,a2i_144h_filtered), add.cell.ids = c('ESL48h','PD48h','PD48h2','a2i72h','ESL72h','ESL72h2','PD72h','a2i96h','ESL120h','PD120h','a2i144h'), project = "allmerged")

#Remove the "scale.data" slot to be sure, normalize data and find variable features collectively
all_objects@assays$RNA@scale.data <- as.matrix(0)
all_objects <- NormalizeData(all_objects)
all_objects <- FindVariableFeatures(all_objects, selection.method = "vst", nfeatures = 3000)

#Find features for integration and remove cell cycle genes according to Seurat lists (mouse orthologues)
features <- SelectIntegrationFeatures(object.list = SplitObject(all_objects, split.by = "orig.ident"), nfeatures = 3000)
features_filtered1 <- setdiff(features, g2m_genes)
features_filtered2 <- setdiff(features_filtered1, s_genes)

#Perform fastMNN integration
all_objects <- RunFastMNN(object.list = SplitObject(all_objects, split.by = "orig.ident"), features = features_filtered2)

#UMAP generation and clustering
all_objects <- RunUMAP(all_objects, reduction = "mnn", dims = 1:50)
all_objects <- FindNeighbors(all_objects, reduction = "mnn", dims = 1:50)
all_objects <- FindClusters(all_objects, resolution = 0.7)

DefaultAssay(all_objects) <- "RNA"

#Set color palettes
colors <- colorRampPalette(brewer.pal(11, "Spectral"))(16) 

blues <- brewer.pal(3, "Blues")
blues = rev(blues) 

#Change order of clusters
levels(x = all_objects) <- c("4","3","2","9","15","0","6","12","10","13","7","8","11","14","1","5")

#make UMAP with colors as in publication
p <- DimPlot(all_objects, reduction = "umap", cols = colors)
AugmentPlot(plot = p)
plot(p)
#UMAP color-coded according to timepoint
p <- DimPlot(all_objects, reduction = "umap", group.by = "time.ident", cols = blues) 
AugmentPlot(plot = p)
plot(p)
#UMAP color-coded according to source mESC medium condition
p <- DimPlot(all_objects, reduction = "umap", group.by = "cond.ident", cols = c("darkseagreen","deepskyblue4","azure4")) #merge.ident #fuse.ident?
AugmentPlot(plot = p)
plot(p)

#Example for gene expression plot on UMAP space
p <- FeaturePlot(all_objects, features = c("Zfp42","Dppa5a","Klf2","Sox2"), order = TRUE, cols = c("grey92", "dodgerblue4"), pt.size = 0.5)
AugmentPlot(plot = p)
plot(p)

#MVG identification, save table as .csv
all_objects_markers <- FindAllMarkers(all_objects, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
write.csv(all_objects_markers, file = "F:\\paper_2\\data_scRNAseq\\1_fastMNNint-all\\allmarkers_onlypos_05_025_07res.csv")

#Dotplot with gene expression
cd_genes2 <- c("Zfp42","Dppa5a","Klf2","Sox2","Apoe","Fabp3","Cd63","Utf1","Pim2","Dnmt3b","Pou3f1","Otx2","Lefty1","Nodal","mt-Cytb","mt-Nd5","Epcam","Cdh1","Cdh2","T","Fgf5","Fgf8","Wnt3","Nkx1-2","Cdx2","Eomes","Mixl1","Gsc","Lhx1","Mesp1","Foxa2","Chrd","Cobl","Bicc1","Tbx6","Tcf15","Notch1","Dll1","Aldh1a2","Meox1","Foxc1","Prrx2","Osr1","Meis2","Hand1","Gata6","Lmo2","Kdr","Cdh5","Etv2","Sox17","Spink1","Amot","Cldn6","Ncam1","Sox1","Sox3","Sox11","Pax3","Crabp1","Tubb3")
result = rev(cd_genes2)
dplot <- DotPlot(all_objects, features = result, cols = "Blues", scale.by = "radius", scale = TRUE) + coord_flip()
dplot + scale_color_distiller(direction = 0)

#Dotplot for the "merge.ident"
all_objects <- SetIdent(all_objects, value = all_objects@meta.data$merge.ident)
levels(x = all_objects) <- c("ESL_48h","PD_48h","a2i_72h","ESL_72h","PD_72h","a2i_96h","ESL120h","PD120h","a2i_144h")
all_objects$merge.ident <- Idents(all_objects)

cd_genes2 <- c("Zfp42","Dppa5a","Klf2","Sox2","Apoe","Fabp3","Cd63","Utf1","Pim2","Dnmt3b","Pou3f1","Otx2","Lefty1","Nodal","mt-Cytb","mt-Nd5","Epcam","Cdh1","Cdh2","T","Fgf5","Fgf8","Wnt3","Nkx1-2","Cdx2","Eomes","Mixl1","Gsc","Lhx1","Mesp1","Foxa2","Chrd","Cobl","Bicc1","Tbx6","Tcf15","Notch1","Dll1","Aldh1a2","Meox1","Foxc1","Prrx2","Osr1","Meis2","Hand1","Gata6","Lmo2","Kdr","Cdh5","Etv2","Sox17","Spink1","Amot","Cldn6","Ncam1","Sox1","Sox3","Sox11","Pax3","Crabp1","Tubb3")
result = rev(cd_genes2)
dplot <- DotPlot(all_objects, features = result, cols = "Blues", scale.by = "size", scale = TRUE) + coord_flip()
dplot + scale_color_distiller(direction = 0)

#Identify MVGs for each source mESC pluripotent state at each developmental timepoint
preCHIR <- subset(x = all_objects, subset = time.ident == "preCHIR")
postCHIR <- subset(x = all_objects, subset = time.ident == "postCHIR")
elo <- subset(x = all_objects, subset = time.ident == "elongated")

postCHIR <- SetIdent(postCHIR, value = postCHIR@meta.data$merge.ident)
postCHIR$merge.ident <- Idents(postCHIR)

preCHIR <- SetIdent(preCHIR, value = preCHIR@meta.data$merge.ident)
preCHIR$merge.ident <- Idents(preCHIR)

elo <- SetIdent(elo, value = elo@meta.data$merge.ident)
elo$merge.ident <- Idents(elo)

DefaultAssay(preCHIR) <- "RNA"
preCHIR_markers <- FindAllMarkers(preCHIR, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.35)
write.csv(preCHIR_markers, file = "F:\\paper_2\\data_scRNAseq\\fastMNNint-all\\preCHIR_onlypos_01_035_mergeident.csv")

DefaultAssay(postCHIR) <- "RNA"
postCHIR_markers <- FindAllMarkers(postCHIR, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.35)
write.csv(postCHIR_markers, file = "F:\\paper_2\\data_scRNAseq\\fastMNNint-all\\postCHIR_onlypos_01_035_mergeident.csv")

DefaultAssay(elo) <- "RNA"
elo_markers <- FindAllMarkers(elo, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.35)
write.csv(elo_markers, file = "F:\\paper_2\\data_scRNAseq\\fastMNNint-all\\elo_onlypos_01_035_mergeident.csv")
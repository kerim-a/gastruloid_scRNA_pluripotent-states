library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)
require(scales)

#Load datasets and set assay
PD_72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_72h_filtered2022_nocc_morereads.rds")
a2i_96h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_96h_filtered_nocc_morereads.rds")
ESL_72h_filtered_new <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\ESL_72h_filtered_new_nocc.rds")
ESL_72h_filtered_old <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\ESL_72h_filtered_old_allcc.rds")

DefaultAssay(PD_72h_filtered) <- "RNA"
DefaultAssay(a2i_96h_filtered) <- "RNA"
DefaultAssay(ESL_72h_filtered_new) <- "RNA"
DefaultAssay(ESL_72h_filtered_old) <- "RNA" 

#Merge if several replicates are present for a given condition
ESL_72h_filtered <- merge(ESL_72h_filtered_new, y = ESL_72h_filtered_old, add.cell.ids = c("ESL_new", "ESL_old"), project = "ESL_72hmerged", merge.data = TRUE)
DefaultAssay(ESL_72h_filtered) <- "RNA" 

#Filtering cells from each germ layer
mesoderm <- subset(ESL_72h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(ESL_72h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(ESL_72h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(ESL_72h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(ESL_72h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(ESL_72h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(ESL_72h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

#Assign a germ layer ident to each cell in the dataset
cells.use3 <- WhichCells(object = mesoderm)
ESL_72h_filtered <- SetIdent(object = ESL_72h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
ESL_72h_filtered <- SetIdent(object = ESL_72h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
ESL_72h_filtered <- SetIdent(object = ESL_72h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
ESL_72h_filtered <- SetIdent(object = ESL_72h_filtered, cells = cells.use4, value = 'pluripotent')

head(x = Idents(object = ESL_72h_filtered))
Idents(object = ESL_72h_filtered)

#Put this into metadata
ESL_72h_filtered$cool.ident <- Idents(ESL_72h_filtered)

#Check the ratios
length(WhichCells(object = ESL_72h_filtered, ident = 'mesoderm'))
2731/4131 #66.11
length(WhichCells(object = ESL_72h_filtered, ident = 'ectoderm'))
94/4131 #2.28
length(WhichCells(object = ESL_72h_filtered, ident = 'endoderm'))
330/4131 #7.99
length(WhichCells(object = ESL_72h_filtered, ident = 'pluripotent'))
972/4131 #23.53

ggplot(ESL_72h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


## Repeat this for the other datasets
# base+1i

mesoderm <- subset(PD_72h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(PD_72h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(PD_72h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(PD_72h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(PD_72h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(PD_72h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(PD_72h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

cells.use3 <- WhichCells(object = mesoderm)
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = cells.use4, value = 'pluripotent')

head(x = Idents(object = PD_72h_filtered))
Idents(object = PD_72h_filtered)

#Check the remaining, unassigned cells manually (scale data before if required)
PD_72h_filtered_markers <- FindAllMarkers(PD_72h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- PD_72h_filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(PD_72h_filtered, features = top10$gene) + NoLegend() + scale_fill_viridis()

rest_pluri <- WhichCells(object = PD_72h_filtered, ident = "0")
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = rest_pluri, value = 'pluripotent')
rest_pluri2 <- WhichCells(object = PD_72h_filtered, ident = "1")
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = rest_pluri2, value = 'pluripotent')
rest_pluri2_2 <- WhichCells(object = PD_72h_filtered, ident = "2")
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = rest_pluri2_2, value = 'pluripotent')
rest_pluri3 <- WhichCells(object = PD_72h_filtered, ident = "3")
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = rest_pluri3, value = 'pluripotent')
rest_pluri4 <- WhichCells(object = PD_72h_filtered, ident = "5")
PD_72h_filtered <- SetIdent(object = PD_72h_filtered, cells = rest_pluri4, value = 'pluripotent')

PD_72h_filtered$cool.ident <- Idents(PD_72h_filtered)

length(WhichCells(object = PD_72h_filtered, ident = 'mesoderm'))
2258/4279 # 52.77
length(WhichCells(object = PD_72h_filtered, ident = 'ectoderm'))
460/4279 # 10.75
length(WhichCells(object = PD_72h_filtered, ident = 'endoderm'))
2/4279 #0.05
length(WhichCells(object = PD_72h_filtered, ident = 'pluripotent'))
1559/4279 #36.43

ggplot(PD_72h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


##base+2i
mesoderm <- subset(a2i_96h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(a2i_96h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(a2i_96h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(a2i_96h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(a2i_96h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(a2i_96h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(a2i_96h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

cells.use3 <- WhichCells(object = mesoderm)
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = cells.use4, value = 'pluripotent')

a2i_96h_filtered_markers <- FindAllMarkers(a2i_96h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- a2i_96h_filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(a2i_96h_filtered, features = top10$gene) + NoLegend() + scale_fill_viridis()

rest_pluri <- WhichCells(object = a2i_96h_filtered, ident = "0")
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = rest_pluri, value = 'pluripotent')
rest_pluri2 <- WhichCells(object = a2i_96h_filtered, ident = "1")
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = rest_pluri2, value = 'pluripotent')
rest_pluri2_2 <- WhichCells(object = a2i_96h_filtered, ident = "2")
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = rest_pluri2_2, value = 'pluripotent')
rest_pluri3 <- WhichCells(object = a2i_96h_filtered, ident = "4")
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = rest_pluri3, value = 'pluripotent')
rest_pluri4 <- WhichCells(object = a2i_96h_filtered, ident = "6")
a2i_96h_filtered <- SetIdent(object = a2i_96h_filtered, cells = rest_pluri4, value = 'pluripotent')

head(x = Idents(object = a2i_96h_filtered))
Idents(object = a2i_96h_filtered)

a2i_96h_filtered$cool.ident <- Idents(a2i_96h_filtered)

length(WhichCells(object = a2i_96h_filtered, ident = 'mesoderm'))
1947/4363 # 44.63 (62)
length(WhichCells(object = a2i_96h_filtered, ident = 'ectoderm'))
584/4363 # 13.39
length(WhichCells(object = a2i_96h_filtered, ident = 'endoderm'))
11/4363 #0.25
length(WhichCells(object = a2i_96h_filtered, ident = 'pluripotent'))
1821/4363 #41.74

ggplot(a2i_96h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


#Plotting results from a .csv via ggplot
getwd()
setwd('C:\\Users\\anlas\\Documents\\Desktop_items')
getwd()

germlayer <- read.csv('R_germ_layer_ratios_postpol.csv')

#Sort order of x axis
germlayer$condition <- factor(germlayer$condition,levels = c("SL", "1i", "a2i"))

gg2 <- ggplot(germlayer, aes(fill=germ_layer, y=value, x=condition)) + 
  geom_bar(position="stack", stat="identity")

gg2 + theme_bw() + scale_fill_viridis(discrete = TRUE, option = "D")
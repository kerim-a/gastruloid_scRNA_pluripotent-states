library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)
require(scales)

#Load datasets
ESL_48h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\ESL_48h_filtered_allcc.rds")
ESL_48h_filtered <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\ESL_48h_filtered_allcc.rds")
a2i_72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_72h_filtered_nocc_morereads.rds")
PD_48h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_48h_forglobalint.rds")
PD_48h_filtered3 <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_48h_forglobalint3.rds")

DefaultAssay(a2i_72h_filtered) <- "RNA"
DefaultAssay(ESL_48h_filtered) <- "RNA"
DefaultAssay(PD_48h_filtered) <- "RNA"
DefaultAssay(PD_48h_filtered3) <- "RNA" 

#Merge into 1 for datasets with mutliple replicates 
PD_48h_filtered_merged <- merge(PD_48h_filtered, y = PD_48h_filtered3, add.cell.ids = c("PD_new", "PD_old"), project = "PD_48hmerged", merge.data = TRUE)
DefaultAssay(PD_48h_filtered_merged) <- "RNA" 

#Apply gene expression filters
mesoderm <- subset(ESL_48h_filtered, subset = Aldh1a2 > 0.2 | Meox1 > 0.2 | Tbx6 > 0.2 | Hand1 > 0.2 | Hand2 > 0.2 | Kdr > 0.2 | Cdx2 > 0.2 | Eomes > 0.2 | Meis1 > 0.2 | Lhx1 > 0.2 | T > 0.2 | Bicc1 > 0.2)

endoderm <- subset(ESL_48h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(ESL_48h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(ESL_48h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(ESL_48h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(ESL_48h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(ESL_48h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

#Designate an ident (i.e. germ layer) for each cell in the object. The order is important
cells.use3 <- WhichCells(object = mesoderm)
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = cells.use4, value = 'pluripotent')

#Find MVGs and check a heatmap to identify thus far unassigned cells
#Scale data if necessary (for heatmap generation)
all.genes <- rownames(ESL_48h_filtered_markers)
ESL_48h_filtered_markers <- ScaleData(ESL_48h_filtered_markers, features = all.genes)

ESL_48h_filtered_markers <- FindAllMarkers(ESL_48h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- ESL_48h_filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(ESL_48h_filtered, features = top10$gene) + NoLegend() + scale_fill_viridis()
write.csv(ESL_48h_filtered_markers, file = "C:\\Users\\anlas\\Desktop\\temp\\ESL_48h.csv")

#Following visual inspection, assign remaining cells to a germ layer (if applicable)
rest_pluri <- WhichCells(object = ESL_48h_filtered, ident = "0")
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = rest_pluri, value = 'pluripotent')
rest_pluri2 <- WhichCells(object = ESL_48h_filtered, ident = "1")
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = rest_pluri2, value = 'pluripotent')
rest_pluri2_2 <- WhichCells(object = ESL_48h_filtered, ident = "3")
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = rest_pluri2_2, value = 'pluripotent')
rest_pluri3 <- WhichCells(object = ESL_48h_filtered, ident = "4")
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = rest_pluri3, value = 'pluripotent')
rest_ecto <- WhichCells(object = ESL_48h_filtered, ident = "6")
ESL_48h_filtered <- SetIdent(object = ESL_48h_filtered, cells = rest_ecto, value = 'ectoderm')

#Specifying a designated metadata slot can help with plotting later
ESL_48h_filtered$cool.ident <- Idents(ESL_48h_filtered)

#This is to check the percentages for each germ layer
length(WhichCells(object = ESL_48h_filtered, ident = 'mesoderm'))
1516/6220 #24.37, i.e. "mesoderm" cells / total cells in object
length(WhichCells(object = ESL_48h_filtered, ident = 'ectoderm'))
283/6220 #4.55
length(WhichCells(object = ESL_48h_filtered, ident = 'endoderm'))
6/6220 #0.1
length(WhichCells(object = ESL_48h_filtered, ident = 'pluripotent'))
4413/6220 #70.95
#Optionally plot/check via ggplot
ggplot(ESL_48h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


### Repeat this for all datasets of interest!

mesoderm <- subset(PD_48h_filtered_merged, subset = Aldh1a2 > 0.2 | Meox1 > 0.2 | Tbx6 > 0.2 | Hand1 > 0.2 | Hand2 > 0.2 | Kdr > 0.2 | Cdx2 > 0.2 | Eomes > 0.2 | Meis1 > 0.2 | Lhx1 > 0.2 | T > 0.2 | Bicc1 > 0.2)

endoderm <- subset(PD_48h_filtered_merged, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(PD_48h_filtered_merged, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(PD_48h_filtered_merged, subset = Sox1 > 0.2)
ectoderm2 <- subset(PD_48h_filtered_merged, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(PD_48h_filtered_merged, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(PD_48h_filtered_merged, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

cells.use3 <- WhichCells(object = mesoderm)
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = c(cells.use,cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = cells.use4, value = 'pluripotent')

head(x = Idents(object = PD_48h_filtered_merged))
Idents(object = PD_48h_filtered_merged)

all.genes <- rownames(PD_48h_filtered_merged)
PD_48h_filtered_merged <- ScaleData(PD_48h_filtered_merged, features = all.genes)

PD_48h_filtered_merged_markers <- FindAllMarkers(PD_48h_filtered_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- PD_48h_filtered_merged_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(PD_48h_filtered_merged, features = top10$gene) + NoLegend() + scale_fill_viridis()

rest_pluri <- WhichCells(object = PD_48h_filtered_merged, ident = "PD_48h")
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = rest_pluri, value = 'pluripotent')
rest_pluri2 <- WhichCells(object = PD_48h_filtered_merged, ident = "1")
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = rest_pluri2, value = 'pluripotent')
rest_pluri2_2 <- WhichCells(object = PD_48h_filtered_merged, ident = "2")
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = rest_pluri2_2, value = 'pluripotent')
rest_pluri3 <- WhichCells(object = PD_48h_filtered_merged, ident = "3")
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = rest_pluri3, value = 'pluripotent')
rest_pluri4 <- WhichCells(object = PD_48h_filtered_merged, ident = "5")
PD_48h_filtered_merged <- SetIdent(object = PD_48h_filtered_merged, cells = rest_pluri4, value = 'pluripotent')

PD_48h_filtered_merged$cool.ident <- Idents(PD_48h_filtered_merged)

length(WhichCells(object = PD_48h_filtered_merged, ident = 'mesoderm'))
37/3892 # 1.0
length(WhichCells(object = PD_48h_filtered_merged, ident = 'ectoderm'))
45/3892 # 1.16
length(WhichCells(object = PD_48h_filtered_merged, ident = 'endoderm'))
0/3892 #0
length(WhichCells(object = PD_48h_filtered_merged, ident = 'pluripotent'))
3787/3892 #97.3

ggplot(PD_48h_filtered_merged@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


##base+2i

mesoderm <- subset(a2i_72h_filtered, subset = Aldh1a2 > 0.2 | Meox1 > 0.2 | Tbx6 > 0.2 | Hand1 > 0.2 | Hand2 > 0.2 | Cdx2 > 0.2 | Eomes > 0.2 | Meis1 > 0.2 | Lhx1 > 0.2 | T > 0.2 | Bicc1 > 0.2)

endoderm <- subset(a2i_72h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05)
endoderm2 <- subset(a2i_72h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(a2i_72h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(a2i_72h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2) 
ectoderm3 <- subset(a2i_72h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(a2i_72h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

cells.use3 <- WhichCells(object = mesoderm)
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = cells.use4, value = 'pluripotent')

all.genes <- rownames(a2i_72h_filtered)
a2i_72h_filtered <- ScaleData(a2i_72h_filtered, features = all.genes)

a2i_72h_filtered_markers <- FindAllMarkers(a2i_72h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- a2i_72h_filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(a2i_72h_filtered, features = top10$gene) + NoLegend() + scale_fill_viridis()

rest_pluri <- WhichCells(object = a2i_72h_filtered, ident = "0")
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = rest_pluri, value = 'pluripotent')
rest_pluri2 <- WhichCells(object = a2i_72h_filtered, ident = "1")
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = rest_pluri2, value = 'pluripotent')
rest_pluri2_2 <- WhichCells(object = a2i_72h_filtered, ident = "2")
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = rest_pluri2_2, value = 'pluripotent')
rest_pluri3 <- WhichCells(object = a2i_72h_filtered, ident = "4")
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = rest_pluri3, value = 'pluripotent')
rest_pluri4 <- WhichCells(object = a2i_72h_filtered, ident = "6")
a2i_72h_filtered <- SetIdent(object = a2i_72h_filtered, cells = rest_pluri4, value = 'pluripotent')


head(x = Idents(object = a2i_72h_filtered))
Idents(object = a2i_72h_filtered)

a2i_72h_filtered$cool.ident <- Idents(a2i_72h_filtered)

length(WhichCells(object = a2i_72h_filtered, ident = 'mesoderm'))
23/3381 # 0.68
length(WhichCells(object = a2i_72h_filtered, ident = 'ectoderm'))
2/3381 # 0.06
length(WhichCells(object = a2i_72h_filtered, ident = 'endoderm'))
0/3381 #0
length(WhichCells(object = a2i_72h_filtered, ident = 'pluripotent'))
3356/3381 #99.26

ggplot(a2i_72h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


#For plotting results in .csv format
getwd()
setwd('C:\\Users\\anlas\\Documents\\Desktop_items')
getwd()

germlayer <- read.csv('R_germ_layer_ratios_prepol.csv')

#Sort order of x axis
germlayer$condition <- factor(germlayer$condition,levels = c("SL", "1i", "a2i"))

gg3 <- ggplot(germlayer, aes(fill=germ_layer, y=value, x=condition)) + 
  geom_bar(position="stack", stat="identity")

gg3 + theme_bw() + scale_fill_viridis(discrete = TRUE, option = "D")
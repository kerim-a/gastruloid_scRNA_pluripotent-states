library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)
require(scales)

#Load datasets
PD_120h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\PD_120h_filtered_nocc_morereads.rds")
ESL_120h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\ESL_120h_filtered_allcc.rds")
a2i_144h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\a2i_144h_filtered_nocc_morereads.rds")

#Set correct assay slot
DefaultAssay(PD_120h_filtered) <- "RNA"
DefaultAssay(ESL_120h_filtered) <- "RNA"
DefaultAssay(a2i_144h_filtered) <- "RNA"

#Perform filtering
mesoderm <- subset(ESL_120h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(ESL_120h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(ESL_120h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(ESL_120h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(ESL_120h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(ESL_120h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(ESL_120h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

#Assign a germ layer identity, order is important
cells.use3 <- WhichCells(object = mesoderm)
ESL_120h_filtered <- SetIdent(object = ESL_120h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)

ESL_120h_filtered <- SetIdent(object = ESL_120h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
ESL_120h_filtered <- SetIdent(object = ESL_120h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
ESL_120h_filtered <- SetIdent(object = ESL_120h_filtered, cells = cells.use4, value = 'pluripotent')

#Put this into a metadata slot
head(x = Idents(object = ESL_120h_filtered))
Idents(object = ESL_120h_filtered)
ESL_120h_filtered$cool.ident <- Idents(ESL_120h_filtered)

#Check results/ratios
length(WhichCells(object = ESL_120h_filtered, ident = 'mesoderm'))
4019/5991 #67.01
length(WhichCells(object = ESL_120h_filtered, ident = 'ectoderm'))
1332/5991 #22.23
length(WhichCells(object = ESL_120h_filtered, ident = 'endoderm'))
453/5991 #7.6
length(WhichCells(object = ESL_120h_filtered, ident = 'pluripotent'))
165/5991 #2.8

ggplot(ESL_120h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


###Repeat this for all other datasets

##base+1i
mesoderm <- subset(PD_120h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(PD_120h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) 
endoderm2 <- subset(PD_120h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) 

ectoderm <- subset(PD_120h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(PD_120h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(PD_120h_filtered, subset = Hes6 > 1 & Tubb3 > 1)

pluripotent <- subset(PD_120h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

cells.use3 <- WhichCells(object = mesoderm)
PD_120h_filtered <- SetIdent(object = PD_120h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)

PD_120h_filtered <- SetIdent(object = PD_120h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
PD_120h_filtered <- SetIdent(object = PD_120h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
PD_120h_filtered <- SetIdent(object = PD_120h_filtered, cells = cells.use4, value = 'pluripotent')

head(x = Idents(object = PD_120h_filtered))
Idents(object = PD_120h_filtered)
PD_120h_filtered$cool.ident <- Idents(PD_120h_filtered)

length(WhichCells(object = PD_120h_filtered, ident = 'mesoderm'))
1929/4819 # 40.02
length(WhichCells(object = PD_120h_filtered, ident = 'ectoderm'))
2426/4819 # 50.34
length(WhichCells(object = PD_120h_filtered, ident = 'endoderm'))
185/4819 #3.84
length(WhichCells(object = PD_120h_filtered, ident = 'pluripotent'))
267/4819 #5.54

ggplot(PD_120h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


#base+2i
mesoderm <- subset(a2i_144h_filtered, subset = Aldh1a2 > 0.1 | Meox1 > 0.1 | Tbx6 > 0.1 | Hand1 > 0.1 | Hand2 > 0.1 | Kdr > 0.1 | Cdx2 > 0.1 | Eomes > 0.1 | Meis1 > 0.1 | Lhx1 > 0.1 | T > 0.1 | Bicc1 > 0.1)

endoderm <- subset(a2i_144h_filtered, subset = Spink1 > 0.2 & Cldn6 > 0.05) # or use vtn instead of cldn6
endoderm2 <- subset(a2i_144h_filtered, subset = Sox17 > 0.1 & Cldn6 > 0.05) # or just sox17?
# endoderm_total <- merge(endoderm, y = endoderm2, add.cell.ids = c("endo1", "endo2"), project = "endo_merge", merge.data = TRUE)
# the merging doesnt work with ident setting for some reason

ectoderm <- subset(a2i_144h_filtered, subset = Sox1 > 0.2)
ectoderm2 <- subset(a2i_144h_filtered, subset = Epha5 > 0.2 & Ncam1 > 0.5 & Sox2 > 0.2)
ectoderm3 <- subset(a2i_144h_filtered, subset = Hes6 > 1 & Tubb3 > 1)
# ectoderm_total <- merge(ectoderm, y = ectoderm2, add.cell.ids = c("ecto1", "ecto2"), project = "ecto_merge", merge.data = TRUE)
# the merging doesnt work with ident setting for some reason

pluripotent <- subset(a2i_144h_filtered, subset = Zfp42 > 0.1 | Dppa5a > 0.5)

#order is important! First meso, then ecto, then endo, then pluri!!!
cells.use3 <- WhichCells(object = mesoderm)
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = cells.use3, value = 'mesoderm')

cells.use2 <- WhichCells(object = ectoderm)
cells.use_2 <- WhichCells(object = ectoderm2)
cells.use_e2 <- WhichCells(object = ectoderm3)
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = c(cells.use2, cells.use_2, cells.use_e2), value = 'ectoderm')

cells.use <- WhichCells(object = endoderm)
cells.use_1 <- WhichCells(object = endoderm2)
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = c(cells.use, cells.use_1), value = 'endoderm')

cells.use4 <- WhichCells(object = pluripotent)
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = cells.use4, value = 'pluripotent')

#Check remaining cells (perform ScaleData before if required)
a2i_144h_filtered_markers <- FindAllMarkers(a2i_144h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- a2i_144h_filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(a2i_144h_filtered, features = top10$gene) + NoLegend() + scale_fill_viridis()

#Assign remaining cells (if feasible)
rest_meso <- WhichCells(object = a2i_144h_filtered, ident = "2")
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = rest_meso, value = 'mesoderm')
rest_meso2 <- WhichCells(object = a2i_144h_filtered, ident = "1")
a2i_144h_filtered <- SetIdent(object = a2i_144h_filtered, cells = rest_meso2, value = 'mesoderm')

head(x = Idents(object = a2i_144h_filtered))
Idents(object = a2i_144h_filtered)
a2i_144h_filtered$cool.ident <- Idents(a2i_144h_filtered)

length(WhichCells(object = a2i_144h_filtered, ident = 'mesoderm'))
1629/5548 # 29.36
length(WhichCells(object = a2i_144h_filtered, ident = 'ectoderm'))
3277/5548 # 59.07
length(WhichCells(object = a2i_144h_filtered, ident = 'endoderm'))
284/5548 #5.12
length(WhichCells(object = a2i_144h_filtered, ident = 'pluripotent'))
307/5548 #5.53

ggplot(a2i_144h_filtered@meta.data, aes(x=cool.ident, fill=merge.ident)) + geom_bar(aes(y = (..count..)/sum(..count..)))


#Plot results from a .csv in ggplot
getwd()
setwd('C:\\Users\\anlas\\Documents\\Desktop_items')
getwd()

germlayer <- read.csv('R_germ_layer_ratios.csv')

#Set order of x-axis
germlayer$condition <- factor(germlayer$condition,levels = c("SL", "1i", "a2i"))

gg1 <- ggplot(germlayer, aes(fill=germ_layer, y=value, x=condition)) + 
  geom_bar(position="stack", stat="identity")

gg1 + theme_bw() + scale_fill_viridis(discrete = TRUE, option = "D")
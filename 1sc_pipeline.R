############
library(devtools)
library(BiocManager)
library(tidyverse)
library(Seurat)
library(cowplot)
library(ggplot2)
library(SingleR)
library(monocle)
library(DESeq2)
library(clusterProfiler)
library(DOSE)
library(pheatmap)
library(org.Hs.eg.db)
library(MAST) 
library(reshape2)
library(SeuratObject)
library(clustree)
library(harmony)
library(celldex)
library(DoubletFinder) 
library(EnhancedVolcano)
library(monocle3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DO.db)
library(R.utils)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(Scillus)
library(infercnv)
library(AnnoProbe)
library(ggthemes)
library(scRNAtoolVis)
library(AUCell)
library(msigdbr)
library(ggrepel)
library(CellChat)
library(ggforce)








###### read data ######
setwd("~/")

folders = list.files('./',pattern='[123456789]$') 
folders

scList = lapply(folders,function(folder){
  CreateSeuratObject(counts = Read10X(folder),
                     project = folder,
                     min.cells = 5, min.features = 300)
})

# remove double cells
for (i in 1:length(scList)) 
{
  # preprocessing
  data = NormalizeData(scList[[i]])
  data = ScaleData(data, verbose = FALSE)
  data = FindVariableFeatures(data, verbose = FALSE)
  data = RunPCA(data, npcs = 30, verbose = FALSE)
  data = RunUMAP(data, reduction = "pca", dims = 1:30)
  # Doublet
  sweep.res.list = paramSweep(data, PCs = 1:30)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn = find.pK(sweep.stats)
  pk_bcmvn = bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop = modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi = round(0.075*nrow(data@meta.data)) 
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
  data = doubletFinder(data, PCs = 1:30, pN = 0.25, pK = pk_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  # save result
  scList[[i]]$doubFind_res = data@meta.data %>% dplyr::select(contains('DF.classifications'))
  scList[[i]]$doubFind_score = data@meta.data %>% dplyr::select(contains('pANN'))
  
  print(i)
  
}


BD_seu = merge(scList[[1]],
                      y = c(scList[[2]],scList[[3]],scList[[4]],scList[[5]],scList[[6]],
                            scList[[7]],scList[[8]],scList[[9]],scList[[10]],scList[[11]],
                            scList[[12]],scList[[13]]),
                      project = "BD")   
BD_seu
table(BD_seu$orig.ident)

## remove BD8
BD_seu = subset(BD_seu, orig.ident %in% c("BD1","BD2","BD3","BD4","BD5","BD6","BD7","BD9","HC1","HC2","HC3","HC4"))
rm(scList)
table(BD_seu$orig.ident)
BD_seu = subset(BD_seu,doubFind_res == "Singlet")




###### add group information ######
metatable = read.csv(file = "~/group_change.csv")
metadata = FetchData(BD_seu, "orig.ident")
metadata = left_join(x = metadata, y = metatable, by = "orig.ident")
rownames(metadata) = rownames(BD_seu@meta.data)
BD_seu = AddMetaData(BD_seu, metadata = metadata)




###### pipeline ######
BD_seu[["percent.rb"]] = PercentageFeatureSet(BD_seu, pattern = "^RP[SL]")
BD_seu[["percent.hb"]] = PercentageFeatureSet(BD_seu, pattern = "^HB[^(P)]")
BD_seu[["percent.mt"]] = PercentageFeatureSet(BD_seu, pattern = "^MT-")
VlnPlot(BD_seu, features = c("nFeature_RNA"),
        ncol = 1,
        group.by = "orig.ident",
        pt.size = 0,raster=FALSE)
VlnPlot(BD_seu, features = c("nCount_RNA"),
        ncol = 1,
        group.by = "orig.ident",
        pt.size = 0,raster=FALSE)
VlnPlot(BD_seu, features = c("percent.mt"),
        ncol = 1,
        group.by = "orig.ident",
        pt.size = 0,raster=FALSE)
VlnPlot(BD_seu, features = c("percent.rb"),
        ncol = 1,
        group.by = "orig.ident",
        pt.size = 0,raster=FALSE)
VlnPlot(BD_seu, features = c("percent.hb"),
        ncol = 1,
        group.by = "orig.ident",
        pt.size = 0,raster=FALSE)
VlnPlot(BD_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"),
        ncol = 3,
        group.by = "orig.ident",
        pt.size = 0)


###
BD_seu = subset(BD_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5 & percent.hb < 1)
print(BD_seu)

VlnPlot(BD_seu, features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2,
        group.by = "orig.ident",
        pt.size = 0,
        raster=FALSE)
VlnPlot(BD_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"),
        ncol = 5,
        group.by = "orig.ident",
        pt.size = 0)
FeatureScatter(BD_seu,feature1 = "nCount_RNA",feature2 = "nFeature_RNA" )





################
##
BD_seu = NormalizeData(BD_seu)
BD_seu = FindVariableFeatures(BD_seu)
DefaultAssay(BD_seu)

## 
BD_seu = ScaleData(BD_seu, features = rownames(BD_seu))
BD_seu = RunPCA(BD_seu, features = VariableFeatures(object = BD_seu),reduction.name = "pca")
ElbowPlot(BD_seu)
DimPlot(BD_seu, reduction = "pca")

##
BD_seu = RunHarmony(BD_seu,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(BD_seu, reduction = "harmony") 
BD_seu = RunUMAP(BD_seu, reduction = "harmony", dims = 1:30, reduction.name = "umap")
DimPlot(BD_seu,reduction = "umap",group.by = "orig.ident", label = T, repel = T)
DimPlot(BD_seu,split.by = "orig.ident",ncol = 4, repel = T, raster = F)
DimPlot(BD_seu,split.by = "group",ncol = 2, repel = T, raster = F)

##
BD_seu = FindNeighbors(BD_seu, reduction = "harmony", dims = 1:30)
BD_seu = FindClusters(BD_seu, resolution = 0.6)  
# BD_seu = FindClusters(BD_seu, resolution = seq(0.1,1.0,0.1))
# clustree(BD_seu)

##
BD_seu = RunUMAP(BD_seu, reduction = "harmony", dims = 1:30)
DimPlot(BD_seu,reduction = "umap",group.by = "seurat_clusters",label = T)

##
DimPlot(BD_seu,split.by = "orig.ident",ncol = 4, repel = T, raster = F)
DimPlot(BD_seu,split.by = "group",ncol = 2, repel = T, raster = F,label = T)
#
save(BD_seu, file = "~/clustering.rda")
#



































 

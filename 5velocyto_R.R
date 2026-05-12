library(Seurat)
library(Matrix)


base::load("~/annotation.rda")

seurat_object = BD_seu
rm(BD_seu,BD_seu.markers)

## (0) remove the LDN
table(seurat_object@meta.data$celltype)
seurat_object = seurat_object[,Idents(seurat_object) %in% c("Neu_01_CXCL8","Neu_02_TANK","Neu_03_S100A12","Neu_04_IFIT1","Neu_05_STAT1","Neu_06_MMP9")]
table(seurat_object@meta.data$celltype)
set.seed(123456)
seurat_object = subset(seurat_object, downsample = 5000)
table(Idents(seurat_object))
table(seurat_object@meta.data$celltype)


## (1) the RNA count matrix
counts = GetAssayData(seurat_object, slot = "counts", assay = "RNA")  
writeMM(counts, file = "~/counts.mtx")   # writeMM: read and write external matrix formats

# 
write.csv(rownames(counts), file = "~/genes.csv", row.names = FALSE)

# 
write.csv(colnames(counts), file = "~/barcodes.csv", row.names = FALSE)


## (2) metadata of cells
metadata = seurat_object@meta.data
write.csv(metadata, file = "~/metadata.csv", row.names = TRUE)


## (3) the reduction coordinates
umap_coords = Embeddings(seurat_object, reduction = "umap")
pca_coords = Embeddings(seurat_object, reduction = "pca")
#tsne_coords = Embeddings(seurat_object, reduction = "tsne")
harmony_coords = Embeddings(seurat_object, reduction = "harmony")

# 
write.csv(umap_coords, file = "~/umap_coords.csv", row.names = TRUE)
write.csv(pca_coords, file = "~/pca_coords.csv", row.names = TRUE)
#write.csv(tsne_coords, "tsne_coords.csv", row.names = TRUE)
write.csv(harmony_coords, file = "~/harmony_coords.csv", row.names = TRUE)














########################## BD vs HC ####################################
base::load("~/annotation.rda")
seurat_object = BD_seu
rm(BD_seu,BD_seu.markers)
## (0)
Idents(seurat_object) = "celltype"
table(seurat_object@meta.data$celltype)
table(seurat_object@meta.data$group)
seurat_object = seurat_object[,Idents(seurat_object) %in% c("Neu_01_CXCL8","Neu_02_TANK","Neu_03_S100A12","Neu_04_IFIT1","Neu_05_STAT1","Neu_06_MMP9")]

# BD or HC
Idents(seurat_object) = "group"
# BD
seurat_object_bd = seurat_object[,Idents(seurat_object) %in% c("BD")]
table(seurat_object_bd@meta.data$celltype)
table(seurat_object_bd@meta.data$group)
# HC
seurat_object_hc = seurat_object[,Idents(seurat_object) %in% c("HC")]
table(seurat_object_hc@meta.data$celltype)
table(seurat_object_hc@meta.data$group)
# 
Idents(seurat_object_bd) = "celltype"
Idents(seurat_object_hc) = "celltype"
set.seed(123456)
seurat_object_bd = subset(seurat_object_bd, downsample = 2000)
seurat_object_hc = subset(seurat_object_hc, downsample = 2000)














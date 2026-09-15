library(Seurat)
library(Matrix)


base::load("~/annotation.rda")
seurat_object = BD_seu
rm(BD_seu,BD_seu.markers)

## remove the LDN
table(seurat_object@meta.data$celltype)
seurat_object = seurat_object[,Idents(seurat_object) %in% c("Neu_01_CXCL8","Neu_02_TANK","Neu_03_S100A12","Neu_04_IFIT1","Neu_05_STAT1","Neu_06_MMP9")]
table(seurat_object@meta.data$celltype)
# 
table(Idents(seurat_object))
Idents(seurat_object) = "group"

bd = seurat_object[,Idents(seurat_object) %in% c("BD")]
Idents(bd) = "celltype"
table(Idents(bd))

hc = seurat_object[,Idents(seurat_object) %in% c("HC")]
Idents(hc) = "celltype"
table(Idents(hc))
##
set.seed(123456)
bd = subset(bd, downsample = 4000)
table(Idents(bd))
table(bd@meta.data$celltype)

hc = subset(hc, downsample = 4000)
table(Idents(hc))
table(hc@meta.data$celltype)

## the RNA count matrix
counts1 = GetAssayData(bd, slot = "counts", assay = "RNA")  
counts2 = GetAssayData(hc, slot = "counts", assay = "RNA") 
writeMM(counts1, file = "~/counts_bd.mtx")   # writeMM: read and write external matrix formats
writeMM(counts2, file = "~/counts_hc.mtx")   # writeMM: read and write external matrix formats

write.csv(rownames(counts1), file = "~/genes_bd.csv", row.names = FALSE)
write.csv(rownames(counts2), file = "~/genes_hc.csv", row.names = FALSE)

write.csv(colnames(counts1), file = "~/barcodes_bd.csv", row.names = FALSE)
write.csv(colnames(counts2), file = "~/barcodes_hc.csv", row.names = FALSE)


## metadata of cells
metadata1 = bd@meta.data
metadata2 = hc@meta.data
write.csv(metadata1, file = "~/metadata_bd.csv", row.names = TRUE)
write.csv(metadata2, file = "~/metadata_hc.csv", row.names = TRUE)


## the reduction coordinates
umap_coords1 = Embeddings(bd, reduction = "umap")
pca_coords1 = Embeddings(bd, reduction = "pca")
#tsne_coords = Embeddings(seurat_object, reduction = "tsne")
harmony_coords1 = Embeddings(bd, reduction = "harmony")

umap_coords2 = Embeddings(hc, reduction = "umap")
pca_coords2 = Embeddings(hc, reduction = "pca")
#tsne_coords = Embeddings(seurat_object, reduction = "tsne")
harmony_coords2 = Embeddings(hc, reduction = "harmony")

# 
write.csv(umap_coords1, file = "~/umap_coords_bd.csv", row.names = TRUE)
write.csv(pca_coords1, file = "~/pca_coords_bd.csv", row.names = TRUE)
#write.csv(tsne_coords, "tsne_coords.csv", row.names = TRUE)
write.csv(harmony_coords1, file = "~/harmony_coords_bd.csv", row.names = TRUE)

write.csv(umap_coords2, file = "~/umap_coords_hc.csv", row.names = TRUE)
write.csv(pca_coords2, file = "~/pca_coords_hc.csv", row.names = TRUE)
#write.csv(tsne_coords, "tsne_coords.csv", row.names = TRUE)
write.csv(harmony_coords2, file = "~/harmony_coords_hc.csv", row.names = TRUE)






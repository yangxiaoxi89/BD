############ seurat  V5
# 

#################################### 
##### scmetabolism

### meta data
base::load("~/annotation.rda")
Idents(BD_seu) = "orig.ident"
BD_test = subset(x = BD_seu, downsample = 3000)
table(BD_test$orig.ident)
table(BD_test$celltype)
dim(BD_test)
Idents(BD_test) = "celltype"
rm(BD_seu.markers,BD_seu)
gc()



### score
DefaultAssay(BD_test) = "RNA"
countexp.seurat = sc.metabolism.Seurat(obj = BD_test,
                                       method = "AUCell",  
                                       imputation = FALSE,
                                       ncores = 2,
                                       metabolism.type = "KEGG")
# countexp.seurat1 = sc.metabolism.Seurat(obj = BD_test,
#                                        method = "AUCell",  
#                                        imputation = FALSE,
#                                        ncores = 2,
#                                        metabolism.type = "REACTOME")  



### 
# 
score = countexp.seurat@assays$METABOLISM$score

# 
BD_test@meta.data = cbind(BD_test@meta.data, t(score))
save(BD_test,countexp.seurat, file = "~/scMetabloism.rda")



### plot
DimPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Glycolysis / Gluconeogenesis",
  dimention.reduction.type = "umap",
  dimention.reduction.run = FALSE,
  size = 1)
DimPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Pentose phosphate pathway",
  dimention.reduction.type = "umap",
  dimention.reduction.run = FALSE,
  size = 1)
DimPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Oxidative phosphorylation",
  dimention.reduction.type = "umap",
  dimention.reduction.run = FALSE,
  size = 1)
DimPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Galactose metabolism",
  dimention.reduction.type = "umap",
  dimention.reduction.run = FALSE,
  size = 1)
DimPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Fructose and mannose metabolism",
  dimention.reduction.type = "umap",
  dimention.reduction.run = FALSE,
  size = 1)



### 
library(ggsci)

# "Glycolysis / Gluconeogenesis",  
# "Pentose phosphate pathway",   
# "Oxidative phosphorylation",  
# "Galactose metabolism",   
# "Fructose and mannose metabolism"   
# "Fatty acid elongation" 
# "Fatty acid biosynthesis"   
# "Fatty acid degradation"


BoxPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Fatty acid degradation",  # 
  phenotype = "celltype",
  ncol = 1) +
  scale_fill_nejm()



# heatmap
# 
df = BD_test@meta.data
# 
avg_df = aggregate(df[,14:ncol(df)],
                   list(df$celltype),
                   mean)
# 
avg_df <- avg_df %>% 
  dplyr::select(1:21) %>%     
  column_to_rownames("Group.1") 
avg_df[1:4,1:4]

avg_df = as.matrix(avg_df)

pheatmap(t(avg_df), 
         show_colnames = T,
         scale = 'row', 
         cluster_rows = T, cutree_rows = 3,
         border_color = "white",   
         color = colorRampPalette(c('#373F89','white',"#ed8b10"))(100),
         cluster_cols = T, cutree_cols = 3,
         treeheight_col = 10,   
         treeheight_row = 10,
         angle_col = "45")






















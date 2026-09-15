############ 
# 
install.packages("Seurat")
packageVersion("Seurat")
packageVersion("SeuratObject")
# 
remove.packages(c("Seurat","SeuratObject"))
# 
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'))
remotes::install_github("satijalab/seurat", ref = "v4.4.0")
install.packages("https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.4.tar.gz", repos = NULL, type = "source")

# 
.libPaths(c("~/SeuratV4", .libPaths()))
library(Seurat)
# 
packageVersion("Seurat")
packageVersion("SeuratObject")




##### 
gc()
Idents(BD_seu) = "orig.ident"
BD_test = subset(x = BD_seu, downsample = 3000)
table(BD_test$orig.ident)
table(BD_test$celltype)
dim(BD_test)
Idents(BD_test) = "celltype"
rm(BD_seu.markers,BD_seu)
gc()



### 
# 
DefaultAssay(BD_test) = "RNA"
countexp.seurat = sc.metabolism.Seurat(obj = BD_test,
                                       method = "AUCell",  
                                       imputation = FALSE,
                                       ncores = 2,
                                       metabolism.type = "KEGG")  

score = countexp.seurat@assays$METABOLISM$score
BD_test@meta.data = cbind(BD_test@meta.data, t(score))
# 
save(BD_test,countexp.seurat, file = "~/scMetabloism.rda")



### 
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
input.pathway = c(
  "Glycolysis / Gluconeogenesis",
  "Pentose phosphate pathway",
  "Oxidative phosphorylation",
  "Galactose metabolism",
  "Fructose and mannose metabolism"
)
# 
input.pathway = rownames(countexp.seurat@assays$METABOLISM$score)[1:20]
library(ggsci)

BoxPlot.metabolism(
  obj = countexp.seurat,
  pathway = "Fatty acid degradation",  
  phenotype = "celltype",
  ncol = 1) +
  scale_fill_nejm()


# 
df = BD_test@meta.data
avg_df = aggregate(df[,14:ncol(df)],
                   list(df$celltype),
                   mean)
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


################## FAO ##################
# 
"Fatty acid degradation" ("hsa00071")  
"Fatty acid metabolism" ("hsa01212")
"Hallmark_Fatty_Acid_Metabolism"

#
FA_degradation_core = c("CPT1A","CPT2","SLC25A20","ACSL1","ACADVL","ACADL","ACADM","ACADS","ACADSB","HADHA",
                         "HADHB","HADH","ECHS1","ECI1","ACAA1","ACAA2", "EHHADH", "ACOX1")

FA_metabolism_core = c("ACSL1","ACSL3","ACSL4","ACSL5","CPT1A","CPT2","SLC25A20","ACADVL","ACADM","ACADS",
                       "HADHA","HADHB","ECHS1","ACAA1","ACAA2","EHHADH","ACOX1","SCD","FADS1","FADS2")

Hallmark_FA_metabolism = c("ACAA1","ACAA2","ACADL","ACADM","ACADS","ACADVL","ACOT2","ACOT8","ACOX1","ACSL1","ACSL4","ACSL5",
                           "CPT1A","CPT2","ECHS1","ECI1","ECI2","EHHADH","HADH","HADHB","HSD17B4","SLC22A5","ETFDH","PPARA")

BD_seu = AddModuleScore(BD_seu, features = list(FA_degradation_core), name = 'Fatty_acid_degradation')
BD_seu = AddModuleScore(BD_seu, features = list(FA_metabolism_core), name = 'Fatty_acid_metabolism')
BD_seu = AddModuleScore(BD_seu, features = list(Hallmark_FA_metabolism), name = 'Hallmark_FA_metabolism')

allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")

#
data = FetchData(BD_seu, vars = c('Hallmark_FA_metabolism1', 'celltype'))
ggplot(data, aes(x = celltype, y = Hallmark_FA_metabolism1, fill = celltype)) +
  geom_violin(trim = FALSE, color = NA) +  
  geom_boxplot(width = 0.2, color = "black", fill = "white", lwd = 0.1) +
  scale_fill_manual(values = allcolour) +
  labs(x = 'Cluster', y = 'Fatty acid metabolism_Hallmark Score') +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")
#
VlnPlot(BD_seu, features = 'Fatty_acid_degradation1', group.by = "celltype", assay = "RNA", pt.size = 0, cols = allcolour) +
  geom_boxplot(width = 0.2, col = "black", fill = "white", lwd = 0.1) + 
  labs(x = 'Cluster',
       y = 'Fatty_acid_degradation Score',
       title = "") +   
  coord_flip() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  NoLegend()  









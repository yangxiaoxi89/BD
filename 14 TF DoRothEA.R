#################################### DoRothEA regulon ###############################
library(Seurat)
library(dorothea)
library(decoupleR)
library(viper)
library(dplyr)
library(ggplot2)

if (!requireNamespace("viper", quietly = TRUE))
  BiocManager::install("viper")

if (!requireNamespace("dorothea", quietly = TRUE))
  BiocManager::install("dorothea")

setwd("~/")
##
gc()
# 
base::load("~/annotation.rda")


#
table(BD_seu$celltype)
set.seed(123456)
n_cells = 5000  
cells_subset = sample(colnames(BD_seu), size = min(n_cells, ncol(BD_seu)), replace = FALSE)
seurat_sub = BD_seu[, cells_subset]


# 
expr_mat = GetAssayData(seurat_sub, assay = "RNA", slot = "data")
dim(expr_mat)
celltype = seurat_sub$celltype
cluster_order <- c(
  "Neu_01_CXCL8",
  "Neu_02_TANK",
  "Neu_03_S100A12",
  "Neu_04_IFIT1",
  "Neu_05_STAT1",
  "Neu_06_MMP9",
  "Neu_07_LTF"
)
# 
pseudo_expr = sapply(
  cluster_order,
  function(cl){
    
    cells <- names(
      celltype[celltype == cl]
    )
    
    Matrix::rowMeans(
      expr_mat[,cells,drop=FALSE]
    )
    
  }
)
colnames(pseudo_expr) = cluster_order
colnames(pseudo_expr)
dim(pseudo_expr)

# 
data(dorothea_hs)
dorothea_net = dorothea_hs %>%
  dplyr::filter(
    confidence %in% c("A","B")
  ) %>%
  dplyr::select(
    source = tf,
    target,
    mor
  )
head(dorothea_net)
class(dorothea_hs)

# 
tf_activity <- run_ulm(
  mat = pseudo_expr,
  network = dorothea_net,
  .source = "source",
  .target = "target",
  .mor = "mor"
)
head(tf_activity)
dim(tf_activity)
# 
tf_mat <- tf_activity %>%
  dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = score
  )

# 
tf_mat = as.data.frame(tf_mat)
rownames(tf_mat) = tf_mat$source
tf_mat$source = NULL
tf_mat = as.matrix(tf_mat)
dim(tf_mat)
tf_interest = c("FOXO1", "CEBPD")
#
tf_interest = tf_interest[tf_interest %in% rownames(tf_mat)]
tf_interest
tf_selected = tf_mat[tf_interest,cluster_order]
tf_selected
pheatmap(
  tf_selected,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color=colorRampPalette(c("navy","white","firebrick3"))(100),
  border_color = NA,
  fontsize = 12,
  angle_col = "45"
)






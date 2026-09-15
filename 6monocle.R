######## trajectories #######
setwd("~/")
##
gc()
base::load("~/annotation.rda")
# monocle3
BD_monocle = BD_seu
rm(BD_seu.markers, BD_seu)
#
BD_monocle = BD_monocle[,Idents(BD_monocle) %in% c("Neu_01_CXCL8","Neu_02_TANK","Neu_03_S100A12","Neu_04_IFIT1","Neu_05_STAT1","Neu_06_MMP9")]
gc()
#
table(Idents(BD_monocle))
Idents(BD_monocle) = "group"

bd = BD_monocle[,Idents(BD_monocle) %in% c("BD")]
Idents(bd) = "celltype"
table(Idents(bd))

hc = BD_monocle[,Idents(BD_monocle) %in% c("HC")]
Idents(hc) = "celltype"
table(Idents(hc))

# 
set.seed(123456)
# 
BD_monocle = bd
BD_monocle = hc
BD_monocle = subset(BD_monocle, downsample = 4000)
table(Idents(BD_monocle))
gc()


# cds
expression_matrix = GetAssayData(BD_monocle, assay = 'RNA',slot = 'counts')
cell_metadata = BD_monocle@meta.data
gene_annotation = data.frame(gene_short_name = rownames(BD_monocle))
rownames(gene_annotation) = rownames(BD_monocle)
cds = new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)
#
cds = preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)   
cds = align_cds(cds, num_dim = 50, alignment_group = "orig.ident")  
cds = reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
plot_cells(cds, color_cells_by = "celltype")

#
cds.embed = cds@int_colData$reducedDims$UMAP
int.embed = Embeddings(BD_monocle, reduction = "umap")
int.embed = int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP = int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype")

# 
cds = cluster_cells(cds)
colnames(colData(cds))
plot_cells(cds, show_trajectory_graph = FALSE)  
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)  
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = F)


### learn_graph() function
cds = learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds = order_cells(cds)
plot_cells(cds, color_cells_by = "celltype", label_cell_groups = FALSE,label_leaves = FALSE,label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,label_leaves = FALSE,label_branch_points = FALSE)
plot_cells(cds,    
           label_cell_groups = F,     
           color_cells_by = "pseudotime",   
           label_branch_points = F,   
           label_roots = F,  
           label_leaves = F,   
           graph_label_size = 0,    
           cell_size = 1,   
           trajectory_graph_color = 'black',   
           trajectory_graph_segment_size = 1)  


####### umap #######
allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")
# 
head(colData(cds))
pseudotime = pseudotime(cds) %>% as.data.frame() 
pseudotime$cell = rownames(pseudotime)  
colnames(pseudotime)[1] = "peu"   
#
cell_types = as.data.frame(colData(cds))[ , "celltype", drop = FALSE]   
cell_types$cell = rownames(cell_types)   
merge = merge(pseudotime, cell_types, by = 'cell')
merge = merge[order(merge$peu), ]

## ggplot
p = ggplot(merge, aes(peu, fill = celltype, color = celltype)) +  
  geom_density(alpha = 0.5, size = 1.0) +      
  scale_fill_manual(values = allcolour) +     
  scale_color_manual(values = allcolour) +   
  theme_classic() +   
  scale_y_continuous(expand = c(0,0)) +    
  theme(axis.line = element_blank(),    
        axis.ticks.y = element_blank(),  
        axis.title = element_blank(),    
        axis.text.y = element_blank(),  
        axis.text.x=element_text(colour='black',size=10)) +  
  labs(title = "Pseudotime")
p

# Ridge map
ggplot(merge, aes(x = peu, y = celltype, fill = celltype)) +   
  geom_density_ridges(scale = 1.5, alpha = 0.7) +   
  scale_y_discrete(position = 'right') +   
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size = 15)) +   
  scale_x_continuous(position = 'top') + 
  scale_fill_manual(values = allcolour)




####### pseudotime gene
Track_genes = graph_test(cds, neighbor_graph = "principal_graph", cores = 10)   
save(Track_genes, file = "~/trackgenes_bd.rda")

# 
track_genes = row.names(Track_genes[,c(5,2,3,4,1,6)] %>% dplyr::filter(q_value <= 0.01 & morans_I > 0.2))
write.csv(Track_genes,"~/Trajectory_genes.csv", row.names = F)


### heatmap
mat = pre_pseudotime_matrix(cds_obj = cds, gene_list = track_genes)   
mat = as.data.frame(mat)
head(mat[1:5,1:5])

## kmeans
ck = ClusterGVis::clusterData(obj = mat,
                              cluster.method = "kmeans",
                              cluster.num = 3)

## enrich
enrich = enrichCluster(object = ck,
                       OrgDb = org.Hs.eg.db,
                       type = "BP",
                       pvalueCutoff = 0.05,
                       topn = 10, 
                       seed = 123)
write.csv(enrich,"~/Trajectory_genes enrich.txt", row.names = T)


## 
pdf('monocle1.pdf',height = 5, width = 7, onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = sample(rownames(mat),30,replace = F),
           cluster.order = c(3,2,1),
           border = "white",
           ht.col.list = list(col_range = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),
                              col_color = c("#2166ac", "#4393c3", "#92c5de", 
                                            "#d1e5f0", "white", "#fddbc7", 
                                            "#f4a582", "#d6604d", "#b2182b")))
dev.off()






### expression
plotGeneOverPseudotime <- function(data, gene, output_file2, colors_group) {  
  
  counts = monocle3::exprs(cds)   
  phe = pData(cds)    
  dim(counts)
  #z=as.data.frame(counts[gene_to_cluster,])
  z = as.data.frame(counts)
  ac = phe[,c("group","celltype","Size_Factor")] 
  head(ac)
  A = merge(ac,t(z),by = "row.names") 
  rownames(A) = A$Row.names
  A = A[,-1]
  colnames(A)[3] = "Pseudotime"
  
  
  
  # 
  p2 <- ggplot(data = A, aes(x = Pseudotime, y = .data[[gene]])) +  
    geom_point(aes(color = celltype), alpha = 0.7, size = 0.4) +  
    scale_color_npg() + 
    new_scale_color() +   
    geom_smooth(method = "loess", se = FALSE, aes(color = group)) +  
    scale_color_manual(values = colors_group, aesthetics = "colour", name = "group") +  
    facet_wrap(~group) +    
    theme_bw() +  
    theme(  
      panel.background = element_blank(),  
      panel.border = element_blank(),  
      panel.grid = element_blank(),  
      axis.ticks.length = unit(0.1, "lines"),  
      axis.ticks = element_blank(),  
      axis.line = element_line(color = "black", size = 0.5),  
      axis.title = element_text(size = 10)  
    )  
  
  # 
  w = length(unique(data$group)) * 3  
  
  ggsave(p2, filename = output_file2, width = w, height = 3)  
}  

# 
plotGeneOverPseudotime(cds, "STAT5B",  "gene_Pseudotime_group.pdf",  
                       colors_group = c("BD" = "darkred", "HC" = "#223D6C"))




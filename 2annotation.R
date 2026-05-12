setwd("~/")
base::load("~/clustering.rda")



###### SingleR ######
DefaultAssay(BD_seu) = "RNA"
hpca.se = HumanPrimaryCellAtlasData() 
# DICE = DatabaseImmuneCellExpressionData()
# MID = MonacoImmuneData()
##
BD_forSingleR = GetAssayData(BD_seu, slot="data")
BD_forSingleR
clusters = BD_seu@meta.data$seurat_clusters
# hpca.se
pred.hesc = SingleR(test = BD_forSingleR, 
                    ref = hpca.se, 
                    labels = hpca.se$label.main,
                    method = "cluster", 
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
table(pred.hesc$labels)
##
celltype = data.frame(ClusterID = rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 
BD_seu@meta.data$singleR = celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(BD_seu, reduction = "umap", group.by = "singleR")

# 
DotPlot(BD_seu, features = unique(c("FCGR3B","CSF3R","CXCR2","G0S2","CXCL8","S100A8","S100A9")),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")


### extract clusters
BD_seu = subset(BD_seu, idents = c(0,1,2,3,4,5,6))
DimPlot(BD_seu,reduction = "umap",group.by = "seurat_clusters",label = T)


### neutrophil markers
Mature = c("MME","CXCR4","NAMPT","CXCR2","SOD2","CSF3R","FCGR3B","TOB1","IL18RAP","IRF1","G0S2","SIGLEC10","CX3CR1","CD300E","AIM2","GBP2",
           "FCGR1A","CD274","IFITM1","LY96","ISG20","TRIM25","MEF2C","MPEG1","MX2","IL1R2","GBP5","IFITM3")
Immature = c("MMP9","CTSB","ORM1","PADI4","S100A12","IFI16","OAS2","IFI44","RGL4","LGALS3","PLAC8","CD177","TFF3","KLF2","PTGS2","IFIT1","SLC23A",
             "RSAD2","OAS1","CLU","CD163","MT1X","CD14","CST7","LAIR1","SLC2A3","GAPDH","TYROBP","NAIP","NFKBIZ","LGALS1","ISG15","CD14","CTS7","FCGR3B")
Early_Immature = c("CAMP","LTF","LCN2","CRISP3","OLFM4","MS4A4A","CEACAM1","CDACAM8","MMP8","CEBPE","FCGR3B")
Precursors = c("MPO","ELANE","MKI67","TOP2A","DEFA","AZU1","FCGR3B")   
# 
##
DotPlot(BD_seu, features = unique(Precursors),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DotPlot(BD_seu, features = unique(Early_Immature),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DotPlot(BD_seu, features = unique(Immature),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DotPlot(BD_seu, features = unique(Mature),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
# 







###### find DEGs ######
DefaultAssay(BD_seu) = "RNA"
BD_seu.markers = FindAllMarkers(BD_seu, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
write.table(BD_seu.markers,file="~/BD_seu.markers.txt", sep="\t", row.names=F)
gc()

## 
top10markers = BD_seu.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DotPlot(BD_seu, features = unique(top10markers$gene), group.by = "seurat_clusters") + RotatedAxis() +
  scale_x_discrete("") + scale_y_discrete("") 
## 
DotPlot(BD_seu, features = unique(top10markers$gene), group.by = "celltype") + RotatedAxis() +
  scale_x_discrete("") + scale_y_discrete("") 

p = DotPlot(BD_seu, features = unique(top10markers$gene),
            assay = "RNA") +
  coord_flip() +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +    # 轴标签
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend("Percent Expression")) +    # legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
p


#
DoHeatmap(BD_seu, label = F,
          features = as.character(unique(top10markers$gene)),
          slot = "scale.data",
          group.by = "celltype",  
          assay = "RNA",  
          group.colors = c("#00468B99","#ED000099","#42B54099","#0099B499","#925E9F99","#FDAF9199","#AD002A99")) + 
          scale_fill_gradientn(colors = c('#1A5592','white',"#B83D3D")) + 
          theme(axis.text.y = element_text(size = 10))
         




###### marker ######
new.cluster.ids = c("0" = "Neu_01_CXCL8",                       
                    "1" = "Neu_02_TANK",                      
                    "2" = "Neu_03_S100A12",                       
                    "3" = "Neu_04_IFIT1",                       
                    "4" = "Neu_05_STAT1", 
                    "5" = "Neu_06_MMP9",
                    "6" = "Neu_07_LTF")  
# BD_seu@meta.data$celltype = BD_seu@meta.data$seurat_clusters
# levels(BD_seu@meta.data$celltype) = new.cluster.ids 
BD_seu = RenameIdents(BD_seu, new.cluster.ids)                        
BD_seu$celltype = BD_seu@active.ident
#
DimPlot(BD_seu, label = T, pt.size = 1,group.by = "celltype",label.size = 6,repel = T) + xlim(-15, 15)
DimPlot(BD_seu, label = T, pt.size = 1,group.by = "seurat_clusters",label.size = 6,repel = T) + xlim(-15, 15)


### Umap
umap = BD_seu@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = BD_seu@meta.data$celltype)

head(umap)

allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")

p = ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) + 
    geom_point(size = 0.6, alpha = 1) + 
    scale_color_manual(values = allcolour)
p

p2 = p  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white"))
p2

p3 = p2 +         
  theme(
    legend.title = element_blank(), 
    legend.key = element_rect(fill = 'white'), 
    legend.text = element_text(size = 12),
    legend.key.size = unit(1,'cm')) + 
  guides(color = guide_legend(override.aes = list(size = 4))) 
p3

p4 = p3 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) + 2, yend = min(umap$UMAP_2) ),
               colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 2),
               colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) + 1, y = min(umap$UMAP_2) - 0.5, label = "UMAP_1",
           color = "black", size = 4) + 
  annotate("text", x = min(umap$UMAP_1) - 0.5, y = min(umap$UMAP_2) + 1, label = "UMAP_2",
           color = "black", size = 4, angle = 90) 
p4
#





###### Violin plot of signature genes ######
# colour = c("#D2EBC8","#3C77AF","#7DBFA7","#AECDE1","#EE934E","#D1352B","#9B5B33","#F5CFE4","#B383B9","#8FA4AE","#FCED82","#F5D2A8","#BBDD78","#415284")
# uniqumarker = c("CXCR4","CXCR2","S100P","MMP9","S100A12","STAT1","ISG15","IFIT1","RSAD2","FKBP5","","GBP1","GBP5","CR1","CST7",
#                 "NCF1","CCR1","MX1","MME","GPR141","FCN1","PAG1","PADI4","RFLNB","PTMA","MMP8","CYBA")
uniqumarker = c("CXCR2","CXCR4","CYBA","CXCL8","S100P",
                "S100A12","MME","MMP9","RFLNB","PADI4","GPR141","PAG1",
                "CST7","FCN1","CR1","NCF1",
                "MMP8","S100P","PTMA",
                "STAT1","IFIT1","ISG15","RSAD2","MX1","GBP5","GBP1","CCR1","FKBP5","PAG1")  
VlnPlot(BD_seu, features = uniqumarker, 
        stack = T,  
        sort = F, 
        split.by =  "celltype" ,  
        flip = F,
        ) +
  theme(legend.position = "none",   
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle("") +  
  scale_fill_manual(values = allcolour) +  
  scale_color_manual(values = allcolour)
  

########### 
vln.dat = FetchData(BD_seu,c(uniqumarker,"celltype"))
vln.dat$Cell = rownames(vln.dat)
# 
vln.dat.melt = reshape2::melt(vln.dat, id.vars = c("Cell","celltype"), 
                               measure.vars = uniqumarker,
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(celltype,gene) 
# 
vln.dat.melt$celltype = as.factor(vln.dat.melt$celltype)
vln.dat.melt$gene = as.factor(vln.dat.melt$gene)
vln.dat.melt$Expr = as.numeric(vln.dat.melt$Expr)


## ggplot
ggplot(vln.dat.melt, aes(factor(gene), Expr, fill = celltype)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  facet_grid(rows = vars(celltype), scales = "free", switch = "y")

# 
p1 = ggplot(vln.dat.melt, aes(gene, Expr, color = celltype, fill = celltype)) +   
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +   
  scale_y_continuous(expand = c(0, 0), position ="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +    
  facet_grid(rows = vars(celltype), scales = "free", switch = "y") +
  scale_fill_manual(values = allcolour) +   
  scale_color_manual(values = allcolour) +  
  theme_cowplot(font_size = 12) + 
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(6, 6, 0, 6, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) + ylab("Expression Level")
p1






##### subgroup proportion ######
DimPlot(BD_seu,split.by = "group",ncol = 2, repel = T, raster = F,label = T)
# 
group_table = as.data.frame(table(BD_seu@meta.data$group, BD_seu@meta.data$celltype))
names(group_table) = c("group","cluster","CellNumber")
plot_group = ggplot(BD_seu@meta.data,aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = allcolour) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5, size = 16)) + 
  labs(y = "Percentage") +
  theme_classic() +
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
plot_group



##
gc()
save(BD_seu,BD_seu.markers, file = "~/annotation.rda")
















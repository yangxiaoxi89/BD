#################################### 
##### SCENIC

## package
packageVersion("AUCell")  # (minimum 1.2.4)
# [1] '1.26.0'
packageVersion("RcisTarget")  # (minimum 1.0.2)
# [1] '1.24.0'
packageVersion("GENIE3")  # (minimum 1.2.1)
# [1] '1.26.0'





############### 
base::load("~/annotation.rda")
# meta
dim(BD_seu)
table(BD_seu$orig.ident)
# 
Idents(BD_seu) = "orig.ident"
BD_test = subset(x = BD_seu, downsample = 3000)
table(BD_test$orig.ident)
table(BD_test$celltype)
dim(BD_test)
Idents(BD_test) = "celltype"
write.csv(t(as.matrix(BD_test@assays$RNA@counts)), file = "~/for_scenic_count.csv")
# 
save(BD_test, file = "~/SCENIC.rda")







###############################
base::load("~/SCENIC.rda")
rm(BD_seu.markers,BD_seu)
gc()
# 
# 
loom = open_loom("~/SCENIC/sce_SCENIC.loom")   
# exprMat <- get_dgem(loom)   
# exprMat_log <- log2(exprMat+1)    
# 
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons_incidMat[1:10,1:10]
# 
regulons = regulonsToGeneLists(regulons_incidMat)   
class(regulons)
# 
regulonAUC = get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
# 
regulonAucThresholds = get_regulon_thresholds(loom)
# 
positive_regulons = grep("\\(\\+\\)$", rownames(regulonAUC), value = TRUE)
# 
regulonAUC_pos = regulonAUC[positive_regulons, ]
# 
embeddings = get_embeddings(loom)
close_loom(loom)


# SeuratData （BD_test）
DimPlot(BD_test,reduction = "umap",label=T)    
regulonAUC_pos = regulonAUC_pos[,match(colnames(BD_test),colnames(regulonAUC_pos))]
# 
identical(colnames(regulonAUC_pos), colnames(BD_test))  
# [1] TRUE
dim(regulonAUC_pos)


# Regulon Specificity Score RSS
rss = calcRSS(AUC = getAUC(regulonAUC_pos),    
              cellAnnotation = BD_test@meta.data$celltype)  

rss = na.omit(rss)
write.csv(rss, file = "~/rss_score.csv")







###### rank plot ######
B_rss = as.data.frame(rss)  
colnames(B_rss)
celltype = colnames(B_rss)

# 
rssRanklist = list() 

for(i in 1:length(celltype)) 
{
  data_rank_plot = cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))
  colnames(data_rank_plot) = c("TF", "celltype")
  data_rank_plot = na.omit(data_rank_plot)
  data_rank_plot = data_rank_plot[order(data_rank_plot$celltype, decreasing = TRUE),]
  data_rank_plot$rank = seq(1, nrow(data_rank_plot))
  
  top30 = head(data_rank_plot, 30)

  ## plot
  p = ggplot(data_rank_plot, aes(x = rank, y = celltype)) + 
    geom_point(size = 2, shape = 16, color = "#336699", alpha = 0.4) +   
    geom_point(data = data_rank_plot[1:6,], size = 2, color = '#DC050C') +  
    theme_bw() +   
    theme(axis.title = element_text(colour = 'black', size = 10),   
          axis.text = element_text(colour = 'black', size = 10),    
          axis.text.x = element_blank(),   
          axis.ticks.x = element_blank(),    
          panel.grid = element_blank()) +    
    labs(x = 'Regulons Rank', y = 'Specificity Score', title = celltype[i]) +   
    geom_text_repel(data = data_rank_plot[1:6,],   
                    aes(label = TF), color = "black", size = 3, fontface = "italic",
                    arrow = arrow(ends = "first", length = unit(0.02, "npc")),  
                    box.padding = 0.2, point.padding = 0.3,   
                    segment.color = 'black', segment.size = 0.3,   
                    force = 1, max.iter = 3e3)   
  
  rssRanklist[[i]] = p
}

######## 
combined_plot = plot_grid(plotlist = rssRanklist[c(1,2,4,5,7)], ncol = 5) 
combined_plot






###### heatmap of regulon activity ######
# 
cellInfo = BD_test@meta.data
# 
regulonAUC_pos = regulonAUC_pos[onlyNonDuplicatedExtended(rownames(regulonAUC_pos)),] 
# 
regulonActivity_byCellType = sapply(split(rownames(cellInfo), cellInfo$celltype),   
                                     function(cells) rowMeans(getAUC(regulonAUC_pos)[,cells]))
# 
regulonActivity_byCellType_Scaled = t(scale(t(regulonActivity_byCellType), 
                                             center = T, scale = T))
# ComplexHeatmap
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:50,], name = "Regulon activity")



##### sce.regulons.csv for statistics #####
# 
reg = read.csv("~/SCENIC/sce.regulons.csv")
# 
row1 = as.character(reg[1, ])
row2 = as.character(reg[2, ])
# 
new_names = ifelse(nchar(row1) > 0, row1, row2) 
# 
# new_names <- mapply(function(x, y) {
#   if(nchar(x) > 0 & nchar(y) > 0) {
#     paste0(x, "_", y)
#   } else if(nchar(x) > 0) {
#     x
#   } else {
#     y
#   }
# }, row1, row2, USE.NAMES = FALSE)

#
colnames(reg) = new_names
# 
reg = reg[-c(1,2), ]
colnames(reg)

# function
parse_target_info = function(tf, s) 
{
  pattern = "\\('([^']+)',\\s*([0-9\\.eE+-]+)\\)"
  matches = stringr::str_match_all(s, pattern)[[1]]   
  if(nrow(matches) == 0) return(data.frame(TF = character(), Gene = character(), Score = numeric()))  
  
  
  data.frame(
    TF = rep(tf, nrow(matches)),  
    Gene = matches[,2],   
    Score = as.numeric(matches[,3]),  
    stringsAsFactors = FALSE
  )
}
# 
TF_Genes_score = do.call(rbind, 
                         lapply(seq_len(nrow(reg)), 
                                function(i) 
{
  parse_target_info(reg$TF[i], reg$TargetGenes[i])
}))

head(TF_Genes_score)
TF_Genes_score = unique(TF_Genes_score) 
dim((TF_Genes_score))
range(TF_Genes_score$Score)    
table(TF_Genes_score$Score >= 1)
hist(TF_Genes_score$Score)

###
TF_GeneNum = TF_Genes_score %>%
  ungroup() %>%   
  group_by(TF) %>%   
  summarise(Freq = n(), .groups = "drop")   
TF_GeneNum = as.data.frame(TF_GeneNum)
head(TF_GeneNum)
TF_GeneNum$TF_genenum = paste0(TF_GeneNum$TF, "(", TF_GeneNum$Freq, ")")
head(TF_GeneNum)
write.csv(TF_GeneNum, "~/TF_GeneNum_score.csv", row.names = F)



### 
# 
col_fun = colorRamp2(c(-2,0,2), c("#69AADB", "white", "#D45590"))

# heatmap
p_basic = Heatmap(
  regulonActivity_byCellType_Scaled[1:50,],             
  name = "Z-score",           
  col = col_fun,            
  rect_gp = gpar(           
    col = "white",            
    lwd = 1                 
  ),
  # width = ncol(data_scaled) * unit(6, "mm"),  
  # height = nrow(data_scaled) * unit(6, "mm"), 
  row_names_gp = gpar(fontsize = 12),         
  column_names_gp = gpar(fontsize = 12),        
  column_names_rot = 45,                      
  cluster_rows = TRUE,                    
  cluster_columns = FALSE,                    
  heatmap_legend_param = list(              
    title = "Z-score",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 7)
  )
)
# plot
draw(p_basic)

# 
# round heatmap 
p_round = Heatmap(
  regulonActivity_byCellType_Scaled[1:50,],
  name = "Z-score",
  col  = col_fun,
  rect_gp = gpar(col = NA, fill = NA),  
  # 
  cell_fun = function(j, i, x, y, w, h, fill) {
    v = regulonActivity_byCellType_Scaled[1:50,][i, j]  
    grid.roundrect(         
      x, y,                 
      width  = w * 0.92,    
      height = h * 0.92,  
      r = unit(0.55, "snpc"), 
      gp = gpar(fill = col_fun(v), col = NA)  
    )
  },
  use_raster = FALSE,  
  # width  = ncol(data_scaled) * unit(6, "mm"),
  # height = nrow(data_scaled) * unit(6, "mm"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_names_rot = 90
)
# 
draw(p_round)





################## 
# 
Idents(BD_test)
table(BD_test@active.ident)
AUCmatrix = regulonAUC_pos@assays@data@listData$AUC
AUCmatrix = as.data.frame(AUCmatrix)

# rssRanklist[c(1,2,4,5,7)]
# 
tf.gene = c(rssRanklist[[4]][["data"]]$TF[1:10],
            rssRanklist[[5]][["data"]]$TF[1:10],
            rssRanklist[[1]][["data"]]$TF[1:10],
            rssRanklist[[2]][["data"]]$TF[1:10],
            rssRanklist[[7]][["data"]]$TF[1:10])
tf = AUCmatrix[tf.gene,]
# 
BD_test = AddMetaData(BD_test,metadata = data.frame(t(tf)))
# 
# mycolor = c("#27447C","grey","darkred")
# FeaturePlot(BD_test,features = c("LTF..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("STAT5B..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("CEBPD..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("FOXO1..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("CUX1..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("DBP..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("STAT1..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("STAT2..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("IRF7..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("IRF9..." ),pt.size = 1.0,order = T,cols = mycolor)
# FeaturePlot(BD_test,features = c("ETV7..." ),pt.size = 0.8,order = T,cols = mycolor)
# 


# 
col_grad = c("grey", "#27447C", "#B22222")

p = FeaturePlot(BD_test, 
                features = c("LTF...", "STAT5B...", "STAT1...", "ETV7..."),  
                reduction = "umap",
                pt.size = 0.4, 
                order = TRUE,         
                cols = col_grad,
                combine = FALSE)    

# 
p = lapply(p, function(gg) 
{
  gg + theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 9)) +
    labs(color = "Activity")
})

# 
patchwork::wrap_plots(p, ncol = 4)








############
library(SummarizedExperiment)
seurat.data = BD_test
Idents(seurat.data) = "celltype"
regulonsToPlot = c("STAT5B(+)","CEBPD(+)","FOXO1(+)","CUX1(+)","DBP(+)",
                   "LTF(+)","STAT1(+)","STAT2(+)","IRF7(+)","IRF9(+)",
                   "ETV7(+)","GFI1(+)","CEBPE(+)","CEBPA(+)")
regulonsToPlot %in% row.names(regulonAUC_pos)
seurat.data@meta.data = cbind(seurat.data@meta.data ,
                              t(assay(regulonAUC_pos[regulonsToPlot,])))


# Vis
allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")
p = RidgePlot(seurat.data, 
              features = c("STAT5B(+)","FOXO1(+)","CUX1(+)","CEBPD(+)","LTF(+)","ETV7(+)","STAT1(+)","STAT2(+)","IRF7(+)","IRF9(+)"), 
              ncol = 2,
              cols = allcolour) 
p







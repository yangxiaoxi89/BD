setwd("~/")
##
gc()
base::load("~/annotation.rda")



# enrichement analysis
head(BD_seu.markers)
BD_seu.markers = BD_seu.markers[which(BD_seu.markers[,5]<0.05),]

# 
is.data.frame(BD_seu.markers)
BD_seu.markers$celltype = ""
for (i in 1:nrow(BD_seu.markers)) 
{
  if(BD_seu.markers[i,6]=="0")
  {
    BD_seu.markers[i,8] = "Neu_01_CXCL8"
  }else if(BD_seu.markers[i,6]=="1")
  {
    BD_seu.markers[i,8] = "Neu_02_TANK"
  }else if(BD_seu.markers[i,6]=="2")
  {
    BD_seu.markers[i,8] = "Neu_03_S100A12"
  }else if(BD_seu.markers[i,6]=="3")
  {
    BD_seu.markers[i,8] = "Neu_04_IFIT1"
  }else if(BD_seu.markers[i,6]=="4")
  {
    BD_seu.markers[i,8] = "Neu_05_STAT1"
  }else if(BD_seu.markers[i,6]=="5")
  {
    BD_seu.markers[i,8] = "Neu_06_MMP9"
  }else
  {
    BD_seu.markers[i,8] = "Neu_07_LTF"
  }
  print(i)
}

# 
gene_id = bitr(BD_seu.markers$gene, fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

# 
mix = merge(BD_seu.markers, gene_id, by.x ="gene", by.y = "SYMBOL")
head(mix)

# 
deg_gene_clusters = data.frame(ENTREZID = mix$ENTREZID,
                               celltype = mix$celltype)
head(deg_gene_clusters)
table(deg_gene_clusters)
# list
list_deg_gene_clusters = split(deg_gene_clusters$ENTREZID, deg_gene_clusters$celltype)





############
formula_go = compareCluster(
  ENTREZID~celltype,
  data = deg_gene_clusters,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",      ### One of "BP", "MF", and "CC" or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)



formula_kegg = compareCluster(
  ENTREZID~celltype,
  data = deg_gene_clusters,
  fun="enrichKEGG",
  keyType = "kegg",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)

##
go = as.data.frame(formula_go)
dim(go)

kegg = as.data.frame(formula_kegg)
dim(kegg) 
go10 = go %>% group_by(celltype) %>% top_n(-10, qvalue)
kegg10 = kegg %>% group_by(celltype) %>% top_n(-10, qvalue)

write.table(go, file = "~/goenrich.txt", sep="\t", row.names=F)
write.table(kegg, file = "~/keggenrich.txt", sep="\t", row.names=F)





##### Bubble plot #####
getwd()
# 
keggenrich = read.table(file = "~/keggenrich.txt",header = T,row.names = NULL,check.names = F,sep="\t",stringsAsFactors = F)
#
padj = unlist(keggenrich[,13])
summary(padj)
#
keggenrich$p = '' 
keggenrich$p[which(keggenrich$p.adjust >= 5.880e-06 & keggenrich$p.adjust < 6.263e-04)] = '< 6.263e-04'
keggenrich$p[which(keggenrich$p.adjust > 6.263e-04 & keggenrich$p.adjust < 1.132e-02)] = '< 1.132e-02'
keggenrich$p[which(keggenrich$p.adjust > 1.132e-02 & keggenrich$p.adjust < 1.323e-02)] = '< 1.323e-02'
keggenrich$p[which(keggenrich$p.adjust > 1.323e-02 & keggenrich$p.adjust < 2.619e-02)] = '< 2.619e-02'
keggenrich$p[which(keggenrich$p.adjust > 2.619e-02 & keggenrich$p.adjust <= 4.818e-02)] = '<= 4.818e-02'
# 
generatio = unlist(keggenrich[,7])
class(generatio)
generatio = sapply(generatio, function(x) eval(parse(text = x)))
print(generatio)
keggenrich[,7] = generatio
# 
keggenrich[,2] = as.character(keggenrich[,2])
keggenrich[,2]

# 
ggplot(data = keggenrich, aes(x = cluster,y = Description))+
  geom_tile(aes(fill = p),
            color = "white")+
  geom_point(aes(size = GeneRatio), color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        #legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(brewer.pal(9,"Blues")[c(3,5,7,8,9)]))+
  labs(fill = "p.adjust")
# brewer.pal(9,"Blues")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"



### 
getwd()
goenrich = read.table(file = "~/goenrich.txt",header = T,row.names = NULL,check.names = F,sep="\t",stringsAsFactors = F)


# 
generatio = unlist(goenrich[,4])
class(generatio)
generatio = sapply(generatio, function(x) eval(parse(text = x)))  
print(generatio)
goenrich[,4] = generatio

# 
goenrich$Description = fct_inorder(goenrich$Description)


# heatmap
allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")
summary(goenrich$GeneRatio) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03448 0.07692 0.09249 0.11834 0.15254 0.37500 

ggplot(goenrich, aes(x = GeneRatio, y = Description, colour = Cluster)) +
  scale_colour_manual(name = "Cluster", values = allcolour) + 
  scale_shape_manual(name = "", values = c(1,1,25,25)) +
  facet_wrap(~Cluster, ncol = 7) +
  geom_point(size = 4, stroke = 0, shape = 16) +
  theme_bw() + xlim(0,0.4) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 0),  
        axis.text.y = element_text(size = 11, color = "black", face = "bold"),  
        axis.title = element_text(size = 10, face = "bold"), text = element_text(size = 15)) +  
  geom_vline(aes(xintercept = 0.12), colour="black", size = 0.1, linetype = "dashed") + 
  xlab("Enrichment Ratio") + ylab("GO term") +  
  # theme(axis.text.y = element_text(size = 10)) +  
  scale_y_discrete(position = "right") +
  theme(legend.direction = "horizontal", legend.position = "top") 














##########################################################
###### GSVA ######
gc()
Idents(BD_seu) = "celltype"
# 
expr = AverageExpression(BD_seu, assays = "RNA", slot = "data")[[1]]
# 
expr = expr[rowSums(expr)>0,]
# 
expr = as.matrix(expr)

##
# 
genesets = msigdbr(species = "Homo sapiens", category = "H") 
genesets = subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets = split(genesets$gene_symbol, genesets$gs_name)
params = gsvaParam(exprData = expr,
                   geneSets = genesets,
                   minSize = 10)
gsva_res = gsva(params)

## 
gsva.df = data.frame(Genesets = rownames(gsva_res), gsva_res, check.names = F)

#####
gsva.df = gsva.df[,-1] 
gsub("_"," ",substr(rownames(gsva.df)[1], 10, nchar(rownames(gsva.df)[1])))
rownames(gsva.df) = gsub("_"," ",substr(rownames(gsva.df), 10, nchar(rownames(gsva.df))))
# 
rownames(gsva.df) = str_to_sentence(str_to_lower(rownames(gsva.df))) 
rows_to_keep = unique(unlist(lapply(gsva.df, function(col) {
  order(col, decreasing = TRUE)[1:3]
})))


filtered_gsva.df = gsva.df[rows_to_keep, ]
# 
write.table(filtered_gsva.df, file = "~/filtered_gsva_hallmark.txt", sep="\t", row.names = T)

## complexheatmap
allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")

df = data.frame(colnames(filtered_gsva.df))
colnames(df) = 'cell_type'
top_anno = HeatmapAnnotation(df = df,   
                             border = F,  
                             show_annotation_name = F,
                             gp = gpar(col = NA), 
                             col = list(cell_type = c("Neu_01_CXCL8" = "#27447C",
                                                      "Neu_02_TANK" = "#73ABCF",
                                                      "Neu_03_S100A12" = "#C72228",
                                                      "Neu_04_IFIT1" = "#9EAAD1",
                                                      "Neu_05_STAT1" = "#168676",
                                                      "Neu_06_MMP9" = "#F3B169",
                                                      "Neu_07_LTF" = "#B88640")))
# 
gsva_exp = t(scale(t(filtered_gsva.df), scale = T, center = T)) 
gsva_exp[which(gsva_exp > 2)] = 2
# 
Heatmap(gsva_exp,
        #name = "Z-score",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        # show_row_dend = T,   
        # show_column_dend = T,
        heatmap_legend_param = list(title = "GSVA score"),  
        col = colorRampPalette(brewer.pal(11, "PiYG"))(11), 
        border = NA,   
        rect_gp = gpar(col = "white", lwd = 1), 
        top_annotation = top_anno   
        ) 



#############################   volcano map   ##############################
BD_seu.markers1 = FindAllMarkers(BD_seu,
                                 only.pos = F,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)  
colour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")

## 
p = jjVolcano(diffData = BD_seu.markers1,
          log2FC.cutoff = 0.25, 
          size  = 3.5, 
          fontface = 'italic', 
          #aesCol = c('purple','orange'),
          tile.col = colour, 
          #col.type = "adjustP", 
          topGeneN = 5, 
          legend.position = c(0.1,0.9)
)








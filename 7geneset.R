#################################### 
base::load("~/annotation.rda")
azurophil = c("HEXA","PRSS57","PRTN3","CTSG","MPO","CTSC")   
specific = c("LCN2","CYBA","CYBB","NCF1","NCF4","LTF","CAMP")  
gelatinase = c("CEACAM1","MMP8","MMP9","ITGAM","SLC11A1")   
secretory = c("CD63","MME","ITGAM","FUT4","CR2","CYBB","MMP25","SLC11A2","FPR1","SCAMP1","VAMP2","STXBP4","CD93","CR1")    
NETs = c("PADI4","CRISPLD2","DYSF","CAT","S100A8","S100A9","ACTB","ACTG1","ACTN1","LCP1","LYZ","MYH9","S100A12","TKT")
geneset = c(azurophil,specific,gelatinase,secretory,NETs)


### 
gene_cell_exp = AverageExpression(BD_seu,     
                                  features = unique(geneset) ,
                                  group.by = 'celltype',
                                  slot = 'data')  
gene_cell_exp = as.data.frame(gene_cell_exp$RNA)


### complexheatmap
# allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")
df = data.frame(colnames(gene_cell_exp))
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

# scale data
marker_exp = t(scale(t(gene_cell_exp), scale = T, center = T))   

### heatmap
# 1）
row_groups = data.frame(
  Gene = c(azurophil, specific, gelatinase, secretory, NETs),
  Group = c(rep("Azurophil granules", length(azurophil)),
            rep("Specific granules", length(specific)),
            rep("Gelatinase granules", length(gelatinase)),
            rep("Secretory vesicles", length(secretory)),
            rep("NETs formation", length(NETs))
))

# 2）
row_order = row_groups$Gene
marker_exp = marker_exp[row_order, ]  

# 3）
group_colors = c("Azurophil granules" = "#C1E6F3",
                 "Specific granules" = "#CCC9E6",
                 "Gelatinase granules" = "#E0D4CA",
                 "Secretory vesicles" = "#C5DEBA",
                 "NETs formation" = "#EFA9AE")

# 4）
row_anno = rowAnnotation( 
  df = data.frame(Group = row_groups$Group),  
  col = list(Group = group_colors),   
  show_annotation_name = FALSE,   
  gp = gpar(col = NA))  

# 5）
Heatmap(marker_exp,
        #name = "Z-score",
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = T,
        show_row_dend = F,   
        show_column_dend = F,
        #column_title = NULL,
        heatmap_legend_param = list(title = "Z-score"), 
        col = colorRampPalette(c("#11427C", "white", "#C31E1F"))(100),  
        border = NA,  
        rect_gp = gpar(col = "white", lwd = 1), 
        row_names_gp = gpar(fontsize = 12), 
        #column_names_gp = gpar(fontsize = 5),   
        top_annotation = top_anno,   
        left_annotation = row_anno,  
        #row_title = NULL,  
        row_split = row_groups$Group)  




## 
FeaturePlot(object = BD_seu, features = "MPO", pt.size = 1, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) + 
  DarkTheme() + 
  theme(text=element_text(size = 14))+ 
  theme(text=element_text(face = "bold"))+
  theme(legend.text=element_text(size = 8))

marker_sign = c("ELANE","PRTN3","CAMP","LTF","MMP9","MMP8","MME","CR1")
# FeaturePlot(BD_seu,features = marker_sign)













###########################################
# cytokines、MHC molecules、immunophenotypes signatures
base::load("~/annotation.rda")

cytokines = c("CCL5","CXCL8","CCL20","TNFSF13B","CD274","CCR1","CXCR2","CXCR1","CXCR4","CCL3","CXCL2","CCL4","TNFRSF14","CXCL3")   
MHC = c("HLA-DPB1","HLA-DPA1","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-C","HLA-B","HLA-A","HLA-E")   
interferon = c("IRF8","DDIT3","TRAF3","PDE4B","IL1B","ISG15","IRF7","IFIH1","STAT1","CD274")  


# average
gene_cell_exp = AverageExpression(BD_seu,     
                                  features = unique(MHC) ,
                                  group.by = 'celltype',
                                  slot = 'data')
gene_cell_exp = as.data.frame(gene_cell_exp$RNA)


# scale data
gene_cell_tran = t(gene_cell_exp)
gene_cell_trn = gene_cell_tran %>% scale(center = TRUE) %>% as.data.frame() 
gene_cell_tran$celltype = rownames(gene_cell_tran)

# 
gene_cell_long = gene_cell_tran %>% 
  pivot_longer(cols = -celltype, names_to = "gene", values_to = "expression")   
# 
gene_cell_long$group = gene_cell_long$celltype

# ggplot
ggplot(gene_cell_long, aes(x = celltype, y = gene, color = expression, size = abs(expression))) +
  geom_point() +
  scale_color_gradient2(low = "#003366", mid = "white", high = "#990033") +
  scale_size_continuous(range = c(3, 8),
                        name = "Value",
                        breaks = c(min(abs(gene_cell_long$expression)), 
                                   max(abs(gene_cell_long$expression))),
                        labels = c("Low", "High")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        strip.text.x = element_blank()) +
  labs(x = "Cell Type", y = "Gene", color = "Expression", size = "Expression", title = "Interferon-related") +
  # 
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  guides(scolor = guide_colorbar(order = 1),
         size = guide_legend(order = 2)) 














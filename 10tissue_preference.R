#################################### 
##### Ro/e


base::load("~/annotation.rda")
# 
rm(BD_seu.markers)
dim(BD_seu)
data = BD_seu@meta.data
data = data[data$celltype %in% c("Neu_01_CXCL8","Neu_02_TANK","Neu_03_S100A12","Neu_04_IFIT1","Neu_05_STAT1","Neu_06_MMP9"), ]
data = droplevels(data) 
data = data[,c(1,6,13)]


# 
R_oe = calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "celltype",  
                     colname.patient = "orig.ident",  
                     colname.tissue = "group",      
                     method = "chisq",    # "chisq", "fisher", and "freq"
                     min.rowSum = 0) 
R_oe

write.table(R_oe, file = "~/R_oe.txt", sep="\t", row.names = T)


# 
# Heatmap(as.matrix(R_oe),
#         col = colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)),
#                          c("blue", "white", "red")),
#         show_heatmap_legend = TRUE,
#         cluster_rows = TRUE,
#         cluster_columns = TRUE)


col_fun = colorRamp2(c(0.6, 0.95, 1.15, 1.6), c("grey", "#FAFAD2",  "#FF8C00", "#8B0000"))   
Heatmap(as.matrix(R_oe),
        show_heatmap_legend = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 18),  
        column_names_gp = gpar(fontsize = 18), 
        column_names_rot = 45, 
        heatmap_legend_param = list(
          title = "Ro/e Index",  
          at = seq(0.5, 2, by = 0.5), 
          labels = seq(0.5, 2, by = 0.5) 
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- R_oe[i, j]
          if (value > 1.5) {
            label <- "+++"
          } else if (value > 1.0) {
            label <- "++"
          } else if (value >= 0.8) {
            label <- "+"
          } else if (value > 0.5) {
            label <- "+/-"
          } else if (value == 0) {
            label <- "-"
          } else {
            label <- NA  
          }
          grid.text(label, x, y, gp = gpar(fontsize = 30 , col = "black"))  
        }
)












setwd("~/")
##
base::load("~/annotation.rda")


##                  
Specific_granules = c("LCN2", "CYBA", "CYBB", "NCF1", "NCF4", "LTF", "CAMP")   
Gelatinase_granules = c("CEACAM1", "MMP8", "MMP9", "ITGAM", "SLC11A1")   
Secretory_vesicles = c("CD63", "MME", "ITGAM", "FUT4", "CR2", "CYBB", 
                       "MMP25", "SLC11A2", "FPR1", "SCAMP1", "VAMP2", "STXBP4", 
                       "CD93", "CR1L")   
NETs_formation = c("PADI4", "CRISPLD2", "DYSF", "CAT", "S100A8", 
                   "S100A9", "ACTB", "ACTG1", "ACTN1", "LCP1", 
                   "LYZ", "MYH9", "S100A12", "TKT") 
Phagocytosis = c("ABCA1", "ADGRB1", "AIF1", "ARHGAP12", "ARHGAP25", 
                 "BECN1", "BIN2", "CDC42", "CLCN3", "ELMO1",
                 "FCER1G", "GSN", "GULP1",
                 "IGLL1", "ITGAM", "ITGB2", "MARCO", "MEGF10",
                 "MFGE8", "MSR1", "MYH9", "RAC1", "RAC3",
                 "RHOBTB1", "RHOBTB2", "SH3BP1", "SIRPA", "THBS1",
                 "TREM2", "TREML4", "VAMP7", "XKR4", "XKR6",
                 "XKR7", "XKR8", "XKR9")  
Chemotaxis = c("C5AR1", "CCL2", "CCL20", 
               "CCL22", "CCL25", 
               "CCL3", "CCL4", "CCL5", "CKLF", 
               "CSF3R", "CX3CL1", "CXADR", "CXCL1", "CXCL10", "CXCL13", 
               "CXCL2", "CXCL3", "CXCL5", "CXCR1", "CXCR2", "EDN3", 
               "FCER1G", "GBF1", "IFNG", "IL17B", "IL1B", 
               "IL1RN", "ITGA1", "ITGA9", "ITGAM", "ITGB2", "LGALS3", 
               "NCKAP1L", "PDE4B", "PDE4D", "PF4", "PPBP", "PREX1", 
               "PRKCA", "S100A8", "S100A9", "SLC37A4", "SPP1", "SYK", "TGFB2", 
               "TREM1", "VAV1", "VAV3", "XCL1")   
Activation = c("ABR", "PTAFR", "TYROBP", "STX11", "BCR", "SYK", "VAMP7", 
               "DNASE1", "DNASE1L3", "FCER1G", "PRAM1", "ITGB2", "ITGAM", 
               "CD177", "ANXA3", "MYO1F", "CEACAM8", "PRTN3") 
NADPH_oxidase	= c("CYBB", "CYBA", "RAC2", "RAC1", "NCF2", "NCF1", "NCF4") 
Glycolysis = c("HK1", "HK2", "HKDC1", "PFKL", "PFKM", 
               "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", 
               "PGAM2", "ENO1", "ENO2", "PKLR", "PKM")  
ROS_formation = c("CYBA","CYBB","CYP1A1","CYP1A2","CYP1B1","DDAH1","DUOX1","DUOX2","GBF1","HSP90AA1",
                  "MPO","NCF1","NOS1","NOS3","P2RX4","RORA","SLC7A2","SOD1","SOD2","SPR") 

#######
BD_seu = AddModuleScore(BD_seu,features = list(Specific_granules), name = 'Specific_granules')
BD_seu = AddModuleScore(BD_seu,features = list(Gelatinase_granules), name = 'Gelatinase_granules')
BD_seu = AddModuleScore(BD_seu,features = list(Secretory_vesicles), name = 'Secretory_vesicles')
BD_seu = AddModuleScore(BD_seu,features = list(NETs_formation), name = 'NETs_formation')
BD_seu = AddModuleScore(BD_seu,features = list(Phagocytosis), name = 'Phagocytosis')
BD_seu = AddModuleScore(BD_seu,features = list(Chemotaxis), name = 'Chemotaxis')
BD_seu = AddModuleScore(BD_seu,features = list(Activation), name = 'Activation')
BD_seu = AddModuleScore(BD_seu,features = list(NADPH_oxidase), name = 'NADPH_oxidase')
BD_seu = AddModuleScore(BD_seu,features = list(Glycolysis), name = 'Glycolysis')
BD_seu = AddModuleScore(BD_seu,features = list(ROS_formation), name = 'ROS_formation')




## 
# mycolor = c("lightblue","grey","darkred")
allcolour = c("#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640")

# 
# FeaturePlot(BD_seu,features = 'Maturation1',pt.size = 0.5,order = T,cols = mycolor)
VlnPlot(BD_seu, features = 'Chemotaxis1', group.by = "celltype", assay = "RNA", pt.size = 0, cols = allcolour) +
  geom_boxplot(width = 0.2, col = "black", fill = "white", lwd = 0.1) +
  labs(x = 'Cluster',
       y = 'Chemotaxis Score',
       title = "") +  
  coord_flip() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  NoLegend()  



###
# data <- FetchData(BD_seu, vars = c('Specific_granules1', 'celltype'))
# # 
# ggplot(data, aes(x = celltype, y = Specific_granules1, fill = celltype)) +
#   geom_violin(trim = FALSE, color = NA) + 
#   geom_boxplot(width = 0.2, color = "black", fill = "white", lwd = 0.1) +
#   scale_fill_manual(values = allcolour) +
#   labs(x = 'Cluster', y = 'Specific Granules Score') +
#   coord_flip() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
#         legend.position = "none")

  
  






                                                                                                                                                                                                                                                                                                                                                                                       
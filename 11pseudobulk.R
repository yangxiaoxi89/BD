setwd("~/")
# 
base::load("~/annotation.rda")
table(BD_seu$group)
table(BD_seu$orig.ident)
table(BD_seu$celltype)


######################### pseu-bulk ####################################
### raw counts
counts = GetAssayData(BD_seu, assay = "RNA", slot = "counts")
### donor-celltype
meta = BD_seu@meta.data
meta$donor = meta$orig.ident
meta$group = meta$group
meta$celltype = meta$celltype
table(meta$donor, meta$group)
table(meta$donor, meta$celltype)

### donor*celltype pseudobulk
meta$pb_id = paste(meta$donor, meta$celltype, sep = "_")
pb_ids = unique(meta$pb_id)  
pb_ids
pb_counts = sapply(
  pb_ids,
  function(x) 
  {
    cells = rownames(meta)[meta$pb_id == x]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }
)
pb_counts = as.matrix(pb_counts)
dim(pb_counts) 


### pseudobulk metadata
pb_meta = data.frame(pb_id = colnames(pb_counts))
pb_meta$donor = sapply(strsplit(pb_meta$pb_id, "_"), `[`, 1)   
pb_meta = unique(meta[, c("pb_id", "donor", "group", "celltype")])    
rownames(pb_meta) = pb_meta$pb_id
pb_meta = pb_meta[colnames(pb_counts), ]
all(rownames(pb_meta) == colnames(pb_counts))
# 
table(pb_meta$celltype, pb_meta$group)
table(pb_meta$donor, pb_meta$celltype)


### 
library(DESeq2)
celltypes_use = unique(pb_meta$celltype)
DE_results = list()
#
for (ct in celltypes_use) {
  
  message("Running: ", ct)
  idx = pb_meta$celltype == ct
  count_sub = pb_counts[, idx, drop = FALSE]
  meta_sub = pb_meta[idx, , drop = FALSE]
  rownames(meta_sub) = meta_sub$pb_id
  meta_sub$group = factor(meta_sub$group, levels = c("HC", "BD"))
  
  if (
    sum(meta_sub$group == "BD") < 2 ||
    sum(meta_sub$group == "HC") < 2) 
  {
    message("Skipping ", ct, ": insufficient donors")
    next
  }
  
  dds = DESeqDataSetFromMatrix(
    countData = round(count_sub),
    colData = meta_sub,
    design = ~ group)
  
  keep = rowSums(counts(dds) >= 10) >= 2
  dds = dds[keep, ]
  
  dds = DESeq(dds)
  
  res = results(
    dds,
    contrast = c("group", "BD", "HC")
  )
  
  res = as.data.frame(res)
  res$gene = rownames(res)
  res$celltype = ct
  res = res[order(res$padj, na.last = T), ]

  DE_results[[ct]] = res
 
   
}
DE_results_all = do.call(
  rbind,
  DE_results
)
rownames(DE_results_all) = NULL
DE_results_all = DE_results_all[which(DE_results_all[,2]>=0), ]
write.csv(DE_results_all, file = "~/pseudpbulk.csv",row.names = FALSE)




######################### FindMarkers ####################################
table(BD_seu$group, BD_seu$orig.ident)

### 
library(MAST)
celltypes_use = unique(BD_seu@meta.data$celltype)
# 
marker_results = list()
#
for (ct in celltypes_use) 
{
  
  message("Running for: ", ct)
  seu_sub = subset(BD_seu, subset = celltype == ct)
  seu_sub$group = factor(seu_sub$group, levels = c("HC", "BD"))
  
  res = FindMarkers(
    object = seu_sub,
    ident.1 = "BD",
    ident.2 = "HC",
    group.by = "group",
    test.use = "MAST",
    latent.vars = "orig.ident"
  )
  
  res$gene = rownames(res)
  res$celltype = ct
  
  marker_results[[ct]] = res
}
marker_results_all = do.call(
  rbind,
  marker_results
)
rownames(marker_results_all) = NULL
marker_results_all = marker_results_all[which(marker_results_all[,2]>=0), ]
write.csv(marker_results_all, file = "~/findmarkers.csv", row.names = F)
save(marker_results,DE_results, file = "~/pseudobulk.rda")



####################### compare ###############################
bulk_sig = rownames(DE_results[[7]])[DE_results[[7]]$padj < 0.05]   
cell_sig = rownames((marker_results)[[7]])[marker_results[[7]]$p_val_adj < 0.05]
# 
bulk_sig = bulk_sig[!is.na(bulk_sig) & bulk_sig != ""]
cell_sig = cell_sig[!is.na(cell_sig) & cell_sig != ""]


### 
library(VennDiagram)
#
venn.plot = venn.diagram(
  x = list(
    `Pseudo-bulk DEGs` = bulk_sig4,
    `Single-cell DEGs` = cell_sig
  ),
  category.names = c(
    "Pseudo-bulk DEGs",
    "Single-cell DEGs"
  ),
  filename = NULL,
  fill = c("#f93", "#56B"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  margin = 0.1
)
grid.newpage()
grid.draw(venn.plot)


### 
overlap_genes = intersect(bulk_sig,cell_sig)
#  
df_cor = data.frame(
  gene = overlap_genes,
  log2FC_bulk = DE_results[[1]][overlap_genes, "log2FoldChange"],
  log2FC_cell = marker_results[[1]][overlap_genes, "avg_log2FC"]
)
# 
cor_test = cor.test(df_cor$log2FC_bulk, df_cor$log2FC_cell, method = "pearson")
cat("Pearson r =", round(cor_test$estimate, 3),
    "， p-value =", format.pval(cor_test$p.value, digits = 3), "\n")

####
library(ggplot2)
library(ggpubr) 

p = ggplot(df_cor, aes(x = log2FC_bulk, y = log2FC_cell)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red") + 
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
  labs(x = "Log2 Fold Change (Pseudo-bulk)",
       y = "Log2 Fold Change (Single-cell)",
       title = paste("Correlation of Log2FC for Overlapping DEGs (n=", nrow(df_cor), ")", sep="")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3)

p


### consistency
genes_use = c("STAT1", "IFIT1")

## 
seu_sub = subset(BD_seu,subset = celltype == "Neu_04_IFIT1")
# 
p <- VlnPlot(
  seu_sub,
  features = genes_use,
  group.by = "group",
  pt.size = 0,
  same.y = FALSE
) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    linewidth = 0.5
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "none"
  ) +
  labs(
    y = "Expression level"
  )

p

# 
library(edgeR)
# 
ct = "Neu_04_IFIT1"
idx = pb_meta$celltype == ct
count_sub = pb_counts[, idx, drop = FALSE]
meta_sub = pb_meta[idx, , drop = FALSE]
rownames(meta_sub) = meta_sub$pb_id
count_sub = count_sub[, rownames(meta_sub), drop = FALSE]
y = DGEList(counts = round(count_sub))
y = calcNormFactors(y)
logCPM = cpm(y, log = TRUE, prior.count = 1)
pb_plot = as.data.frame(t(logCPM[genes_use, , drop = FALSE]))
pb_plot$pb_id = rownames(pb_plot)
pb_plot <- pb_plot %>%
  left_join(
    meta_sub[, c("pb_id", "group")],
    by = "pb_id"
  ) %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  )

pb_mean <- pb_plot %>%
  group_by(gene, group) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  pb_plot,
  aes(
    x = group,
    y = expression
  )
) +
  
  geom_jitter(
    width = 0.08,
    size = 3,
    alpha = 0.8
  ) +
  
  geom_line(
    data = pb_mean,
    aes(
      x = group,
      y = mean_expression,
      group = gene
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = pb_mean,
    aes(
      x = group,
      y = mean_expression
    ),
    size = 4
  ) +
  
  facet_wrap(
    ~ gene,
    scales = "free_y"
  ) +
  
  theme_classic() +
  
  labs(
    x = NULL,
    y = "log-CPM"
  )







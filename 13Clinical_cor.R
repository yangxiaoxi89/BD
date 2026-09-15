# 
clinc_inf = read.csv(file = "~/R clin_info.csv",header = T,stringsAsFactors = F)

#
colnames(clinc_inf)
# [1] "Sample.ID"                              "pct_MMP9"                               "pct_S100A12"                           
# [4] "pct_IFIT1"                              "pct_STAT1"                              "Disease.activity"                      
# [7] "Beh.et.s.Disease.Current.Activity.Form" "Medication.exposure"                    "Disease.duration"                      
# [10] "Organ.involvement"                      "Systemic.inflammatory.status"  

# 
clinc_inf$Disease.activity_num = ifelse(clinc_inf$Disease.activity == "Active", 1, 0)   
clinc_inf$BDCAF = as.numeric(clinc_inf$Beh.et.s.Disease.Current.Activity.Form)   
clinc_inf$medicine = ifelse(clinc_inf$Used.medication == "Methylprednisolone", 1, 0)   

# activity
cor_mmp9 = cor.test(clinc_inf$pct_MMP9, clinc_inf$Disease.activity_num, method = "spearman")
cor_s100a12 = cor.test(clinc_inf$pct_S100A12, clinc_inf$Disease.activity_num, method = "spearman")
print(cor_mmp9)
print(cor_s100a12)

# bdcaf
cor_mmp9_bdcaf = cor.test(clinc_inf$pct_MMP9, clinc_inf$BDCAF, method = "spearman")
cor_s100a12_bdcaf = cor.test(clinc_inf$pct_S100A12, clinc_inf$BDCAF, method = "spearman")
print(cor_mmp9_bdcaf)
print(cor_s100a12_bdcaf)



# 
# calculate correlations


cor_results <- data.frame(
  Comparison = c(
    "MMP9_DiseaseActivity",
    "S100A12_DiseaseActivity",
    "MMP9_BDCAF",
    "S100A12_BDCAF"
  ),
  rho = NA,
  p_value = NA
)

# Disease activity
res1 <- cor.test(
  clinc_inf$pct_MMP9,
  clinc_inf$Disease.activity_num,
  method = "spearman",
  exact = FALSE
)

res2 <- cor.test(
  clinc_inf$pct_S100A12,
  clinc_inf$Disease.activity_num,
  method = "spearman",
  exact = FALSE
)

# BDCAF
res3 <- cor.test(
  clinc_inf$pct_MMP9,
  clinc_inf$BDCAF,
  method = "spearman",
  exact = FALSE
)

res4 <- cor.test(
  clinc_inf$pct_S100A12,
  clinc_inf$BDCAF,
  method = "spearman",
  exact = FALSE
)

cor_results$rho <- c(
  res1$estimate,
  res2$estimate,
  res3$estimate,
  res4$estimate
)

cor_results$p_value <- c(
  res1$p.value,
  res2$p.value,
  res3$p.value,
  res4$p.value
)

# BH multiple-testing correction
cor_results$adjusted_p <- p.adjust(
  cor_results$p_value,
  method = "BH" 
)

cor_results

#####
# activity
plot_activity <- function(data, gene, rho, p, adjusted_p) {
  
  x <- data[[paste0("pct_", gene)]]
  y <- data$Disease.activity_num
  
  p_label <- ifelse(
    p < 0.001,
    "P < 0.001",
    paste0("P = ", sprintf("%.3f", p))
  )
  
  adj_p_label <- ifelse(
    adjusted_p < 0.001,
    "FDR < 0.001",
    paste0("FDR = ", sprintf("%.3f", adjusted_p))
  )
  
  ggplot(data, aes(
    x = factor(
      Disease.activity_num,
      levels = c(0, 1),
      labels = c("Inactive", "Active")
    ),
    y = x
  )) +
    
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      linewidth = 0.5
    ) +
    
    geom_jitter(
      width = 0.12,
      size = 2.5,
      alpha = 0.85
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = max(x, na.rm = TRUE) * 0.95,
      label = paste0(
        "Spearman \u03c1 = ", round(rho, 2),
        "\n", p_label,
        "\n", adj_p_label
      ),
      hjust = 0.5,
      size = 3.5
    ) +
    
    labs(
      x = "Disease activity",
      y = paste0(gene, "+ neutrophils (%)")
    ) +
    
    theme_classic() +
    
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
}


p_MMP9_activity <- plot_activity(
  clinc_inf,
  "MMP9",
  cor_results$rho[1],
  cor_results$p_value[1],
  cor_results$adjusted_p[1]
)

p_S100A12_activity <- plot_activity(
  clinc_inf,
  "S100A12",
  cor_results$rho[2],
  cor_results$p_value[2],
  cor_results$adjusted_p[2]
)

# bdcaf
plot_bdcaf <- function(data, gene, rho, p, adjusted_p) {
  
  x <- data[[paste0("pct_", gene)]]
  y <- data$BDCAF
  
  p_label <- ifelse(
    p < 0.001,
    "P < 0.001",
    paste0("P = ", sprintf("%.3f", p))
  )
  
  adj_p_label <- ifelse(
    adjusted_p < 0.001,
    "FDR < 0.001",
    paste0("FDR = ", sprintf("%.3f", adjusted_p))
  )
  
  ggplot(data, aes(
    x = x,
    y = y
  )) +
    
    geom_point(size = 3) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.7
    ) +
    
    geom_text(
      aes(label = Sample.ID),
      vjust = -0.8,
      size = 3
    ) +
    
    annotate(
      "text",
      x = max(x, na.rm = TRUE) * 0.75,
      y = max(y, na.rm = TRUE) * 0.90,
      label = paste0(
        "Spearman \u03c1 = ", round(rho, 2),
        "\n", p_label,
        "\n", adj_p_label
      ),
      size = 3.5
    ) +
    
    labs(
      x = paste0(gene, "+ neutrophils (%)"),
      y = "BDCAF score"
    ) +
    
    theme_classic() +
    
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
}

p_MMP9_bdcaf <- plot_bdcaf(
  clinc_inf,
  "MMP9",
  cor_results$rho[3],
  cor_results$p_value[3],
  cor_results$adjusted_p[3]
)

p_S100A12_bdcaf <- plot_bdcaf(
  clinc_inf,
  "S100A12",
  cor_results$rho[4],
  cor_results$p_value[4],
  cor_results$adjusted_p[4]
)


activity_all <- (
  p_MMP9_activity |
    p_S100A12_activity
)

activity_all


bdcaf_all <- (
  p_MMP9_bdcaf |
    p_S100A12_bdcaf
)

bdcaf_all







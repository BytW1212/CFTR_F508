# ==============================================================================
# Project: Transcriptomic Profiling of CFTR F508del Mutation
# Description: Differential expression analysis of GSE15568 using limma.
# ==============================================================================

#1. Load required packages
library(GEOquery)
library(limma)
library(hgu133a.db)
library(ggplot2)

#2. Download and extract data
gse <- getGEO("GSE15568", destdir = ".", getGPL = FALSE)
gse_data <- gse[[1]]

exp_matrix <- exprs(gse_data)
clinical_data <- pData(gse_data)

#3. Group samples
clinical_data$group <- ifelse(grepl("non-CF", clinical_data$title), "non-CF", "D508")
clinical_data$group <- factor(clinical_data$group, levels = c("non-CF", "D508"))

print("Sample group amount：")
print(table(clinical_data$group))

#4. Differential expression analysis (limma)
design <- model.matrix(~ 0 + clinical_data$group)

colnames(design) <- c("Control", "CF")

fit <- lmFit(exp_matrix, design)
contrast_matrix <- makeContrasts(CF_vs_Control = CF - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)

#5. Probe annotation
diff_results$Gene_Symbol <- mapIds(hgu133a.db,
                                   keys = rownames(diff_results),
                                   column = "SYMBOL",
                                   keytype = "PROBEID",
                                   multiVals = "first")

diff_results_clean <- diff_results[!is.na(diff_results$Gene_Symbol), ]

#6. Data visualization
diff_results_clean$Significance <- "Not Significant"
diff_results_clean$Significance[diff_results_clean$logFC > 1 & diff_results_clean$adj.P.Val < 0.05] <- "Up-regulated"
diff_results_clean$Significance[diff_results_clean$logFC < -1 & diff_results_clean$adj.P.Val < 0.05] <- "Down-regulated"

volcano_plot <- ggplot(diff_results_clean, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = "blue", 
                                "Not Significant" = "grey", 
                                "Up-regulated" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: CF vs Non-CF Control",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

print(volcano_plot)
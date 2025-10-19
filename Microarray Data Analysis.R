#########################################################
#  Differential Gene Expression Analysis - GSE13911
#  Author: Mahmoud Ahmed Abohassan
#  Date: October 2025
#########################################################

#--------------------------------------------------------
# 1. Setup Working Directory and Folders
#--------------------------------------------------------
setwd("G:/GSE13911")

# Create folders if they do not exist
dir.create("Results", showWarnings = FALSE)
dir.create("Result_Plots", showWarnings = FALSE)

# Confirm working directory
getwd()

#--------------------------------------------------------
# 2. Install and Load Required Packages
#--------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "GEOquery", "AnnotationDbi", "hgu133plus2.db"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

#--------------------------------------------------------
# 3. Load GSE13911 Dataset
#--------------------------------------------------------
# Download expression data from GEO
gse <- getGEO("GSE13911", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])
phenoData <- pData(gse[[1]])

#--------------------------------------------------------
# 4. Probe-to-Gene Mapping
#--------------------------------------------------------
probe_ids <- rownames(exprSet)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Combine expression and gene annotation
expr_df <- as.data.frame(exprSet) %>%
  rownames_to_column("PROBEID") %>%
  mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  filter(!is.na(SYMBOL))

# Collapse multiple probes per gene by averaging
expr_avg <- limma::avereps(
  as.matrix(expr_df[, -c(1, ncol(expr_df))]),
  ID = expr_df$SYMBOL
)

#--------------------------------------------------------
# 5. Prepare Design Matrix
#--------------------------------------------------------
# Identify sample groups
group_labels <- as.factor(ifelse(grepl("Normal", phenoData$title, ignore.case = TRUE),
                                 "Normal", "Cancer"))

design <- model.matrix(~0 + group_labels)
colnames(design) <- levels(group_labels)

#--------------------------------------------------------
# 6. Differential Expression Analysis using Limma
#--------------------------------------------------------
fit <- lmFit(expr_avg, design)
contrast_matrix <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, coef = "Cancer_vs_Normal", number = Inf, adjust.method = "BH")

#--------------------------------------------------------
# 7. Identify Upregulated and Downregulated Genes
#--------------------------------------------------------
deg_results$threshold <- ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "Not Significant")
)

up_genes <- subset(deg_results, threshold == "Upregulated")
down_genes <- subset(deg_results, threshold == "Downregulated")

# Save DEG results
write.csv(deg_results, "Results/All_DEGs.csv", row.names = TRUE)
write.csv(up_genes, "Results/Upregulated_DEGs.csv", row.names = TRUE)
write.csv(down_genes, "Results/Downregulated_DEGs.csv", row.names = TRUE)

#--------------------------------------------------------
# 8. Volcano Plot
#--------------------------------------------------------
png("Result_Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot - GSE13911",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value",
       color = "Regulation")

dev.off()

#--------------------------------------------------------
# 9. Heatmap of Top 25 DEGs
#--------------------------------------------------------
top25_genes <- rownames(deg_results[order(deg_results$adj.P.Val), ])[1:25]
heatmap_data <- expr_avg[top25_genes, ]

png("G:/GSE13911/Result_Plots/heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 25 Differentially Expressed Genes - GSE13911"
)
dev.off()
dim(heatmap_data)
head(heatmap_data)
#--------------------------------------------------------
# 10. Short Summary
#--------------------------------------------------------
cat("------------------------------------------------------------\n")
cat("Summary:\n")
cat("Multiple probes mapping to the same gene were averaged using limma::avereps().\n")
cat("Contrast performed: Cancer vs Normal gastric tissue.\n")
cat("Upregulated genes: ", nrow(up_genes), "\n")
cat("Downregulated genes: ", nrow(down_genes), "\n")
cat("All results saved in G:/GSE13911/Results and plots in G:/GSE13911/Result_Plots\n")
cat("------------------------------------------------------------\n")

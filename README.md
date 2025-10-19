🧬 Differential Gene Expression Analysis — GSE13911
📖 Description

This repository contains R code and results for analyzing the GSE13911 microarray dataset.
The goal is to identify differentially expressed genes (DEGs) between gastric cancer and normal gastric tissue samples using the Limma pipeline.

⚙️ Workflow

Downloaded dataset GSE13911 using GEOquery.

Annotated probe IDs to gene symbols using hgu133plus2.db.

Averaged multiple probes mapping to the same gene.

Performed differential expression analysis (cancer vs normal) using limma.

Created:

Volcano plot showing up/downregulated genes

Heatmap for the top 25 DEGs

📊 Results Summary

Accession ID: GSE13911

Annotation package: hgu133plus2.db

Comparison: cancer_vs_normal

Upregulated genes: ~640

Downregulated genes: ~520

Duplicate probes were handled by averaging their expression values per gene.

📂 Outputs
File	Description
Results/DEG_results_complete.csv	All DEGs with logFC and p-values
Results/Upregulated_DEGs.csv	Upregulated genes
Results/Downregulated_DEGs.csv	Downregulated genes
Result_Plots/volcano_plot_GSE13911.png	Volcano plot
Result_Plots/heatmap_top25_DEGs.png	Heatmap of top 25 DEGs
🧰 Required Packages
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tibble)

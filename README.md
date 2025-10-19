Differential Gene Expression Analysis — GSE13911
Project Overview

This project performs a Differential Gene Expression (DGE) analysis on the publicly available microarray dataset GSE13911, which investigates gene expression differences between gastric cancer and normal gastric tissue samples.
The workflow follows a standard bioinformatics pipeline using the Limma package in R.

⚙️ Workflow Summary

Dataset Import:
The dataset GSE13911 was downloaded using the GEOquery package.

Annotation:
Probe IDs were mapped to gene symbols using the annotation package
hgu133plus2.db.

Handling Duplicate Probes:
Multiple probes mapped to the same gene were averaged to obtain a single representative expression value for each gene.

Differential Expression Analysis:
The limma package was used to identify differentially expressed genes (DEGs) between cancer and normal samples.

Visualization:

A Volcano Plot was generated to display significantly upregulated and downregulated genes.

A Heatmap of the top 25 DEGs was created to visualize gene expression clustering.

Result Export:
DEG results were saved as CSV files, and plots were exported as PNG images in the Result_Plots folder.

📊 Results Summary

Dataset Accession ID: GSE13911

Annotation Package Used: hgu133plus2.db

Comparison Performed: cancer_vs_normal

Metric	Value
Upregulated Genes	~640
Downregulated Genes	~520
Adjusted p-value cutoff	< 0.05
Log2 Fold Change cutoff	> 1

Explanation:

Several probes mapped to the same gene due to targeting different gene regions or isoforms.

Duplicate probes were resolved by averaging their expression values.

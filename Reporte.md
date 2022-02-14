### Reporte 

# Differential expression in Esophageal Adenocarcinoma on different stages 

### Introduction

Recount3 tool in R was used to select a determined set of sequences belonging to a differential expression experiment (SRP119465). This sequences were later on treated to be used in different statistical graphics that could lead to a better understanding of differential expression in conditions like disease stage of Esophageal Adenocarcinoma.

### Results

Differential expression between the different stages of the disease (IIA nad IIIB) and the assigned gene proportion, where stage IIA can be spotted to have a different expression [Disease_Group_vs_AGP.png] (https://github.com/Rodel-OL/rnaseq_2022_proyecto/blob/master/Disease_Group_vs_AGP.png)

A simple plot to compare the expression from different stages using the assigned gene prop field, that can be better seen in the first graph [Plot_assigned_gene_prop_vs_d_stage] (https://github.com/Rodel-OL/rnaseq_2022_proyecto/blob/master/Plot_assigned_gene_prop_vs_d_stage.png)

log fold change and average log-expression for the "disease" attribute [Average_log_expression_disease.png] (https://github.com/Rodel-OL/rnaseq_2022_proyecto/blob/master/Average_log_expression_disease.png)

Volcano plot with 3 remarked genes, from which CLDN3 seems downregulated, while both BEX3 and GPNMB are upregulated [Volcano_val_3.png] (https://github.com/Rodel-OL/rnaseq_2022_proyecto/blob/master/Volcano_val_3.png)

Volcano plot with 5 remarked genes, where 2 are downregulated (CLDN4, LGALS4) and 3 are upregulated (BEX3, GPNMB, KRT5) [Volcano_5_val.png] (https://github.com/Rodel-OL/rnaseq_2022_proyecto/blob/master/Volcano_5_val.png)

### Libraries used

recount3

limma

edgeR

ggplot2

Code lines used to get data from recount3 

> folder > R
>
> file > 02_recount3.R

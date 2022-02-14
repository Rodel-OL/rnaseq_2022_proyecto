## load recount3
library("recount3")

human_projects <- available_projects()

## Data from SRP119465 was used

proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## Before sending, search for and select "SRP119465"
## Verifying only one line was sent
stopifnot(nrow(proj_info_interactive) == 1)
## Creating RSE object
rse_gene_interactive <- create_rse(proj_info_interactive)
## Adding assay counts
assay(rse_gene_interactive, "counts") <- compute_read_counts(rse_gene_interactive)
## Using SRA attributes
rse_gene_SRP119465 <- expand_sra_attributes(rse_gene_interactive)
## What attributes are added from sra
colData(rse_gene_SRP119465)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP119465)))
]

## Change characters to factor and use of as.numeric for later use
rse_gene_SRP119465$sra_attribute.age <- as.numeric(rse_gene_SRP119465$sra_attribute.age)
rse_gene_SRP119465$sra_attribute.disease <- factor(rse_gene_SRP119465$sra_attribute.disease)
rse_gene_SRP119465$sra_attribute.disease_stage <- factor(rse_gene_SRP119465$sra_attribute.disease_stage)
rse_gene_SRP119465$sra_attribute.sex <- factor(rse_gene_SRP119465$sra_attribute.sex)
rse_gene_SRP119465$sra_attribute.tissue <- factor(rse_gene_SRP119465$sra_attribute.tissue)

## Summary
summary(as.data.frame(colData(rse_gene_SRP119465)[
  ,
  grepl("^sra_attribute.[age|disease|disease_stage|sex|tissue]", colnames(colData(rse_gene_SRP119465)))
]))

## Gene information Quality, use of asigned_gene_prop field
rse_gene_SRP119465$assigned_gene_prop <- rse_gene_SRP119465$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP119465$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP119465$assigned_gene_prop)

## Plotting disease stage and assigned gene prop
with(colData(rse_gene_SRP119465), plot(assigned_gene_prop, sra_attribute.disease_stage))

## Histogram for bad quality gene information from samples
hist(rse_gene_SRP119465$assigned_gene_prop)

gene_means <- rowMeans(assay(rse_gene_SRP119465, "counts"))
summary(gene_means)
## Since we got 0.01 at 1st Qu. we may eliminate genes with less
## Test without changes
rse_gene_SRP119465_unfiltered <- rse_gene_SRP119465
## Elimination
rse_gene_SRP119465 <- rse_gene_SRP119465[gene_means > 0.01, ]
## Percentage of retained genes [72.61%]
round(nrow(rse_gene_SRP119465) / nrow(rse_gene_SRP119465_unfiltered) * 100, 2)

## Data normalization
## Use of edgeR
library("edgeR")
dge <- DGEList(
  counts = assay(rse_gene_SRP119465, "counts"),
  genes = rowData(rse_gene_SRP119465)
)
dge <- calcNormFactors(dge)

## Differential Expression

library("ggplot2")

## Disease stage Group plotting
ggplot(as.data.frame(colData(rse_gene_SRP119465)), aes(y = assigned_gene_prop, x = sra_attribute.disease_stage)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Disease Group")

# library(ExploreModelMatrix)
#
# vd <- VisualizeDesign(
#   sampleData = as.data.frame(colData(rse_gene_SRP119465)),
#   designFormula = ~ assigned_gene_prop + assigned_gene_prop:sra_attribute.disease_stage + assigned_gene_prop:sra_attribute.disease + assigned_gene_prop:sra_attribute.age,
#   textSizeFitted = 4
# )
# cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

## Need to add a statistical model first, to then use limma
mod <- model.matrix(~ sra_attribute.disease + sra_attribute.disease_stage + assigned_gene_prop,
                    data = colData(rse_gene_SRP119465)
)
colnames(mod)

## Use of limma and plotting
library("limma")
vGene <- voom(dge, mod, plot = TRUE)
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP119465),
  sort.by = "none"
)

## Use of volcano plot
volcanoplot(eb_results, coef = 2, highlight = 5, names = de_results$gene_name)
## Downregulated[CLDN3, LGALS4] Upregulated[BEX3, KRT5, GPNMB]
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
## Downregulated[CLDN3] Upregulated[BEX3, GPNMB]

library("iSEE")
iSEE::iSEE(rse_gene_SRP119465)

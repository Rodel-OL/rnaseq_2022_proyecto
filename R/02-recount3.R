## load recount3
library("recount3")

human_projects <- available_projects()

## SRP119465

proj_info_interactive <- interactiveDisplayBase::display(human_projects)

## Verifying only one line was sent
stopifnot(nrow(proj_info_interactive) == 1)
## Creating RSE object
rse_gene_interactive <- create_rse(proj_info_interactive)
## Adding assay counts
assay(rse_gene_interactive, "counts") <- compute_read_counts(rse_gene_interactive)
## Using SRA attributes
rse_gene_SRP119465 <- expand_sra_attributes(rse_gene_interactive)
## What attributes are added
colData(rse_gene_SRP119465)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP119465)))
]

## Change characters to something "plottable"
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

## Plotting
with(colData(rse_gene_SRP119465), plot(assigned_gene_prop, sra_attribute.disease))
#with(colData(rse_gene_SRP119465), plot(assigned_gene_prop, recount_qc.intron_sum_%))

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
## Use of shiny?
## Need to add a statistical model first, to then use limma

library("iSEE")
iSEE::iSEE(rse_gene_SRP119465)

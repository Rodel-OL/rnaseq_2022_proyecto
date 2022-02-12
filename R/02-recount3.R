## Cargar recount3
library("recount3")

human_projects <- available_projects()

SRP133965

proj_info_interactive <- interactiveDisplayBase::display(human_projects)

## Verifying only one line was sent
stopifnot(nrow(proj_info_interactive) == 1)
## Creating RSE object
rse_gene_interactive <- create_rse(proj_info_interactive)

assay(rse_gene_interactive, "counts") <- compute_read_counts(rse_gene_interactive)

rse_gene_SRP133965 <- expand_sra_attributes(rse_gene_interactive)
colData(rse_gene_interactive)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_interactive)))
]

library("iSEE")
iSEE::iSEE(rse_gene_SRP133965)

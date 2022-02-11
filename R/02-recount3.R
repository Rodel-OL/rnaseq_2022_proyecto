## Cargar recount3
library("recount3")

human_projects <- available_projects()

SRP169065

proj_info <- subset(
  human_projects,
  project == "SRP169065" & project_type == "data_sources"
)

rse_gene_SRP169065 <- create_rse(proj_info)

assay(rse_gene_SRP169065, "counts") <- compute_read_counts(rse_gene_SRP169065)

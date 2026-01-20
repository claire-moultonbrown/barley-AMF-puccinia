##############################################################################

# 02_go_profiler

##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(gprofiler2)
})

proj <- "/path/to/project"
deg <- read_csv(file.path(proj, "results_edgeR_deg_amf.csv"))

sig_genes <- deg %>%
  filter(is_deg) %>%
  pull(gene_id)

# Adjust organism if you are using Ensembl Plants or custom identifiers.
# If gene IDs are MorexV3 gene IDs, you may need to map them first.
gp <- gost(
  query = sig_genes,
  organism = "hordeum_vulgare",
  correction_method = "fdr"
)

saveRDS(gp, file.path(proj, "gprofiler_amf_gost.rds"))
write_csv(gp$result, file.path(proj, "gprofiler_amf_results.csv"))

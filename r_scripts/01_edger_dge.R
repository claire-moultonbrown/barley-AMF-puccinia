##############################################################################

# 01_edger_dge

##############################################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(tidyverse)
})

proj <- "/path/to/project"
count_dir <- file.path(proj, "counts", "all_counts")
meta <- read_csv("data_raw/metadata.csv"))

# Read HTSeq count tables, keep gene counts, drop HTSeq summary rows
read_htseq <- function(f) {
  x <- read_tsv(f, col_names = c("gene_id", "count"), show_col_types = FALSE)
  x %>%
    filter(!str_starts(gene_id, "__")) %>%
    mutate(count = as.integer(count))
}

files <- file.path(count_dir, paste0(meta$sample_id, "_counts.txt"))
names(files) <- meta$sample_id

counts_list <- lapply(files, read_htseq)
gene_ids <- counts_list[[1]]$gene_id

count_mat <- do.call(cbind, lapply(counts_list, `[[`, "count"))
rownames(count_mat) <- gene_ids
colnames(count_mat) <- names(files)

# edgeR object
group <- factor(meta$Treatment, levels = c("C", "A", "P", "AP"))
y <- DGEList(counts = count_mat, group = group)

# Filter CPM < 1 (keep genes expressed in enough samples)
keep <- rowSums(cpm(y) >= 1) >= 1
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalisation
y <- calcNormFactors(y, method = "TMM")

# Design matrix (2x2 factorial to match Methods)
design <- model.matrix(~ AMF * Puccinia, data = meta)

# Estimate dispersion and fit GLM
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Example contrasts: main effects and interaction
res_amf <- glmQLFTest(fit, coef = "AMF")
res_puc <- glmQLFTest(fit, coef = "Puccinia")
res_int <- glmQLFTest(fit, coef = "AMF:Puccinia")

# Add logCPM and apply thresholds: FDR < 0.05, FC > 2, logCPM > 2
as_deg <- function(res) {
  tt <- topTags(res, n = Inf)$table %>%
    rownames_to_column("gene_id") %>%
    as_tibble()
  tt %>%
    mutate(
      FDR = p.adjust(PValue, method = "BH"),
      FC = 2^abs(logFC),
      is_deg = (FDR < 0.05) & (FC > 2) & (logCPM > 2)
    )
}

deg_amf <- as_deg(res_amf)
deg_puc <- as_deg(res_puc)
deg_int <- as_deg(res_int)

write_csv(deg_amf, file.path(proj, "results_edgeR_deg_amf.csv"))
write_csv(deg_puc, file.path(proj, "results_edgeR_deg_puccinia.csv"))
write_csv(deg_int, file.path(proj, "results_edgeR_deg_interaction.csv"))

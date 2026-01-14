# ##############################################################################

# Figure 4 – transcriptomics

# ##############################################################################

# ==============================
# Setup
# ==============================
rm(list = ls())

dir_raw       <- "data_raw"
dir_processed <- "data_processed"
dir_figures   <- "figures"
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures,   showWarnings = FALSE, recursive = TRUE)

# ==============================
# Load libraries
# ==============================
# Data handling and plotting (ggplot-based figures)
library(tidyverse)

# Multi-panel figure assembly
library(cowplot)

# Heatmap (ComplexHeatmap-based)
library(ComplexHeatmap)
library(circlize)   # for colorRamp2()
library(grid)       # for gpar(), unit()

# ==============================
# Read and format data
# ==============================
volcano     <- read.csv(file.path(dir_raw, "fig4a_transcriptomics_volcano.csv"))
amf_induced <- read.csv(file.path(dir_raw, "fig4b_transcriptomics_amf_induced_go_terms.csv"))
degs        <- read.csv(file.path(dir_raw, "fig4c_transcriptomics_expression_by_treatment.csv"))
anno        <- read.csv(file.path(dir_raw, "fig4c_transcriptomics_cluster_annotation.csv"))
clusters    <- read.csv(file.path(dir_raw, "fig4d_transcriptomics_enriched_cluster_go_terms.csv"))

cols <- c("Regulation", "term_name")
amf_induced[cols] <- lapply(amf_induced[cols], factor)

clusters$Cluster <- as.factor(clusters$Cluster)

# ==============================
# Plot parameters
# ==============================

theme_claire <- theme(
  panel.background = element_rect(fill = NA),
  legend.title     = element_blank(),
  axis.text        = element_text(size = 7,  face = "bold", colour = "black"),
  axis.title       = element_text(size = 10, face = "bold", colour = "black"),
  legend.text      = element_text(size = 14, face = "bold"),
  panel.border     = element_blank(),
  axis.line        = element_line(colour = "black", linewidth = 0.5),
  legend.key       = element_blank()
)

theme_bubble <- theme(
  panel.background = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  axis.text.x      = element_blank(),
  axis.text        = element_text(size = 7, face = "bold", colour = "black"),
  axis.title.y     = element_blank(),
  strip.placement  = "outside",
  strip.text       = element_text(angle = 0, size = 8, face = "bold", colour = "black"),
  strip.background = element_rect(fill = "grey", colour = "black", linewidth = 0.5),
  legend.title     = element_text(size = 8, face = "bold", colour = "black"),
  legend.text      = element_text(size = 7, face = "bold"),
  legend.key       = element_blank()
)

bubble_gradient_cols <- c("blue", "yellow", "red")

make_bubble <- function(df, y_col, size_col, colour_col, facet_formula, colour_name, size_name, theme_obj) {
  ggplot(df, aes(x = "", y = reorder(.data[[y_col]], .data[[size_col]], FUN = max))) +
    geom_point(aes(size = .data[[size_col]], colour = .data[[colour_col]]), alpha = 0.9) +
    facet_grid(facet_formula, scales = "free_y", space = "free_y") +
    scale_colour_gradientn(colours = bubble_gradient_cols, name = colour_name) +
    scale_size_continuous(name = size_name) +
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw() + theme_obj
}

# ==============================
# Figure 4A volcano plot
# ==============================

# thresholds
fc_thr <- 2
p_thr  <- 0.05

volcano2 <- volcano %>%
  mutate(
    log2p = -log2(pvalue),
    status = case_when(
      pvalue < p_thr & logFC >=  fc_thr ~ "Up",
      pvalue < p_thr & logFC <= -fc_thr ~ "Down",
      TRUE                            ~ "NS"
    ))

fig4A <- ggplot(volcano2, aes(x = logFC, y = log2p)) +
  geom_point(aes(colour = status), alpha = 0.5, size = 1) +
  scale_colour_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  labs(
    x = bquote(log[2]*"FC"),
    y = bquote(log[2]*"p-value"),
    colour = "Regulation"
  ) +
  theme_claire +
  guides(colour = "none")
fig4A

# ==============================
# Figure 4B bubble plot of AMF induced GO terms
# ==============================

amf_induced2 <- amf_induced %>%
  mutate(
    neglog10padj = -log10(padj),
    Regulation = factor(Regulation,
      levels = c("Upregulated with AMF", "Downregulated with AMF"))
  )

fig4B <- make_bubble(
  amf_induced2,
  y_col         = "term_name",
  size_col      = "gene_number",
  colour_col    = "neglog10padj",
  facet_formula = Regulation ~ .,
  colour_name   = "-log10(adj p-value)",
  size_name     = "Count",
  theme_obj     = theme_bubble + 
    theme(strip.text = element_text(angle = 270))
)
fig4B

# ==============================
# Figure 4C heatmap of DEGs
# ==============================

# ------------------------------
# tidy
# ------------------------------
mat <- degs %>%
  dplyr::select(Gene.ID, C, A, P, AP) %>%
  mutate(across(c(C, A, P, AP), ~ suppressWarnings(as.numeric(.x)))) %>%
  distinct(Gene.ID, .keep_all = TRUE) %>%
  column_to_rownames("Gene.ID") %>%
  as.matrix()

# Force column order C A P AP
mat <- mat[, c("C", "A", "P", "AP"), drop = FALSE]

row_anno <- anno %>%
  dplyr::select(Gene.ID, Cluster) %>%
  distinct(Gene.ID, .keep_all = TRUE) %>%
  mutate(
    Cluster = ifelse(is.na(Cluster) | Cluster == "", "Unknown", Cluster),
    Cluster = factor(Cluster)
  ) %>%
  column_to_rownames("Gene.ID")

# Match and fill any missing rows
row_anno <- row_anno[rownames(mat), , drop = FALSE]
row_anno$Cluster <- addNA(row_anno$Cluster); levels(row_anno$Cluster)[is.na(levels(row_anno$Cluster))] <- "Unknown"

# ------------------------------
# Transform + row scaling
# ------------------------------
mat_log <- log2(mat + 1)
mat_scaled <- t(scale(t(mat_log)))
mat_scaled[is.nan(mat_scaled)] <- NA

# Optional: drop rows that became all NA after scaling
keep <- rowSums(is.finite(mat_scaled)) > 0
mat_scaled <- mat_scaled[keep, , drop = FALSE]

# ------------------------------
# Row clustering + cut into 7 clusters
# ------------------------------
k_clusters <- 7
row_hc <- hclust(dist(mat_scaled), method = "complete")
cl <- cutree(row_hc, k = k_clusters)

cluster_df <- tibble(
  Gene.ID = names(cl),
  Cluster = paste0("c", cl)
) %>%
  arrange(as.integer(sub("^c", "", Cluster)), Gene.ID)

write_csv(cluster_df, 'data_processed/degs_7clusters_assignments.csv')

# ------------------------------
# Heatmap colours
# ------------------------------
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#D73027")
)

cluster_cols <- c(
  "Cluster 1"      = "#1B9E77",
  "Cluster 2"      = "#D95F02",
  "Cluster 3"      = "#7570B3",
  "Cluster 4"      = "#E7298A",
  "Cluster 5"      = "#66A61E",
  "Cluster 6"      = "#E6AB02",
  "Cluster 7"      = "#A6CEE3"
)

# ------------------------------
# Row annotation for cluster
# ------------------------------

ha <- rowAnnotation(
  df = row_anno,
  col = list(Cluster = cluster_cols),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 14),
  annotation_legend_param = list(
    Cluster = list(
      title_gp  = gpar(fontsize = 18, fontface = "bold"),
      labels_gp = gpar(fontsize = 16)
    )
  )
)

# ------------------------------
# Heatmap with dendrogram + visible gaps between 7 clusters
# ------------------------------
ht <- Heatmap(
  mat_scaled,
  name = "Row z-score\nlog2(expr + 1)",
  col = col_fun,
  cluster_rows = row_hc,
  cluster_columns = FALSE,
  row_split = k_clusters,
  row_gap = unit(2, "mm"),
  row_title = NULL,
  left_annotation = ha,
  show_row_names = FALSE,
  rect_gp = gpar(col = NA),
  column_names_gp = gpar(fontsize = 18, fontface = "bold"),
  column_names_rot = 0,
  heatmap_legend_param = list(
    title_gp   = gpar(fontsize = 18, fontface = "bold"),
    labels_gp  = gpar(fontsize = 16),
    legend_gap = unit(8, "mm")
  )
)

# Render heatmap to a temporary PNG
heatmap_png <- file.path(dir_figures, "fig4c_heatmap.png")

png(heatmap_png, width = 2200, height = 2600, res = 300)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Bring heatmap back in as an image panel for cowplot
fig4C <- ggdraw() + draw_image(heatmap_png, scale = 1)
fig4C


# ==============================
# Figure 4D bubble plot of clustered DEG GO terms
# ==============================

clusters2 <- clusters %>%
  mutate(neglog10padj = -log10(pvalue))

fig4D <- make_bubble(
  clusters2,
  y_col         = "GO.term",
  size_col      = "number",
  colour_col    = "neglog10padj",
  facet_formula = ~ Cluster,
  colour_name   = "-log10(adj p-value)",
  size_name     = "Count",
  theme_obj     = theme_bubble)
fig4D

# ==============================
# Assemble multi-panel figure
# ==============================
fig4_multiplot <- ggdraw() +
  draw_plot(fig4A, x = 0.05, y = 0.70, width = 0.27, height = 0.30) +
  draw_plot(fig4B, x = 0.50, y = 0.50, width = 0.48, height = 0.48) +
  draw_plot(fig4C, x = 0.00, y = 0.00, width = 0.40, height = 0.70) +
  draw_plot(fig4D, x = 0.40, y = 0.00, width = 0.60, height = 0.48) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)", "(D)"),
    size  = 16,
    x     = c(0, 0.4, 0, 0.4),
    y     = c(1, 1, 0.77, 0.5 )
  )
fig4_multiplot

# ==============================
# Save
# ==============================
ggsave(
  filename = file.path(dir_figures, "Figure_4.png"),
  plot     = fig4_multiplot,
  width    = 26,
  height   = 18,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

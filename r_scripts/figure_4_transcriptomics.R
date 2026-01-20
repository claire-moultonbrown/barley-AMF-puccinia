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
library(tidyverse)
library(cowplot)

# Heatmap (pheatmap-based, as per your template)
library(RColorBrewer)
library(pheatmap)
library(dendextend)

# ==============================
# Read and format data
# ==============================
volcano     <- read.csv(file.path(dir_raw, "fig4a_transcriptomics_volcano.csv"))
amf_induced <- read.csv(file.path(dir_raw, "fig4b_transcriptomics_amf_induced_go_terms.csv"))
degs        <- read.csv(file.path(dir_raw, "fig4c_transcriptomics_expression_by_treatment.csv"))
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

get_pheatmap_tree_row <- function(ph_obj) {
  # pheatmap return type varies by version / context; handle list + S4
  if (!is.null(tryCatch(ph_obj$tree_row, error = function(e) NULL))) return(ph_obj$tree_row)
  if (!is.null(tryCatch(ph_obj[["tree_row"]], error = function(e) NULL))) return(ph_obj[["tree_row"]])
  if (methods::is(ph_obj, "pheatmap") && !is.null(tryCatch(ph_obj@tree_row, error = function(e) NULL))) return(ph_obj@tree_row)
  stop("Could not extract tree_row from pheatmap object (unexpected pheatmap return type).")
}


# ==============================
# Figure 4A volcano plot
# ==============================
#log2 fold change >1 so fc > 2
fc_thr <- 1
p_thr  <- 0.05

volcano2 <- volcano %>%
  mutate(
    log2p = -log2(pvalue),
    status = case_when(
      pvalue < p_thr & logFC >=  fc_thr ~ "Up",
      pvalue < p_thr & logFC <= -fc_thr ~ "Down",
      TRUE                             ~ "NS"
    )
  )

fig4A <- ggplot(volcano2, aes(x = logFC, y = log2p)) +
  geom_point(aes(colour = status), alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
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
    Regulation = factor(Regulation, levels = c("Upregulated with AMF", "Downregulated with AMF"))
  )

fig4B <- make_bubble(
  amf_induced2,
  y_col         = "term_name",
  size_col      = "gene_number",
  colour_col    = "neglog10padj",
  facet_formula = Regulation ~ .,
  colour_name   = "-log10(adj p-value)",
  size_name     = "Count",
  theme_obj     = theme_bubble + theme(strip.text = element_text(angle = 270))
)
fig4B

# ==============================
# Figure 4C heatmap of DEGs (pheatmap method)
# ==============================

# ------------------------------
# Build matrix (C/A/P/AP)
# ------------------------------
mat <- degs %>%
  dplyr::select(Gene.ID, C, A, P, AP) %>%
  mutate(across(c(C, A, P, AP), ~ suppressWarnings(as.numeric(.x)))) %>%
  distinct(Gene.ID, .keep_all = TRUE) %>%
  column_to_rownames("Gene.ID") %>%
  as.matrix()

mat <- mat[, c("C", "A", "P", "AP"), drop = FALSE]

# ------------------------------
# pheatmap settings to match your template
# ------------------------------
labels_col <- c("C", "A", "P", "AP")

hm_cols   <- rev(brewer.pal(n = 9, name = "RdBu"))
hm_breaks <- c(-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2)

# Cluster colours (7 clusters, like your earlier approach)
ann_colors <- list(
  Cluster = c(
    "C1" = "#CE3D32FF",
    "C2" = "#FF14FF",
    "C3" = "#749B58FF",
    "C4" = "#466983FF",
    "C5" = "#BA6338FF",
    "C6" = "#5DB1DDFF",
    "C7" = "#8022FF"
  )
)

# ------------------------------
# Draw heatmap (saved to file)
# ------------------------------
heatmap_png <- file.path(dir_figures, "fig4c_heatmap.png")

# Note:
# - scale = "row" gives row z-scores (as in your template)
# - clustering_method = "ward.D2" matches your template and your previous ComplexHeatmap choice
# - cutree_rows = 7 splits dendrogram into 7 clusters (exported below)
# - annotation_row uses clusters derived from the pheatmap tree (so no annotation file needed)

# We generate the heatmap first WITHOUT annotation_row, extract clusters, then redraw with annotation.
tmp_hm <- pheatmap::pheatmap(
  mat,
  color                    = hm_cols,
  breaks                   = hm_breaks,
  scale                    = "row",
  cluster_rows             = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method        = "ward.D2",
  cluster_cols             = FALSE,
  show_rownames            = FALSE,
  show_colnames            = TRUE,
  labels_col               = labels_col,
  fontsize                 = 12,
  fontsize_col             = 12,
  cellwidth                = 50,
  border_color             = NA,
  legend                   = TRUE,
  annotation_names_col     = FALSE
)

# ------------------------------
# Extract clusters from dendrogram + save
# ------------------------------
k_clusters <- 7
row_tree <- get_pheatmap_tree_row(tmp_hm)
tree_k <- sort(cutree(row_tree, k = k_clusters))

cluster_df <- tibble(
  Gene.ID = names(tree_k),
  Cluster = paste0("C", as.integer(tree_k))
) %>%
  arrange(as.integer(sub("^C", "", Cluster)), Gene.ID)

write_csv(cluster_df, file.path(dir_processed, "fig4c_deg_clusters_cutree_k7.csv"))

# Build row annotation from the clustering (no external annotation file)
row_anno <- cluster_df %>%
  column_to_rownames("Gene.ID")

# Ensure row order matches matrix rows
row_anno <- row_anno[rownames(mat), , drop = FALSE]

# ------------------------------
# Redraw heatmap with row annotation and save to PNG
# ------------------------------
# Render heatmap to a PNG (pheatmap draws directly to the active device)
heatmap_png <- file.path(dir_figures, "fig4c_heatmap.png")

png(heatmap_png, width = 2200, height = 2600, res = 300)

pheatmap(
  mat,
  color                    = hm_cols,
  breaks                   = hm_breaks,
  scale                    = "row",
  cluster_rows             = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method        = "ward.D2",
  cluster_cols             = FALSE,
  show_rownames            = FALSE,
  show_colnames            = TRUE,
  labels_col               = labels_col,
  fontsize                 = 12,
  fontsize_col             = 12,
  cellwidth                = 50,
  border_color             = NA,
  legend                   = TRUE,
  annotation_row           = row_anno,
  annotation_names_row     = FALSE,
  annotation_colors        = ann_colors
)

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
  y_col         = "GO_term",
  size_col      = "number",
  colour_col    = "neglog10padj",
  facet_formula = ~ Cluster,
  colour_name   = "-log10(adj p-value)",
  size_name     = "Count",
  theme_obj     = theme_bubble
)
fig4D

# ==============================
# Assemble multi-panel figure
# ==============================
fig4_multiplot <- ggdraw() +
  draw_plot(fig4A, x = 0.05, y = 0.70, width = 0.27, height = 0.30) +
  draw_plot(fig4B, x = 0.50, y = 0.50, width = 0.48, height = 0.48) +
  draw_plot(fig4C, x = 0.00, y = 0.00, width = 0.44, height = 0.70) +
  draw_plot(fig4D, x = 0.40, y = 0.00, width = 0.60, height = 0.48) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)", "(D)"),
    size  = 16,
    x     = c(0, 0.4, 0, 0.4),
    y     = c(1, 1, 0.7, 0.5)
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

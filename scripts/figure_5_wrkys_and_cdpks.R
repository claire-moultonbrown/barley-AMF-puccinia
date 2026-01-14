##############################################################################

# Figure 5 – WRKYs and CDPKs

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
# Libraries
# ==============================
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ==============================
# Read data
# ==============================
wrky    <- readr::read_csv(file.path(dir_raw, "fig5a_transcriptomics_wrky_expression.csv"), show_col_types = FALSE)
wrky_go <- readr::read_csv(file.path(dir_raw, "fig5b_transcriptomics_wrky_go_terms.csv"),   show_col_types = FALSE)
cdpk    <- readr::read_csv(file.path(dir_raw, "fig5c_transcriptomics_cdpk_expression.csv"), show_col_types = FALSE)

# ==============================
# Global settings for plotting
# ==============================
expr_cols <- c("C", "A", "P", "AP")

# heatmap colour function (row z-scores)
hm_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027"))

# bubble palette
bubble_gradient_cols <- c("blue", "yellow", "red")

# legend text sizes (match Fig 5A)
legend_title_gp <- gpar(fontsize = 18, fontface = "bold")
legend_labels_gp <- gpar(fontsize = 16, fontface = "plain")

# shared bubble theme
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
  legend.key       = element_blank(),
  legend.position  = "right",
  legend.box       = "horizontal",
  legend.direction = "vertical"
)

# ==============================
# Helper functions
# ==============================
# make the expression matrix for heatmaps
prep_expr_matrix <- function(df, id_col, expr_cols) {
  df %>%
    dplyr::select(all_of(c(id_col, expr_cols))) %>%
    mutate(across(all_of(expr_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    column_to_rownames(id_col) %>%
    as.matrix() %>%
    (\(m) m[, expr_cols, drop = FALSE])()
}

# calculate expression scaled by row
row_zscore_log2p1 <- function(mat) {
  m <- log2(mat + 1)
  z <- t(scale(t(m)))
  z[is.nan(z)] <- NA
  keep <- rowSums(is.finite(z)) > 0
  z[keep, , drop = FALSE]
}

# render the heatmaps
render_heatmap_png <- function(ht, filename, width = 2200, height = 2600, res = 300) {
  png(filename, width = width, height = height, res = res)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  ggdraw() + draw_image(filename, scale = 1)
}

# make bubble plots
make_bubble <- function(df, y_col, size_col, colour_col, colour_name, size_name, theme_obj) {
  ggplot(df, aes(x = "", y = reorder(.data[[y_col]], .data[[size_col]], FUN = max))) +
    geom_point(aes(size = .data[[size_col]], colour = .data[[colour_col]]), alpha = 0.9) +
    scale_colour_gradientn(colours = bubble_gradient_cols, name = colour_name) +
    scale_size_continuous(name = size_name) +
    guides(size = guide_legend(order = 1), colour = guide_colourbar(order = 2)) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw() + theme_obj
}

# ==============================
# Figure 5A – WRKY heatmap (Role + Class)
# ==============================
mat_wrky <- prep_expr_matrix(wrky, id_col = "Name", expr_cols = expr_cols)
mat_wrky_scaled <- row_zscore_log2p1(mat_wrky)
row_hc_wrky <- hclust(dist(mat_wrky_scaled), method = "complete")

row_anno_wrky <- wrky %>%
  dplyr::select(Name, Role, Class) %>%
  distinct(Name, .keep_all = TRUE) %>%
  mutate(
    Role  = ifelse(is.na(Role)  | Role  == "" | Role  == "NA", "Unknown", Role),
    Class = ifelse(is.na(Class) | Class == "" | Class == "NA", "Unknown", Class),
    Role  = factor(Role),
    Class = factor(Class, levels = c("I", "II", "III"))
  ) %>%
  column_to_rownames("Name")

row_anno_wrky <- row_anno_wrky[rownames(mat_wrky_scaled), , drop = FALSE]

# palettes
role_cols <- c("Abiotic" = "yellow", 
               "Biotic"  = "#42BE65", 
               "N/A"     = "grey60", 
               "Unknown" = "grey80")
role_cols <- role_cols[levels(row_anno_wrky$Role)]

class_cols <- c("I"       = "#C6DBEF", 
                "II"      = "#6BAED6", 
                "III"     = "#2171B5", 
                "Unknown" = "grey80")

# if Unknown never occurs, this harmlessly drops it
class_cols <- class_cols[levels(droplevels(row_anno_wrky$Class))]

# Role closest to heatmap: df columns ordered Class then Role
row_anno_wrky <- row_anno_wrky[, c("Class", "Role"), drop = FALSE]

ha_wrky <- rowAnnotation(
  df  = row_anno_wrky,
  col = list(Class = class_cols, Role = role_cols),
  gp  = gpar(col = "black", lwd = 0.6),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 14),
  annotation_legend_param = list(
    Class = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp),
    Role  = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp)
  )
)

ht_wrky <- Heatmap(
  mat_wrky_scaled,
  name = "Row z-score\nlog2(expr + 1)",
  col  = hm_col_fun,
  cluster_rows = row_hc_wrky,
  cluster_columns = FALSE,
  column_order = expr_cols,
  left_annotation = ha_wrky,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_rot = 0,
  rect_gp = gpar(col = "black"),
  heatmap_legend_param = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp)
)

fig5A <- render_heatmap_png(
  ht_wrky,
  filename = file.path(dir_figures, "fig5a_heatmap.png")
)

# ==============================
# Figure 5B – WRKY GO bubble plot
# ==============================
wrky_go2 <- wrky_go %>%
  mutate(neglog10padj = -log10(p_value))

fig5B <- make_bubble(
  wrky_go2,
  y_col = "GO_name",
  size_col = "number",
  colour_col = "neglog10padj",
  colour_name = "-log10(adj p-value)",
  size_name = "Count",
  theme_obj = theme_bubble
)

# ==============================
# Figure 5C – CDPK heatmap (single annotation, no title)
# ==============================
mat_cdpk <- prep_expr_matrix(cdpk, id_col = "Name", expr_cols = expr_cols)
mat_cdpk_scaled <- row_zscore_log2p1(mat_cdpk)
row_hc_cdpk <- hclust(dist(mat_cdpk_scaled), method = "complete")

row_anno_cdpk <- cdpk %>%
  dplyr::select(Name, Annotation) %>%
  distinct(Name, .keep_all = TRUE) %>%
  mutate(
    Annotation = ifelse(is.na(Annotation) | Annotation == "" | Annotation == "NA", "Unknown", Annotation),
    Annotation = factor(
      Annotation,
      levels = c("Calmodulin-like protein", "Calmodulin", "Calcium-dependent protein kinase",
                 "CDPK-related protein", "Unknown"),
      labels = c("Calmodulin-like\nprotein", "Calmodulin", "Calcium-dependent\nprotein kinase",
                 "CDPK-related\nprotein", "Unknown")
    )
  ) %>%
  column_to_rownames("Name")

row_anno_cdpk <- row_anno_cdpk[rownames(mat_cdpk_scaled), , drop = FALSE]

anno_cols <- c(
  "Calmodulin-like\nprotein"           = "#ffb000",
  "Calmodulin"                         = "#fe6100",
  "Calcium-dependent\nprotein kinase"  = "#dc267f",
  "CDPK-related\nprotein"              = "#785ef0",
  "Unknown"                            = "grey80"
)

ha_cdpk <- rowAnnotation(
  df  = row_anno_cdpk,
  col = list(Annotation = anno_cols),
  show_annotation_name = FALSE,
  gp = gpar(col = "black", lwd = 0.6),
  annotation_legend_param = list(
    Annotation = list(
      title = NULL,
      labels_gp = legend_labels_gp
    )
  )
)

ht_cdpk <- Heatmap(
  mat_cdpk_scaled,
  name = "Row z-score\nlog2(expr + 1)",
  col  = hm_col_fun,
  cluster_rows = row_hc_cdpk,
  cluster_columns = FALSE,
  column_order = expr_cols,
  left_annotation = ha_cdpk,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_rot = 0,
  rect_gp = gpar(col = "black"),
  heatmap_legend_param = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp)
)

fig5C <- render_heatmap_png(
  ht_cdpk,
  filename = file.path(dir_figures, "fig5c_heatmap.png")
)

# ==============================
# Assemble multi-panel figure
# ==============================
fig5_multiplot <- ggdraw() +
  draw_plot(fig5A, x = 0.02, y = 0.30, width = 0.48, height = 0.70) +
  draw_plot(fig5B, x = 0.02, y = 0.00, width = 0.50, height = 0.28) +
  draw_plot(fig5C, x = 0.52, y = 0.15, width = 0.48, height = 0.70) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)"),
    size  = 16,
    x     = c(0, 0, 0.5),
    y     = c(1, 0.3, 0.85)
  )
fig5_multiplot

ggsave(
  filename = file.path(dir_figures, "Figure_5.png"),
  plot     = fig5_multiplot,
  width    = 26,
  height   = 18,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

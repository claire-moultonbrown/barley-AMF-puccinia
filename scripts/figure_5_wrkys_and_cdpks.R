##############################################################################

# Figure 5 – WRKYs and CDPKs

##############################################################################

# ==============================
# Setup
# ==============================
rm(list = ls())

dir_raw       <- "data_raw"
dir_figures   <- "figures"
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# ==============================
# Libraries
# ==============================
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(scales)
library(grid)

# ==============================
# Read data
# ==============================
wrky    <- readr::read_csv(file.path(dir_raw, "fig5a_transcriptomics_wrky_expression.csv"), show_col_types = FALSE)
wrky_go <- readr::read_csv(file.path(dir_raw, "fig5b_transcriptomics_wrky_go_terms.csv"),   show_col_types = FALSE)
cdpk    <- readr::read_csv(file.path(dir_raw, "fig5c_transcriptomics_cdpk_expression.csv"), show_col_types = FALSE)

# ==============================
# Global settings (shared across heatmaps to match legend sizes)
# ==============================
expr_cols  <- c("C", "A", "P", "AP")
labels_col <- c("C", "A", "P", "AP")

hm_cols   <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
hm_breaks <- c(-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2)

hm_legend_breaks <- c(-2, -1, 0, 1, 2)
hm_legend_labels <- c("-2", "-1", "0", "1", "2")

hm_fontsize     <- 9
hm_fontsize_col <- 9
hm_cellwidth    <- 30
hm_border_col   <- "black"

# Bubble palette (UNCHANGED)
bubble_gradient_cols <- c("blue", "yellow", "red")

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
# Helpers
# ==============================
prep_expr_matrix <- function(df, id_col, expr_cols) {
  df %>%
    dplyr::select(all_of(c(id_col, expr_cols))) %>%
    mutate(across(all_of(expr_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    column_to_rownames(id_col) %>%
    as.matrix() %>%
    (\(m) m[, expr_cols, drop = FALSE])()
}

normalise_unknown <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Na", "na", "N/A", "n/a", "Unknown", "unknown", "UNK", "unk")] <- "Unknown"
  x
}

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

# pheatmap -> gtable (silent), for cowplot placement and legend extraction
pheatmap_gtable <- function(mat,
                            annotation_row = NULL,
                            annotation_colors = NULL,
                            legend = TRUE,
                            annotation_legend = TRUE) {
  
  ph <- pheatmap::pheatmap(
    mat,
    color                    = hm_cols,
    breaks                   = hm_breaks,
    scale                    = "row",
    cluster_rows             = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method        = "ward.D2",
    cluster_cols             = FALSE,
    show_rownames            = TRUE,
    show_colnames            = TRUE,
    labels_col               = if (ncol(mat) == length(labels_col)) labels_col else colnames(mat),
    fontsize                 = hm_fontsize,
    fontsize_col             = hm_fontsize_col,
    cellwidth                = hm_cellwidth,
    border_color             = hm_border_col,
    legend                   = legend,
    annotation_row           = annotation_row,
    annotation_colors        = annotation_colors,
    annotation_names_row     = FALSE,
    annotation_legend        = annotation_legend,
    legend_breaks            = hm_legend_breaks,
    legend_labels            = hm_legend_labels,
    silent                   = TRUE
  )
  
  ph$gtable
}

extract_grob_by_name <- function(gt, name_pattern) {
  idx <- grep(name_pattern, gt$layout$name)
  if (length(idx) == 0) return(NULL)
  gt$grobs[[idx[1]]]
}

grob_to_gg <- function(g) {
  cowplot::ggdraw() + cowplot::draw_grob(g)
}

# ==============================
# Panel A – WRKY heatmap (Class then Role; drop Unknown)
# ==============================
mat_wrky <- prep_expr_matrix(wrky, id_col = "Name", expr_cols = expr_cols)

anno_wrky <- wrky %>%
  dplyr::select(Name, Role, Class) %>%
  distinct(Name, .keep_all = TRUE) %>%
  mutate(
    Role  = normalise_unknown(Role),
    Class = normalise_unknown(Class)
  ) %>%
  filter(Role != "Unknown", Class != "Unknown")

keep_wrky <- intersect(rownames(mat_wrky), anno_wrky$Name)
mat_wrky  <- mat_wrky[keep_wrky, , drop = FALSE]
anno_wrky <- anno_wrky %>% filter(Name %in% rownames(mat_wrky))

annotation_row_wrky <- data.frame(
  Class = factor(anno_wrky$Class, levels = c("I", "II", "III")),
  Role  = factor(anno_wrky$Role,  levels = c("Abiotic", "Biotic", "N/A")),
  row.names   = anno_wrky$Name,
  check.names = FALSE
)

# Force column order explicitly
annotation_row_wrky <- annotation_row_wrky[rownames(mat_wrky), c("Class", "Role"), drop = FALSE]

# Colours (your updated purples, and your Role colours)
class_cols <- c("I" = "#dadaeb", "II" = "#9e9ac8", "III" = "#6a51a3")
role_cols  <- c("Abiotic" = "yellow", "Biotic" = "#00d65c", "N/A" = "grey80")

ann_colors_wrky <- list(Class = class_cols, Role = role_cols)

gt_A <- pheatmap_gtable(
  mat = mat_wrky,
  annotation_row = annotation_row_wrky,
  annotation_colors = ann_colors_wrky,
  legend = TRUE,
  annotation_legend = TRUE
)
fig5A <- grob_to_gg(gt_A)

# ==============================
# Panel B – bubble plot (UNCHANGED)
# ==============================
wrky_go2 <- wrky_go %>%
  mutate(neglog10padj = -log10(p_value))

fig5B <- make_bubble(
  wrky_go2,
  y_col       = "GO_name",
  size_col    = "number",
  colour_col  = "neglog10padj",
  colour_name = "-log10(adj p-value)",
  size_name   = "Count",
  theme_obj   = theme_bubble
)

# ==============================
# Panel C – CDPK heatmap (annotation legend above; heatmap legend right; drop Unknown)
# ==============================
mat_cdpk <- prep_expr_matrix(cdpk, id_col = "Name", expr_cols = expr_cols)

anno_cdpk <- cdpk %>%
  dplyr::select(Name, Annotation) %>%
  distinct(Name, .keep_all = TRUE) %>%
  mutate(Annotation = normalise_unknown(Annotation)) %>%
  filter(Annotation != "Unknown")

keep_cdpk <- intersect(rownames(mat_cdpk), anno_cdpk$Name)
mat_cdpk  <- mat_cdpk[keep_cdpk, , drop = FALSE]
anno_cdpk <- anno_cdpk %>% filter(Name %in% rownames(mat_cdpk))

anno_cdpk$Annotation <- factor(
  anno_cdpk$Annotation,
  levels = c(
    "Calmodulin-like protein",
    "Calmodulin",
    "Calcium-dependent protein kinase",
    "CDPK-related protein"
  ),
  labels = c(
    "Calmodulin-like protein",
    "Calmodulin",
    "Calcium-dependent protein kinase",
    "CDPK-related protein"
  )
)

annotation_row_cdpk <- data.frame(
  Annotation = anno_cdpk$Annotation,
  row.names   = anno_cdpk$Name,
  check.names = FALSE
)
annotation_row_cdpk <- annotation_row_cdpk[rownames(mat_cdpk), , drop = FALSE]

cdpk_levels <- levels(droplevels(annotation_row_cdpk$Annotation))
cdpk_cols <- scales::hue_pal()(length(cdpk_levels))
names(cdpk_cols) <- cdpk_levels
ann_colors_cdpk <- list(Annotation = cdpk_cols)

# Annotation legend only (no heatmap legend)
gt_annleg <- pheatmap_gtable(
  mat = mat_cdpk,
  annotation_row = annotation_row_cdpk,
  annotation_colors = ann_colors_cdpk,
  legend = FALSE,
  annotation_legend = TRUE
)
ann_leg_grob <- extract_grob_by_name(gt_annleg, "annotation_legend")
if (is.null(ann_leg_grob)) ann_leg_grob <- extract_grob_by_name(gt_annleg, "legend")
fig5C_legend <- if (!is.null(ann_leg_grob)) grob_to_gg(ann_leg_grob) else cowplot::ggdraw()

# Heatmap body with heatmap legend ON; annotation legend OFF
gt_body <- pheatmap_gtable(
  mat = mat_cdpk,
  annotation_row = annotation_row_cdpk,
  annotation_colors = ann_colors_cdpk,
  legend = TRUE,
  annotation_legend = FALSE
)
fig5C_body <- grob_to_gg(gt_body)

# Stack legend above body; legend smaller
fig5C <- ggdraw()+
  draw_plot(fig5C_legend, x = 0.25, y = 0.75, width = 0.25, height = 0.25)+
  draw_plot(fig5C_body,   x = 0.25,    y = 0,    width = 0.5,  height = 0.85)

fig5C

# ==============================
# Assemble multi-panel figure
# ==============================
fig5_multiplot <- ggdraw() +
  draw_plot(fig5A, x = 0.02, y = 0.30, width = 0.48, height = 0.70) +
  draw_plot(fig5B, x = 0.02, y = 0.00, width = 0.50, height = 0.28) +
  draw_plot(fig5C, x = 0.55, y = 0.05, width = 0.42, height = 0.85) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)"),
    size  = 16,
    x     = c(0.00, 0.00, 0.6),
    y     = c(1.00, 0.30, 0.92)
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

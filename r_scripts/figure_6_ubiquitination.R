# ##############################################################################

# Figure 6 – Ubiquitination

# ##############################################################################

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
library(tidyverse)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(multcompView)

# ==============================
# Inputs and output filepaths
# ==============================
f_western_png <- file.path(dir_figures, "fig6a_ubiquitination_western.png")
f_western_csv <- file.path(dir_raw,     "fig6a_ubiquitination_western_data.csv")

f_expr_csv    <- file.path(dir_raw,     "fig6b_ubiquitin_expression.csv")
f_anno_csv    <- file.path(dir_raw,     "fig6b_ubiquitin_annotation.csv")

f_heatmap_png <- file.path(dir_figures, "fig6b_heatmap.png")
f_out_fig6    <- file.path(dir_figures, "Figure_6.png")

# ==============================
# Shared heatmap style (match Figures 4/5)
# ==============================
# Important: labels_col must match the column order you plot (value_cols below)
labels_col <- c("C", "P", "A", "AP")

hm_cols   <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
hm_breaks <- c(-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2)

hm_legend_breaks <- c(-2, -1, 0, 1, 2)
hm_legend_labels <- c("-2", "-1", "0", "1", "2")

hm_fontsize     <- 9
hm_fontsize_col <- 14   # column label size
hm_cellwidth    <- 30
hm_border_col   <- "black"

# ==============================
# Helpers
# ==============================
scale_rows <- function(x) {
  x_scaled <- t(scale(t(x)))
  x_scaled[is.nan(x_scaled)] <- NA
  x_scaled
}

get_main_pvals <- function(mod) {
  tab <- stats::anova(mod)
  rn <- trimws(rownames(tab))
  
  pick <- function(term) {
    i <- which(rn == term)
    if (length(i) != 1) return(NA_real_)
    as.numeric(tab[i, "Pr(>F)"])
  }
  
  list(
    p_amf = pick("AMF"),
    p_puc = pick("Puccinia"),
    p_int = pick("AMF:Puccinia")
  )
}

fmt_p_stars <- function(p, digits = 3) {
  if (is.na(p)) return("p = NA")
  stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.1) "." else ""
  p_txt <- formatC(p, format = "f", digits = digits)
  paste0("p = ", p_txt, stars)
}

anova_label_from_model <- function(mod) {
  p <- get_main_pvals(mod)
  paste0(
    "A: ",   fmt_p_stars(p$p_amf), "\n",
    "P: ",   fmt_p_stars(p$p_puc), "\n",
    "AxP: ", fmt_p_stars(p$p_int)
  )
}

# ---- pheatmap -> gtable (silent), so we can assemble with cowplot ----
pheatmap_gtable <- function(mat,
                            annotation_row = NULL,
                            annotation_colors = NULL,
                            legend = TRUE,
                            annotation_legend = TRUE,
                            cellwidth = hm_cellwidth,
                            fontsize = hm_fontsize,
                            fontsize_col = hm_fontsize_col,
                            border_color = hm_border_col) {
  
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
    # Use labels_col that matches the plotted order, so headers stay centred
    labels_col               = if (ncol(mat) == length(labels_col)) labels_col else colnames(mat),
    fontsize                 = fontsize,
    fontsize_col             = fontsize_col,
    cellwidth                = cellwidth,
    border_color             = border_color,
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

# ---- Panel B: pheatmap heatmap -> PNG -> cowplot panel ----
make_pheatmap_panel <- function(expr_df,
                                anno_df,
                                name_col = "Name",
                                value_cols = c("C", "P", "A", "AP"),
                                anno_cols = c("Function", "Type"),
                                anno_palettes = list(),
                                outfile_png,
                                width_px = 2200,
                                height_px = 2600,
                                res_dpi = 300) {
  
  # expression matrix
  mat <- expr_df %>%
    dplyr::select(dplyr::all_of(c(name_col, value_cols))) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(value_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::distinct(.data[[name_col]], .keep_all = TRUE) %>%
    tibble::column_to_rownames(name_col) %>%
    as.matrix()
  
  # enforce the desired plotted order
  mat <- mat[, value_cols, drop = FALSE]
  
  # row annotation dataframe
  row_anno <- anno_df %>%
    dplyr::select(dplyr::all_of(c(name_col, anno_cols))) %>%
    dplyr::distinct(.data[[name_col]], .keep_all = TRUE) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(anno_cols), ~ ifelse(is.na(.x) | .x == "" | .x == "NA", "Unknown", .x))) %>%
    tibble::column_to_rownames(name_col)
  
  # match order to matrix rows (and drop anything missing)
  row_anno <- row_anno[rownames(mat), , drop = FALSE]
  
  # coerce to factors with stable level order
  for (ac in anno_cols) row_anno[[ac]] <- factor(row_anno[[ac]])
  
  # ensure palettes cover all levels
  ann_colors <- list()
  for (ac in anno_cols) {
    levs <- levels(row_anno[[ac]])
    pal  <- anno_palettes[[ac]]
    
    if (is.null(pal)) {
      pal <- setNames(grDevices::hcl.colors(length(levs), "Set 2"), levs)
    }
    if ("Unknown" %in% levs && !("Unknown" %in% names(pal))) {
      pal <- c(pal, Unknown = "grey80")
    }
    pal <- pal[levs]
    if (any(is.na(pal))) {
      miss <- levs[is.na(pal)]
      stop("Palette for '", ac, "' is missing levels: ", paste(miss, collapse = ", "))
    }
    ann_colors[[ac]] <- pal
  }
  
  # build gtable then draw to PNG
  gt <- pheatmap_gtable(
    mat = mat,
    annotation_row = row_anno,
    annotation_colors = ann_colors,
    legend = TRUE,
    annotation_legend = TRUE
  )
  
  png(outfile_png, width = width_px, height = height_px, res = res_dpi)
  grid::grid.newpage()
  grid::grid.draw(gt)
  dev.off()
  
  cowplot::ggdraw() + cowplot::draw_image(outfile_png, scale = 1)
}

# ==============================
# Panel styling (A)
# ==============================
treatment_colours <- c(
  "Control"  = "white",
  "AMF"      = "#AEECEF",
  "Puccinia" = "#6D9DC5",
  "Both"     = "#068D9D"
)

theme_claire <- theme(
  panel.background = element_blank(),
  legend.title     = element_blank(),
  axis.title.x     = element_blank(),
  axis.text.y      = element_text(size = 7,  face = "bold", colour = "black"),
  axis.text.x      = element_text(size = 10, face = "bold", colour = "black"),
  axis.title.y     = element_text(size = 10, face = "bold"),
  axis.line        = element_line(colour = "black", linewidth = 0.5),
  legend.key       = element_blank(),
  aspect.ratio     = 1
)

# ##############################################################################
# PANEL 6A – Western blot + quantification
# ##############################################################################
ubiquitin <- read.csv(f_western_csv, stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Control", "Puccinia", "AMF", "Both")),
    Treatment_initial = factor(Treatment_initial, levels = c("C", "P", "A", "AP")),
    AMF      = factor(AMF),
    Puccinia = factor(Puccinia)
  )

control_mean <- ubiquitin %>%
  dplyr::filter(Treatment == "Control") %>%
  dplyr::summarise(mu = mean(ubiquitin.actin, na.rm = TRUE)) %>%
  dplyr::pull(mu)

ubiquitin <- ubiquitin %>%
  dplyr::mutate(
    fc      = ubiquitin.actin / control_mean,
    log2_fc = log2(fc)
  )

#Ubiquitin ANOVA
ubiquitin_aov     <- aov(log2_fc ~ AMF * Puccinia, data = ubiquitin)
anova_label_6A    <- anova_label_from_model(ubiquitin_aov)

#Tukey HSD
# one-way model for letters
ubiquitin_aov_trt <- aov(log2_fc ~ Treatment, data = ubiquitin)
tuk <- TukeyHSD(ubiquitin_aov_trt, "Treatment")
cld <- multcompView::multcompLetters4(ubiquitin_aov_trt, tuk)
letters_df <- tibble::tibble(
  Treatment  = names(cld$Treatment$Letters),
  sig_letter = unname(cld$Treatment$Letters)
)

#make summary dataframe for plotting ubiquitin signal inc. Tukey HSD letters
ub_sum <- ubiquitin %>%
  group_by(Treatment, Treatment_initial) %>%
  summarise(
    mean = mean(fc, na.rm = TRUE),
    n    = sum(!is.na(fc)),
    sd   = sd(fc, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "Treatment")

#Add western blot image
fig6A_upper <- ggdraw() + draw_image(f_western_png, scale = 1)

#Draw bar chart of ubiquitin signal from Western blot with mean and SE plotted
fig6A_lower <- ggplot(ub_sum, aes(x = Treatment_initial, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7, linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.5) +
  geom_text(aes(y = mean + se + 0.25, label = sig_letter), size = 3) +
  annotate("text", x = 3.65, y = Inf, label = anova_label_6A, size = 2.5, vjust = 0.99) +
  ylab("Ubiquitin/Actin\n(fold change vs Control)") +
  theme_claire +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  scale_fill_manual(values = treatment_colours) +
  theme(legend.position = "none")

fig6A <- ggdraw() +
  draw_plot(fig6A_upper, x = 0.00, y = 0.58, width = 1.00, height = 0.42) +
  draw_plot(fig6A_lower, x = 0.20, y = 0.00, width = 0.50, height = 0.70)

fig6A

# ##############################################################################
# PANEL 6B – Heatmap
# ##############################################################################
expr <- read.csv(f_expr_csv, stringsAsFactors = FALSE, check.names = FALSE)
anno <- read.csv(f_anno_csv, stringsAsFactors = FALSE, check.names = FALSE)

# keep labels_col consistent with plotted order
labels_col <- c("C", "P", "A", "AP")

function_cols <- c(
  "E1"      = "#ff9289",
  "E2"      = "#d9abff",
  "E3"      = "#d5c91a",
  "Unknown" = "grey80"
)

type_cols <- c(
  "ATG"     = "#feb221",
  "RBR"     = "#18e583",
  "RING"    = "#40cdff",
  "RING-H2" = "#ff98ff",
  "U-box"   = "#87d91c",
  "UBA"     = "#19f0ff",
  "UBC"     = "#1be8cd",
  "UEV"     = "#ea7ecc",
  "Unknown" = "grey80"
)

fig6B <- make_pheatmap_panel(
  expr_df        = expr,
  anno_df        = anno,
  name_col       = "Name",
  value_cols     = c("C", "P", "A", "AP"),
  anno_cols      = c("Function", "Type"),
  anno_palettes  = list(Function = function_cols, Type = type_cols),
  outfile_png    = f_heatmap_png
)
fig6B

# ##############################################################################
# Combine Figure 6 (A + B) and save
# ##############################################################################
fig6_multiplot <- ggdraw() +
  draw_plot(fig6A, x = 0.02, y = 0.00, width = 0.5, height = 1) +
  draw_plot(fig6B, x = 0.50, y = 0.00, width = 0.6, height = 1) +
  draw_plot_label(
    label = c("(A)", "(B)"),
    size  = 16,
    x     = c(0.00, 0.50),
    y     = c(1.00, 1.00)
  )

ggsave(
  filename = f_out_fig6,
  plot     = fig6_multiplot,
  width    = 22,
  height   = 16,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

fig6_multiplot

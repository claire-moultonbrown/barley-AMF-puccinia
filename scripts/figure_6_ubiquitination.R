# ##############################################################################
# Figure 6 – Ubiquitination
# ##############################################################################

# ==============================
# Setup
# ==============================
rm(list = ls())

dir_raw       <- "data_raw"
dir_figures   <- "figures"
dir.create(dir_figures,   showWarnings = FALSE, recursive = TRUE)

# ==============================
# Libraries
# ==============================
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(circlize)   # colorRamp2()
library(grid)       # gpar(), unit()

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
# Helpers
# ==============================
scale_rows <- function(x) {
  x_scaled <- t(scale(t(x)))
  x_scaled[is.nan(x_scaled)] <- NA
  x_scaled
}

get_main_pvals <- function(mod) {
  # Works reliably for aov/lm objects
  tab <- stats::anova(mod)

  # Make rownames consistent for matching
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

# ==============================
# Heatmap panel function (ComplexHeatmap -> PNG -> cowplot panel)
# ==============================
make_heatmap_panel <- function(expr_df,
                               anno_df,
                               name_col = "Name",
                               value_cols = c("C", "P", "A", "AP"),
                               anno_cols = c("Function", "Type"),
                               row_split_k = 4,
                               outfile_png,
                               heatmap_name = "Row z-score\nlog2(expr + 1)",
                               col_fun = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027")),
                               anno_palettes = list(),
                               show_row_names = TRUE,
                               row_name_size = 7,
                               col_name_size = 14,
                               col_name_face = "bold",
                               rect_border = "black",
                               rect_lwd = 0.5,
                               row_gap_mm = 2,
                               width_px = 2200,
                               height_px = 2600,
                               res_dpi = 300) {
  
  # ---- expression matrix ----
  mat <- expr_df %>%
    dplyr::select(dplyr::all_of(c(name_col, value_cols))) %>%
    mutate(across(all_of(value_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    distinct(.data[[name_col]], .keep_all = TRUE) %>%
    column_to_rownames(name_col) %>%
    as.matrix()
  
  mat <- mat[, value_cols, drop = FALSE]
  
  # ---- row annotation dataframe ----
  row_anno <- anno_df %>%
    dplyr::select(dplyr::all_of(c(name_col, anno_cols))) %>%
    distinct(.data[[name_col]], .keep_all = TRUE) %>%
    mutate(across(all_of(anno_cols), ~ ifelse(is.na(.x) | .x == "", "Unknown", .x))) %>%
    column_to_rownames(name_col)
  
  # match order to matrix rows
  row_anno <- row_anno[rownames(mat), , drop = FALSE]
  
  # coerce to factors
  for (ac in anno_cols) row_anno[[ac]] <- factor(row_anno[[ac]])
  
  # ---- scale ----
  mat_scaled <- scale_rows(log2(mat + 1))
  
  # ---- clustering ----
  row_hc <- hclust(dist(mat_scaled))
  
  # ---- palettes: ensure every level is mapped ----
  # Build a col list for rowAnnotation()
  anno_col_list <- list()
  for (ac in anno_cols) {
    levs <- levels(row_anno[[ac]])
    
    pal <- anno_palettes[[ac]]
    if (is.null(pal)) {
      # fallback palette if user didn't provide one (not ideal, but safe)
      pal <- setNames(grDevices::hcl.colors(length(levs), "Set 2"), levs)
    }
    
    # ensure Unknown exists
    if (!("Unknown" %in% names(pal)) && ("Unknown" %in% levs)) {
      pal <- c(pal, Unknown = "grey80")
    }
    
    # reorder + check missing
    pal <- pal[levs]
    if (any(is.na(pal))) {
      missing_levels <- levs[is.na(pal)]
      stop("Palette for '", ac, "' is missing levels: ", paste(missing_levels, collapse = ", "))
    }
    anno_col_list[[ac]] <- pal
  }
  
  ha <- rowAnnotation(
    df  = row_anno,
    col = anno_col_list,
    annotation_name_side = "bottom",
    annotation_name_gp   = gpar(fontface = "bold", fontsize = 10)
  )
  
  ht <- Heatmap(
    mat_scaled,
    name = heatmap_name,
    col  = col_fun,
    cluster_rows    = row_hc,
    cluster_columns = FALSE,
    column_order    = value_cols,
    row_split       = row_split_k,
    row_gap         = unit(row_gap_mm, "mm"),
    left_annotation = ha,
    rect_gp         = gpar(col = rect_border, lwd = rect_lwd),
    show_row_names  = show_row_names,
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = row_name_size),
    column_names_gp = gpar(fontsize = col_name_size, fontface = col_name_face),
    column_names_rot = 0
  )
  
  png(outfile_png, width = width_px, height = height_px, res = res_dpi)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  ggdraw() + draw_image(outfile_png, scale = 1)
}

# ==============================
# Panel styling
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
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Puccinia", "AMF", "Both")),
    Treatment_initial = factor(Treatment_initial, levels = c("C", "P", "A", "AP")),
    AMF      = factor(AMF),
    Puccinia = factor(Puccinia)
  )

control_mean <- ubiquitin %>%
  filter(Treatment == "Control") %>%
  summarise(mu = mean(ubiquitin.actin, na.rm = TRUE)) %>%
  pull(mu)

ubiquitin <- ubiquitin %>%
  mutate(
    fc      = ubiquitin.actin / control_mean,
    log2_fc = log2(fc)
  )

# ANOVA label generated from model (no hard-coding)
ubiquitin_aov   <- aov(log2_fc ~ AMF * Puccinia, data = ubiquitin)
anova_label_6A  <- anova_label_from_model(ubiquitin_aov)

ub_sum <- ubiquitin %>%
  group_by(Treatment, Treatment_initial) %>%
  summarise(
    mean = mean(fc, na.rm = TRUE),
    n    = sum(!is.na(fc)),
    sd   = sd(fc, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(sig_letter = c("bc", "a", "c", "ab"))

fig6A_upper <- ggdraw() + draw_image(f_western_png, scale = 1)

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

panel_6A <- ggdraw() +
  draw_plot(fig6A_upper, x = 0.00, y = 0.58, width = 1.00, height = 0.42) +
  draw_plot(fig6A_lower, x = 0.20, y = 0.00, width = 0.50, height = 0.70)
panel_6A

# ##############################################################################
# PANEL 6B – Heatmap (wrapped in function)
# ##############################################################################

expr <- read.csv(f_expr_csv, stringsAsFactors = FALSE, check.names = FALSE)
anno <- read.csv(f_anno_csv, stringsAsFactors = FALSE, check.names = FALSE)

function_cols <- c(
  "E1"      = "#0072B2",
  "E2"      = "#E69F00",
  "E3"      = "#009E73",
  "Unknown" = "grey80"
)

type_cols <- c(
  "ATG"     = "#1B9E77",
  "RBR"     = "#D95F02",
  "RING"    = "#7570B3",
  "RING-H2" = "#E7298A",
  "U-box"   = "#66A61E",
  "UBA"     = "#E6AB02",
  "UBC"     = "#A6761D",
  "UEV"     = "#A6CEE3",
  "Unknown" = "grey80"
)

panel_6B <- make_heatmap_panel(
  expr_df       = expr,
  anno_df       = anno,
  name_col      = "Name",
  value_cols    = c("C", "P", "A", "AP"),       # keep your current order
  anno_cols     = c("Function", "Type"),
  row_split_k   = 4,
  outfile_png   = f_heatmap_png,
  anno_palettes = list(Function = function_cols, Type = type_cols),
  show_row_names = TRUE
)

# ##############################################################################
# Combine Figure 6 (A + B) and save
# ##############################################################################

fig6_multiplot <- plot_grid(
  panel_6A, panel_6B,
  ncol = 2,
  rel_widths = c(1, 1.2),
  labels = c("(A)", "(B)"),
  label_size = 18,
  label_x = c(0.02, 0.00),
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

ggsave(
  filename = f_out_fig6,
  plot     = fig6_multiplot,
  width    = 26,
  height   = 16,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

fig6_multiplot

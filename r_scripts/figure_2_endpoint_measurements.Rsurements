# ##############################################################################

# Figure 2 - Endpoint measurements of AMF colonisation, Puccinia infection,
# and plant biomass

# ##############################################################################

# ==============================
# Setup
# ==============================
rm(list = ls())

dir_raw     <- "data_raw"
dir_figures <- "figures"
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# ==============================
# Libraries
# ==============================
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(cowplot)
library(emmeans)
library(multcomp)

# ==============================
# Read data
# ==============================
roots_pic <- file.path(dir_figures, "fig2a_root_annotated.png")

barley <- readr::read_csv(file.path(dir_raw, "fig2_endpoint_measurements.csv"),
                          show_col_types = FALSE) %>%
  mutate(
    across(c(AMF, Puccinia, Treatment, Measurement), as.factor),
    Treatment = factor(Treatment, levels = treat_levels)
  )

# ==============================
# Global settings
# ==============================
treat_levels <- c("Control", "AMF", "Puccinia", "Both")
treat_labels <- c("C", "A", "P", "AP")
names(treat_labels) <- treat_levels

treatment_colours <- c(
  "Control"  = "white",
  "AMF"      = "#AEECEF",
  "Puccinia" = "#6D9DC5",
  "Both"     = "#068D9D"
)

theme_claire <- theme(
  panel.background = element_rect(fill = NA),
  legend.title     = element_blank(),
  axis.title.x     = element_blank(),
  axis.text.y      = element_text(size = 13, face = "bold", colour = "black"),
  axis.text.x      = element_text(size = 18, face = "bold", colour = "black"),
  axis.title.y     = element_text(size = 18, face = "bold"),
  legend.text      = element_text(size = 22, face = "bold"),
  panel.border     = element_blank(),
  axis.line        = element_line(colour = "black", linewidth = 1),
  legend.key       = element_blank(),
  aspect.ratio     = 1
)

# ==============================
# Helper functions
# ==============================

# ==============================
# Helper: nice p-value formatting

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3g", p)
}

# ==============================
# Helper: get p-value stars

fmt_p_stars <- function(p) {
  if (is.na(p)) return("p = NA")
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05,  "*", "")))
  paste0("p = ", formatC(p, format = "g", digits = 3), stars)
}

# ==============================
# Helper: get p-values
get_main_pvals <- function(mod) {
  tab <- anova(mod)
  
  getp <- function(term) {
    i <- which(rownames(tab) == term)
    if (length(i) != 1) return(NA_real_)
    tab[i, "Pr(>F)"]
  }
  
  list(
    p_amf = getp("AMF"),
    p_puc = getp("Puccinia"),
    p_int = getp("AMF:Puccinia")
  )
}


# ==============================
# Helper: build one bar panel from raw data
# - runs ANOVA (Value ~ AMF * Puccinia)
# - extracts Tukey letters for each Treatment (AMF*Puccinia combo)
# - plots mean +- SE with letters + p-values

make_endpoint_panel <- function(df, ylab_text, y_limits) {
  stopifnot(all(c("Value", "AMF", "Puccinia", "Treatment") %in% names(df)))
  
  # Summary stats for bars
  df_sum <- df %>%
    group_by(Treatment) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      n    = sum(!is.na(Value)),
      sd   = sd(Value, na.rm = TRUE),
      se   = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      Treatment = factor(Treatment, levels = treat_levels),
      Treatment_initial = unname(treat_labels[as.character(Treatment)])
    )
  
  # ANOVA
  aov_fit <- aov(Value ~ AMF * Puccinia, data = df)
  
  pvals <- get_main_pvals(aov_fit)
  
  p_amf <- pvals$p_amf
  p_puc <- pvals$p_puc
  p_int <- pvals$p_int
  
  stats_label <- paste0(
    "A: ",   fmt_p_stars(p_amf), "\n",
    "P: ",   fmt_p_stars(p_puc), "\n",
    "AxP: ", fmt_p_stars(p_int)
  )
  
  # Tukey letters on the AMF*Puccinia combinations
  em <- emmeans(aov_fit, ~ AMF * Puccinia)
  cld_tbl <- multcomp::cld(em, Letters = letters, decreasing = TRUE) %>%
    as.data.frame()
  
  # Map letters back onto Treatment order via AMF/Puccinia pattern
  # Assumes Treatment encodes the 4 combinations (Control, AMF, Puccinia, Both)
  # Control: AMF=0, Puccinia=0
  # AMF:     AMF=1, Puccinia=0
  # Puccinia:AMF=0, Puccinia=1
  # Both:    AMF=1, Puccinia=1
  # If your factor labels differ, this still works because it matches by the factor levels in df.
  combo_key <- df %>%
    distinct(Treatment, AMF, Puccinia) %>%
    mutate(Treatment = factor(Treatment, levels = treat_levels))
  
  cld_map <- combo_key %>%
    left_join(
      cld_tbl %>% transmute(AMF, Puccinia, sig_letter = .group),
      by = c("AMF", "Puccinia")
    ) %>%
    arrange(Treatment) %>%
    mutate(sig_letter = gsub(" ", "", sig_letter))
  
  df_sum <- df_sum %>%
    left_join(cld_map %>% dplyr::select(Treatment, sig_letter), by = "Treatment")
  
  # Plot
  y_pad <- (y_limits[2] - y_limits[1]) * 0.05
  letter_y <- df_sum$mean + df_sum$se + y_pad
  
  ggplot(df_sum, aes(x = Treatment_initial, y = mean, fill = Treatment)) +
    geom_col(colour = "black", width = 0.7, linewidth = 1) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
    geom_text(aes(y = letter_y, label = sig_letter), size = 5) +
    annotate("text",
             x = 2.5,
             y = y_limits[2] * 0.88,
             label = stats_label,
             size = 3.5,
             hjust = 0
    ) +
    ylab(ylab_text) +
    scale_y_continuous(expand = c(0, 0), limits = y_limits) +
    scale_fill_manual(values = treatment_colours) +
    theme_claire +
    theme(legend.position = "none")
}

# ==============================
# Figure 2A: AMF stained root image panel
# ==============================

fig2A <- ggdraw() + draw_image(roots_pic, scale = 1)

# ==============================
# Figure 2B–D panels (B=AMF abundance, C=Puccinia abundance, D=Plant biomass)
# ==============================
amf_df     <- barley %>% filter(Measurement == "AMF")
pucc_df    <- barley %>% filter(Measurement == "PhEF2")
biomass_df <- barley %>% filter(Measurement == "Mass")

fig2B <- make_endpoint_panel(amf_df,
                             ylab_text = "AMF\nabundance",
                             y_limits = c(0, 25))
fig2B

fig2C <- make_endpoint_panel(pucc_df,
                             ylab_text = "Puccinia\nabundance",
                             y_limits = c(0, 35))
fig2C

fig2D <- make_endpoint_panel(biomass_df, 
                             ylab_text = "Plant\nbiomass (g)",  
                             y_limits = c(0, 3.1))
fig2D


fig2C

# ==============================
# Legend (single object)
# ==============================
leg_plot <- ggplot(barley %>% filter(Measurement == "AMF"),
                   aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot() +
  theme_claire +
  scale_fill_manual(values = treatment_colours) +
  theme(legend.position = "right")
leg <- get_legend(leg_plot)

# ==============================
# Assemble multi-panel figure
# ==============================
fig2_multiplot <- ggdraw() +
  draw_plot(fig2A, x = 0.10, y = 0.50, width = 0.46, height = 0.48) +
  draw_plot(fig2B, x = 0.00, y = 0.00, width = 0.33, height = 0.48) +
  draw_plot(fig2C, x = 0.33, y = 0.00, width = 0.33, height = 0.48) +
  draw_plot(fig2D, x = 0.66, y = 0.00, width = 0.33, height = 0.48) +
  draw_plot(leg,   x = 0.50, y = 0.50, width = 0.50, height = 0.50) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)", "(D)"),
    size  = 20,
    x     = c(0, 0, 0.33, 0.66),
    y     = c(1, 0.5, 0.5, 0.5)
  )
fig2_multiplot

# ==============================
# Save
# ==============================
ggsave(
  filename = file.path(dir_figures, "Figure_2.png"),
  plot     = fig2_multiplot,
  width    = 28,
  height   = 15,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

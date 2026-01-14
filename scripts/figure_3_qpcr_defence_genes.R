# ##############################################################################

# Figure 3 – qPCR of defence genes

# ##############################################################################

rm(list = ls())

# ==============================
# Directories
# ==============================
dir_raw       <- "data_raw"
dir_figures   <- "figures"
dir.create(dir_figures,   showWarnings = FALSE, recursive = TRUE)

# ==============================
# Libraries
# ==============================
library(tidyverse)
library(cowplot)
library(emmeans)
library(multcomp)

# ==============================
# Read + tidy data
# ==============================
barley <- read.csv(file.path(dir_raw, "fig3_qpcr_defence_genes.csv"))

barley <- barley %>%
  mutate(
    across(c(AMF, Puccinia, Treatment, Gene), as.factor),
    Treatment = factor(Treatment, levels = c("Control", "AMF", "Puccinia", "Both")),
    Treatment_initial = recode(Treatment,
                               "Control"  = "C",
                               "AMF"      = "A",
                               "Puccinia" = "P",
                               "Both"     = "AP"
    ),
    Treatment_initial = factor(Treatment_initial, levels = c("C", "A", "P", "AP"))
  )

# ==============================
# Helpers
# ==============================

#--------------------
# p stars function
fmt_p_stars <- function(p) {
  if (is.na(p)) return("p = NA")
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", "")))
  # scientific for tiny p, normal otherwise
  p_txt <- if (p < 0.001) format(p, scientific = TRUE, digits = 2) else format(round(p, 3), nsmall = 3)
  paste0("p = ", p_txt, stars)
}

#--------------------
# p values function
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

#--------------------
# Per-gene stats + letters function (letters up a bit, ANOVA down so it stays inside)
get_gene_stats <- function(df_gene, letter_frac = 0.15, anova_frac = 0.90) {
  
  mod <- aov(Value ~ AMF * Puccinia, data = df_gene)
  
  # Tukey on the 4 combinations
  em  <- emmeans(mod, ~ AMF * Puccinia)
  cld <- cld(em, Letters = letters, decreasing = TRUE) %>%
    as.data.frame() %>%
    mutate(
      Treatment_initial = case_when(
        AMF == "AM-" & Puccinia == "Puccinia-" ~ "C",
        AMF == "AM+" & Puccinia == "Puccinia-" ~ "A",
        AMF == "AM-" & Puccinia == "Puccinia+" ~ "P",
        AMF == "AM+" & Puccinia == "Puccinia+" ~ "AP",
        TRUE ~ NA_character_
      ),
      Treatment_initial = factor(Treatment_initial, levels = c("C", "A", "P", "AP")),
      sig_letter = gsub(" ", "", .group)
    )
  
  # ANOVA label (once per gene)
  pvals <- get_main_pvals(mod)
  anova_label <- paste0(
    "A: ",   fmt_p_stars(pvals$p_amf), "\n",
    "P: ",   fmt_p_stars(pvals$p_puc), "\n",
    "AxP: ", fmt_p_stars(pvals$p_int)
  )
  
  gene_name <- unique(df_gene$Gene)[1]
  
  ymin_panel  <- min(df_gene$Value, na.rm = TRUE)
  ymax_panel  <- max(df_gene$Value, na.rm = TRUE)
  panel_range <- ymax_panel - ymin_panel
  if (!is.finite(panel_range) || panel_range == 0) panel_range <- 1
  
  # Letters: a bit higher above the upper box (Q3)
  letters_df <- df_gene %>%
    group_by(Treatment_initial) %>%
    summarise(
      q3 = unname(quantile(Value, 0.75, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      y_letter = q3 + letter_frac * panel_range,
      Gene     = gene_name
    ) %>%
    left_join(
      cld %>%
        transmute(
          Treatment_initial,
          sig_letter
        ),
      by = "Treatment_initial"
    )
  
  # ANOVA text: place within the panel using a fraction of the range
  # anova_frac = 0.90 means 90% of the way from ymin to ymax (inside the plot)
  anova_df <- tibble(
    Gene      = gene_name,
    x_anova   = 1.25,
    y_anova   = ymin_panel + anova_frac * panel_range,
    anova_txt = anova_label
  )
  
  list(letters = letters_df, anova = anova_df)
}


# ==============================
# Compute letters + ANOVA labels (per gene)
# ==============================
stats_list <- barley %>%
  group_split(Gene) %>%
  map(get_gene_stats)

letters_df <- map_dfr(stats_list, "letters")
anova_df   <- map_dfr(stats_list, "anova")

# ==============================
# Plot parameters
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
  axis.text.y      = element_text(size = 13, face = "bold", colour = "black"),
  axis.text.x      = element_text(size = 13, face = "bold", colour = "black"),
  strip.text       = element_text(size = 16, face = "bold"),
  axis.title.y     = element_text(size = 18, face = "bold"),
  legend.text      = element_text(size = 22, face = "bold"),
  axis.line        = element_line(colour = "black", linewidth = 1),
  legend.key       = element_blank(),
  aspect.ratio     = 1
)

# ==============================
# Figure 3
# ==============================
fig3 <- ggplot(barley, aes(x = Treatment_initial, y = Value, fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(Gene ~ ., scales = "free_y") +
  ylab("Relative expression (2^-ΔCt)") +
  scale_fill_manual(values = treatment_colours) +
  theme_claire +
  theme(legend.position = "none") +
  geom_text(
    data = letters_df,
    aes(x = Treatment_initial, y = y_letter, label = sig_letter),
    inherit.aes = FALSE,
    size = 5
  ) +
  geom_text(
    data = anova_df,
    aes(x = x_anova, y = y_anova, label = anova_txt),
    inherit.aes = FALSE,
    size = 3.5,
    hjust = 0
  )

# ==============================
# Panel labels
# ==============================
fig3_labelled <- ggdraw() +
  draw_plot(fig3) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)", "(D)"),
    size  = 20,
    x     = c(0.1, 0.55, 0.1, 0.55),
    y     = c(1, 1, 0.52, 0.52)
  )

fig3_labelled

# ==============================
# Save
# ==============================
ggsave(
  filename = file.path(dir_figures, "Figure_3.png"),
  plot     = fig3_labelled,
  width    = 15,
  height   = 15,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

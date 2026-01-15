# ##############################################################################
# Figure 3 – qPCR of defence genes (means ± SE)
# ##############################################################################

rm(list = ls())

# ==============================
# Directories
# ==============================
dir_raw     <- "data_raw"
dir_figures <- "figures"
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# ==============================
# Libraries
# ==============================
library(cowplot)
library(emmeans)
library(multcomp)
library(tidyverse)

# ==============================
# Read + tidy data
# ==============================
barley <- read.csv(file.path(dir_raw, "fig3_qpcr_defence_genes.csv"))

barley <- barley %>%
  mutate(
    across(c(AMF, Puccinia, Treatment, Gene), as.factor),
    AMF = factor(AMF, levels = c("AM-", "AM+")),
    Puccinia = factor(Puccinia, levels = c("Puccinia-", "Puccinia+")),
    Treatment = factor(Treatment, levels = c("Control", "AMF", "Puccinia", "Both")),
    Treatment_initial = recode(
      Treatment,
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

# p-value + stars
p_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  p_txt <- if (p < 0.001) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    format(round(p, 3), nsmall = 3)
  }
  paste0("p = ", p_txt, p_stars(p))
}

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
# Per-gene stats (ANOVA + letters + labels)
# ==============================
get_gene_stats <- function(df_gene,
                           response = "Value",
                           letter_frac = 0.08,
                           anova_frac  = 0.65,   # same relative height across facets
                           adjust_method = "sidak") {
  
  mod <- aov(
    reformulate(c("AMF", "Puccinia", "AMF:Puccinia"),
                response = response),
    data = df_gene
  )
  
  # Panel range (keep yrng for letter spacing)
  y <- df_gene[[response]]
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  yrng <- ymax - ymin
  if (!is.finite(yrng) || yrng == 0) yrng <- 1
  if (!is.finite(ymax) || ymax == 0) ymax <- 1
  
  # ANOVA label
  pvals <- get_main_pvals(mod)
  anova_txt <- paste0(
    "A: ",   fmt_p(pvals$p_amf), "\n",
    "P: ",   fmt_p(pvals$p_puc), "\n",
    "AxP: ", fmt_p(pvals$p_int)
  )
  
  # ANOVA label y-position: aligned across facets (relative to ymax)
  anova_df <- tibble(
    Gene      = unique(df_gene$Gene),
    x_anova   = 1.25,
    y_anova   = Inf,
    anova_txt = anova_txt
  )
  
  
  # Letters (CLD)
  em <- emmeans(mod, ~ AMF * Puccinia)
  
  cld_tbl <- cld(
    em,
    Letters = letters,
    adjust  = adjust_method
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      Treatment_initial = case_when(
        AMF == "AM-" & Puccinia == "Puccinia-" ~ "C",
        AMF == "AM+" & Puccinia == "Puccinia-" ~ "A",
        AMF == "AM-" & Puccinia == "Puccinia+" ~ "P",
        AMF == "AM+" & Puccinia == "Puccinia+" ~ "AP"
      ),
      Treatment_initial = factor(Treatment_initial, levels = c("C", "A", "P", "AP")),
      sig_letter = gsub(" ", "", .group)
    ) %>%
    dplyr::select(Treatment_initial, sig_letter)
  
  # Letter y-position: above mean ± SE, with a small buffer scaled by panel range
  letters_df <- df_gene %>%
    group_by(Treatment_initial) %>%
    summarise(
      mean = mean(.data[[response]], na.rm = TRUE),
      se   = sd(.data[[response]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      Gene     = unique(df_gene$Gene),
      y_letter = mean + se + letter_frac * yrng
    ) %>%
    left_join(cld_tbl, by = "Treatment_initial")
  
  list(letters = letters_df, anova = anova_df)
}


# ==============================
# Summary stats for plotting
# ==============================
summary_df <- barley %>%
  group_by(Gene, Treatment, Treatment_initial) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se   = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ==============================
# Letters + ANOVA per gene
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
  axis.line        = element_line(colour = "black", linewidth = 1),
  legend.key       = element_blank(),
  aspect.ratio     = 1
)

# ==============================
# Figure 3 – bar plot (mean ± SE)
# ==============================
fig3 <- ggplot(summary_df,
               aes(x = Treatment_initial, y = mean, fill = Treatment)) +
  geom_col(colour = "black", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  facet_wrap(Gene ~ ., scales = "free") +
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
    hjust = 0.4,
    vjust = 1.2 #move anova label down (bigger = further down)
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
    y     = c(0.99, 0.99, 0.5, 0.5)
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


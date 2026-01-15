# ##############################################################################

# Figure S1 – Biomass panels (Exp. 4): Root mass, Shoot mass, Root:shoot ratio

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
barley <- read.csv(file.path(dir_raw, "meta_data.csv"))

barley <- barley %>%
  mutate(
    across(c(Inoculum, Puccinia, Treatment), as.factor),
    Inoculum = factor(Inoculum, levels = c("AM-", "AM+")),
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

#--------------------
# Significance stars from p
p_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

#--------------------
# p-value formatting: decimals; scientific only for very small p
# Adds stars for significant p-values
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  p_txt <- if (p < 0.001) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    format(round(p, 3), nsmall = 3)
  }
  paste0("p = ", p_txt, p_stars(p))
}

#--------------------
# Extract main p-values from two-way ANOVA
get_main_pvals <- function(mod) {
  tab <- anova(mod)
  
  getp <- function(term) {
    i <- which(rownames(tab) == term)
    if (length(i) != 1) return(NA_real_)
    tab[i, "Pr(>F)"]
  }
  
  list(
    p_amf = getp("Inoculum"),
    p_puc = getp("Puccinia"),
    p_int = getp("Inoculum:Puccinia")
  )
}

#--------------------
# Upper whisker (upper adjacent value) for a boxplot:
# max(x[x <= Q3 + 1.5*IQR])
upper_whisker <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  q3  <- unname(stats::quantile(x, 0.75, na.rm = TRUE))
  iqr <- unname(stats::IQR(x, na.rm = TRUE))
  cap <- q3 + 1.5 * iqr
  max(x[x <= cap], na.rm = TRUE)
}

#--------------------
# Stats for a given response column:
# - two-way ANOVA
# - letters from Inoculum x Puccinia combinations mapped to C/A/P/AP
# - positions for letters + ANOVA text inside panel
# Letters are placed just above the upper whisker per Treatment.
get_panel_stats <- function(df, response_col,
                            letter_frac = 0.04,
                            anova_frac = 0.90,
                            adjust_method = "sidak") {
  
  # Two-way factorial ANOVA
  mod <- aov(
    reformulate(c("Inoculum", "Puccinia", "Inoculum:Puccinia"), response = response_col),
    data = df
  )
  
  # Panel range for placement
  y <- df[[response_col]]
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  yrng <- ymax - ymin
  if (!is.finite(yrng) || yrng == 0) yrng <- 1
  
  # ANOVA label (p-values + stars)
  pvals <- get_main_pvals(mod)
  anova_label <- paste0(
    "A: ",   fmt_p(pvals$p_amf), "\n",
    "P: ",   fmt_p(pvals$p_puc), "\n",
    "AxP: ", fmt_p(pvals$p_int)
  )
  
  anova_df <- tibble(
    Panel     = response_col,
    x_anova   = 3.25,
    y_anova   = ymin + anova_frac * yrng,
    anova_txt = anova_label
  )
  
  # Letters from factorial combinations (Inoculum x Puccinia), mapped to C/A/P/AP
  em  <- emmeans(mod, ~ Inoculum * Puccinia)
  
  cld_tbl <- cld(
    em,
    Letters    = letters,
    decreasing = TRUE,
    adjust     = adjust_method
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      Treatment_initial = case_when(
        Inoculum == "AM-" & Puccinia == "Puccinia-" ~ "C",
        Inoculum == "AM+" & Puccinia == "Puccinia-" ~ "A",
        Inoculum == "AM-" & Puccinia == "Puccinia+" ~ "P",
        Inoculum == "AM+" & Puccinia == "Puccinia+" ~ "AP",
        TRUE ~ NA_character_
      ),
      Treatment_initial = factor(Treatment_initial, levels = c("C", "A", "P", "AP")),
      sig_letter = gsub(" ", "", .group)
    ) %>%
    dplyr::select(Treatment_initial, sig_letter)
  
  # Letter y positions: just above the upper whisker (per Treatment)
  letters_df <- df %>%
    group_by(Treatment_initial) %>%
    summarise(
      whisker_top = upper_whisker(.data[[response_col]]),
      .groups = "drop"
    ) %>%
    mutate(
      Panel    = response_col,
      y_letter = whisker_top + letter_frac * yrng
    ) %>%
    left_join(cld_tbl, by = "Treatment_initial") %>%
    dplyr::select(Panel, Treatment_initial, y_letter, sig_letter)
  
  list(mod = mod, letters = letters_df, anova = anova_df)
}

# ==============================
# Long-format data for faceting
# ==============================
panel_map <- tibble(
  Panel = c("Root..g.", "Shoot..g.", "Root.shoot.ratio"),
  Panel_label = c("Root mass (g)", "Shoot mass (g)", "Root:shoot ratio")
)

barley_long <- barley %>%
  pivot_longer(
    cols = all_of(panel_map$Panel),
    names_to = "Panel",
    values_to = "Value"
  ) %>%
  left_join(panel_map, by = "Panel") %>%
  mutate(
    Panel = factor(Panel, levels = panel_map$Panel),
    Panel_label = factor(Panel_label, levels = panel_map$Panel_label)
  )

# ==============================
# Compute letters + ANOVA labels per panel
# ==============================
stats_list <- panel_map$Panel %>%
  set_names() %>%
  map(function(resp) get_panel_stats(barley, response_col = resp))

letters_df <- map_dfr(stats_list, "letters")
anova_df   <- map_dfr(stats_list, "anova")

# Join panel labels for facets
letters_df <- letters_df %>%
  left_join(panel_map, by = c("Panel" = "Panel")) %>%
  mutate(Panel_label = factor(Panel_label, levels = panel_map$Panel_label))

anova_df <- anova_df %>%
  left_join(panel_map, by = c("Panel" = "Panel")) %>%
  mutate(Panel_label = factor(Panel_label, levels = panel_map$Panel_label))

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
  axis.text.y      = element_text(size = 12, face = "bold", colour = "black"),
  axis.text.x      = element_text(size = 14, face = "bold", colour = "black"),
  strip.text       = element_text(size = 16, face = "bold"),
  axis.title.y     = element_text(size = 18, face = "bold"),
  legend.text      = element_text(size = 22, face = "bold"),
  axis.line        = element_line(colour = "black", linewidth = 1),
  legend.key       = element_blank(),
  aspect.ratio     = 1
)

# ==============================
# Figure S1 – faceted boxplots with letters + ANOVA text
# ==============================
figS1 <- ggplot(barley_long, aes(x = Treatment_initial, y = Value, fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(Panel_label ~ ., scales = "free_y") +
  ylab(NULL) +
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
# Panel labels (A/B/C)
# ==============================
figS1_labelled <- ggdraw() +
  draw_plot(figS1) +
  draw_plot_label(
    label = c("(A)", "(B)", "(C)"),
    size  = 20,
    x     = c(0.06, 0.37, 0.68),
    y     = c(0.99, 0.99, 0.99)
  )

figS1_labelled

# ==============================
# Save
# ==============================
ggsave(
  filename = file.path(dir_figures, "Figure_S1_biomass_panels.png"),
  plot     = figS1_labelled,
  width    = 30,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

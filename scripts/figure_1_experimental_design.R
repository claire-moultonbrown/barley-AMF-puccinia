# ##############################################################################

# Figure 1 - Experimental design

# ##############################################################################
#script only included for completeness

# ==============================
# Clear environment
# ==============================
rm(list = ls())

# ==============================
# Define directories
# ==============================
dir_figures   <- "figures"

# ==============================
# Load libraries
# ==============================
library(magick)
library(cowplot)

# ==============================
# Fig 1: Experimental design image
# ==============================
exp_design_pic <- file.path(dir_figures, "fig1_experimental_plan_fig.png")

fig1 <- ggdraw() +
  draw_image(exp_design_pic, scale = 1)
fig1

# ==============================
# Save figure (PNG)
# ==============================
ggsave(
  filename = file.path(dir_figures, "Figure_1.png"),
  plot     = fig1,
  width    = 28,
  height   = 15,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

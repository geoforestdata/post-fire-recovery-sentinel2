# Builds the overview figures used at the top of the README.
#
# Usage (from the repository root, after 01_build_dataset.R has written data/):
#   Rscript R/make_hero_figure.R            # reads data/
#   Rscript R/make_hero_figure.R some/dir   # reads another folder
#
# Outputs:
#   figures/hero.png            overview for the README (3 panels)
#   figures/social_preview.png  1280 x 640 image for GitHub's social preview

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
source("R/functions.R")

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"

dat <- load_data(data_dir)
pts <- dat$pts
r <- dat$raster
dir.create("figures", showWarnings = FALSE)

# Same logistic model as in analysis.Rmd
logit_fit <- glm(recovered ~ dNBR + NDVI_pre + elevation + slope,
                 data = pts, family = binomial)
grid <- tibble(
  dNBR = seq(min(pts$dNBR), max(pts$dNBR), length.out = 100),
  NDVI_pre = median(pts$NDVI_pre),
  elevation = median(pts$elevation),
  slope = median(pts$slope)
)
pr <- predict(logit_fit, grid, type = "link", se.fit = TRUE)
grid <- grid |>
  mutate(p = plogis(pr$fit),
         lower = plogis(pr$fit - 1.96 * pr$se.fit),
         upper = plogis(pr$fit + 1.96 * pr$se.fit))

base_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        legend.key.height = grid::unit(0.9, "cm"))

map_theme <- theme(plot.title = element_text(face = "bold", size = 11),
                   legend.key.height = grid::unit(0.9, "cm"))

p1 <- plot_raster_layer(r, "dNBR", "Fire severity (dNBR)", option = "magma") + map_theme
p2 <- plot_raster_layer(r, "NDVI_3", "NDVI, year +3", option = "viridis") + map_theme

p3 <- ggplot(grid, aes(dNBR, p)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#B4443C", alpha = 0.2) +
  geom_line(colour = "#B4443C", linewidth = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Modelled probability of recovery",
       x = "dNBR", y = "Probability") +
  base_theme

titre <- plot_annotation(
  title = "Post-fire vegetation recovery from Sentinel-2",
  subtitle = "Severity, year +3 NDVI and modelled probability of recovery",
  theme = theme(plot.title = element_text(face = "bold", size = 16),
                plot.subtitle = element_text(colour = "grey30", size = 11))
)

ggsave("figures/hero.png", p1 + p2 + p3 + titre,
       width = 14, height = 4.8, dpi = 200, bg = "white")

ggsave("figures/social_preview.png", p1 + p2 + titre,
       width = 12.8, height = 6.4, dpi = 100, bg = "white")

cat("Figures written to figures/hero.png and figures/social_preview.png\n")

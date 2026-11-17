# Helper functions for analysis.Rmd

# Three severity classes from common dNBR thresholds (USGS / FIREMON):
# below 0.27 (low or unburned), 0.27 to 0.66 (moderate), 0.66 and above (high).
# These thresholds are generic: they have not been calibrated for the boreal forest,
# and dNBR is computed here with the year +1 summer (extended assessment).
# Three classes rather than five keep enough pixels in each class.
classify_severity <- function(dnbr) {
  cut(
    dnbr,
    breaks = c(-Inf, 0.27, 0.66, Inf),
    labels = c("Low", "Moderate", "High"),
    right = FALSE,
    ordered_result = FALSE
  )
}

# Reads the dataset written by 01_build_dataset.R and prepares the analysis table.
# Cleaning steps: drop missing values, drop pixels with pre-fire NDVI below
# `min_ndvi_pre` (water and non-vegetated surfaces), add severity class, recovery
# index (year +3 NDVI / pre-fire NDVI) and a binary "recovered" flag.
# Point coordinates are projected to the raster CRS (metres).
load_data <- function(data_dir = "data", min_ndvi_pre = 0.3, recovery_threshold = 0.9) {
  pts_raw <- readr::read_csv(file.path(data_dir, "sample_points.csv"),
                             show_col_types = FALSE)
  r <- terra::rast(file.path(data_dir, "raster_stack.tif"))

  pts_sf <- sf::st_as_sf(pts_raw, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    sf::st_transform(sf::st_crs(terra::crs(r)))
  xy <- sf::st_coordinates(pts_sf)

  pts <- pts_raw |>
    dplyr::mutate(x = xy[, 1], y = xy[, 2]) |>
    tidyr::drop_na() |>
    dplyr::filter(NDVI_pre >= min_ndvi_pre) |>
    dplyr::mutate(
      severity = classify_severity(dNBR),
      RI3 = NDVI_3 / NDVI_pre,
      recovered = as.integer(RI3 >= recovery_threshold)
    )

  list(pts = pts, raster = r, n_raw = nrow(pts_raw))
}

# Map of one layer of a SpatRaster with ggplot2.
# Large rasters are aggregated to keep the rendering light.
plot_raster_layer <- function(r, layer, title, option = "viridis", limits = NULL) {
  if (terra::ncell(r) > 3e5) {
    r <- terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / 3e5)))
  }
  d <- as.data.frame(r[[layer]], xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"

  # Outline of the burned area, taken from the non-missing cells of the raster
  footprint <- terra::ifel(is.na(r[[layer]]), NA, 1)
  outline <- sf::st_as_sf(terra::as.polygons(footprint, dissolve = TRUE))

  # Scale bar (a round length close to 25% of the map width) and north arrow
  xr <- range(d$x)
  yr <- range(d$y)
  w <- diff(xr)
  h <- diff(yr)
  lengths_m <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
  bar <- max(lengths_m[lengths_m <= 0.25 * w], 100)
  bar_label <- if (bar >= 1000) paste(bar / 1000, "km") else paste(bar, "m")
  x0 <- xr[1] + 0.04 * w
  y0 <- yr[1] + 0.05 * h
  xn <- xr[2] - 0.05 * w

  ggplot2::ggplot() +
    ggplot2::geom_raster(data = d, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_sf(data = outline, fill = NA, colour = "grey20", linewidth = 0.3) +
    ggplot2::annotate("segment", x = x0, xend = x0 + bar, y = y0, yend = y0,
                      linewidth = 1, colour = "grey15") +
    ggplot2::annotate("text", x = x0 + bar / 2, y = y0 + 0.035 * h,
                      label = bar_label, size = 3, colour = "grey15") +
    ggplot2::annotate("segment", x = xn, xend = xn, y = yr[2] - 0.17 * h, yend = yr[2] - 0.05 * h,
                      arrow = grid::arrow(length = grid::unit(0.18, "cm")),
                      linewidth = 0.7, colour = "grey15") +
    ggplot2::annotate("text", x = xn, y = yr[2] - 0.02 * h, label = "N",
                      size = 3.2, fontface = "bold", colour = "grey15") +
    ggplot2::scale_fill_viridis_c(option = option, limits = limits, name = layer) +
    ggplot2::coord_sf(crs = sf::st_crs(terra::crs(r)), datum = NA, expand = FALSE) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

# Global Moran's I for a vector of values (for example residuals),
# with a k-nearest-neighbours neighbourhood.
moran_global <- function(coords, values, k = 8) {
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  lw <- spdep::nb2listw(nb, style = "W")
  spdep::moran.test(values, lw)
}

# Extracts I and p from a moran.test object as a list
summarize_moran <- function(m) {
  list(I = unname(m$estimate["Moran I statistic"]), p = m$p.value)
}

# Formats a p-value
format_p <- function(p) {
  if (p < 0.001) "p < 0.001" else paste0("p = ", format(round(p, 3), nsmall = 3))
}

# Area under the ROC curve (AUC) from the Mann-Whitney statistic
auc_roc <- function(prob, observed) {
  n1 <- sum(observed == 1)
  n0 <- sum(observed == 0)
  w <- suppressWarnings(
    stats::wilcox.test(prob[observed == 1], prob[observed == 0], exact = FALSE)$statistic
  )
  unname(w / (n1 * n0))
}

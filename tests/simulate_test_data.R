# SYNTHETIC DATA, FOR TESTING THE CODE ONLY.
# This script generates files that mimic the structure of the 01_build_dataset.R
# export so that analysis.Rmd can be checked from start to finish.
# Results obtained with these data have no ecological meaning and must never be
# published.
#
# Usage (from the repository root):
#   Rscript tests/simulate_test_data.R
#   rmarkdown::render("analysis.Rmd", params = list(data_dir = "tests/synthetic_data"))

library(terra)
library(sf)

set.seed(1)
out_dir <- "tests/synthetic_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 150
r0 <- rast(nrows = n, ncols = n, xmin = 400000, xmax = 400000 + n * 20,
           ymin = 5330000, ymax = 5330000 + n * 20, crs = "EPSG:32617")

field <- function(sd = 1, window = 15) {
  v <- r0
  values(v) <- rnorm(ncell(v), sd = sd)
  focal(v, w = window, fun = "mean", na.rm = TRUE)
}
norm01 <- function(v) (v - global(v, "min", na.rm = TRUE)[[1]]) /
  (global(v, "max", na.rm = TRUE)[[1]] - global(v, "min", na.rm = TRUE)[[1]])

dnbr      <- -0.05 + 1.0 * norm01(field())
elevation <- 300 + 120 * norm01(field(window = 25))
slope     <- 2 + 10 * norm01(field(window = 9))
ndvi_pre  <- 0.65 + 0.15 * norm01(field(window = 21))

noise <- function() 0.03 * field(window = 5) * 10
ndvi_1 <- ndvi_pre * (1 - 0.7 * dnbr) + noise()
ndvi_2 <- ndvi_pre * (1 - 0.5 * dnbr) + noise()
ndvi_3 <- ndvi_pre * (1 - 0.35 * dnbr) + 0.0002 * (elevation - 300) + noise()

stack <- c(dnbr, ndvi_pre * 0 + 0.2, ndvi_pre, ndvi_1, ndvi_2, ndvi_3, elevation, slope)
names(stack) <- c("dNBR", "NBR_pre", "NDVI_pre", "NDVI_1", "NDVI_2", "NDVI_3",
                  "elevation", "slope")
# Irregular "fire perimeter": everything outside is NA, as in the real raster stack
dist_c <- r0
xy_c <- xyFromCell(r0, seq_len(ncell(r0)))
values(dist_c) <- sqrt((xy_c[, 1] - mean(xy_c[, 1]))^2 + (xy_c[, 2] - mean(xy_c[, 2]))^2)
dist_c <- dist_c / global(dist_c, "max")[[1]]
inside <- ((1 - dist_c) + 0.6 * (norm01(field(window = 31)) - 0.5)) > 0.45
stack <- mask(stack, ifel(inside, 1, NA))

writeRaster(stack, file.path(out_dir, "raster_stack.tif"), overwrite = TRUE)

pts <- spatSample(stack, size = 1500, method = "random", as.points = TRUE,
                  na.rm = TRUE)
pts_ll <- project(pts, "EPSG:4326")
d <- as.data.frame(pts)
xy <- crds(pts_ll)
d$lon <- xy[, 1]
d$lat <- xy[, 2]
write.csv(d, file.path(out_dir, "sample_points.csv"), row.names = FALSE)
cat("Synthetic files written to", out_dir, "\n")

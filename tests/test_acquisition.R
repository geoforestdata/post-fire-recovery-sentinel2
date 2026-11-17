# Offline test of the compositing and dataset logic (no network, no rstac needed).
# Builds fake Sentinel-2 scenes as local GeoTIFFs and swaps the two network functions.
# For testing only; outputs go to tests/tmp_acq (git-ignored).
suppressMessages({library(terra); library(dplyr)})
source("R/acquisition.R")
set.seed(1)
out <- "tests/tmp_acq"; unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)

# Fire perimeter (lon/lat polygon near Rouyn-Noranda), saved as gpkg
poly <- vect("POLYGON((-79.10 48.20, -79.00 48.20, -79.00 48.26, -79.10 48.26, -79.10 48.20))", crs = "EPSG:4326")
writeVector(poly, file.path(out, "perim.gpkg"), overwrite = TRUE)

# Scenes live in UTM 17N with 10 m pixels, slightly larger than the perimeter
aoi <- project(poly, "EPSG:32617")
e <- ext(aoi) + 2000
base <- rast(e, resolution = 10, crs = "EPSG:32617")
xy <- crds(base, na.rm = FALSE)
gradient <- (xy[, 1] - min(xy[, 1])) / diff(range(xy[, 1]))

fake_scene <- function(year, k, cloudy = FALSE) {
  # true reflectance: pre-fire healthy, burnt in the fire year, recovering afterwards
  burn <- 0.8 * gradient                       # severity gradient west to east
  rec  <- pmin(1, 0.35 * (year - FIRE) )       # recovery fraction
  dmg  <- if (year < FIRE) 0 else burn * (1 - rec)
  nir  <- 0.35 - 0.20 * dmg + rnorm(length(dmg), 0, 0.005)
  red  <- 0.04 + 0.08 * dmg + rnorm(length(dmg), 0, 0.003)
  sw   <- 0.12 + 0.25 * dmg + rnorm(length(dmg), 0, 0.005)
  scl  <- rep(4, length(dmg))
  if (cloudy) { i <- sample(length(dmg), length(dmg) * 0.3); nir[i] <- 0.9; red[i] <- 0.9; sw[i] <- 0.9; scl[i] <- 9 }
  off <- if (year >= 2022) 1000 else 0          # baseline 04.00 offset
  mk <- function(v, nm) {
    r <- base; values(r) <- round(v * 10000) + off
    f <- file.path(out, sprintf("%d_%d_%s.tif", year, k, nm))
    writeRaster(r, f, datatype = "INT2U", overwrite = TRUE); f
  }
  rs <- base; values(rs) <- scl
  fs <- file.path(out, sprintf("%d_%d_scl.tif", year, k)); writeRaster(rs, fs, datatype = "INT1U", overwrite = TRUE)
  tibble(id = sprintf("fake_%d_%d", year, k), datetime = sprintf("%d-07-%02dT15:00:00Z", year, 10 + k),
         cloud = if (cloudy) 30 else 5, baseline = if (year >= 2022) "05.10" else "03.00",
         red = mk(red, "red"), nir = mk(nir, "nir"), swir22 = mk(sw, "sw"), scl = fs)
}
FIRE <- 2019
scene_db <- lapply(FIRE + (-1:3), function(y) bind_rows(fake_scene(y, 1), fake_scene(y, 2, TRUE), fake_scene(y, 3)))
names(scene_db) <- FIRE + (-1:3)
search_fake <- function(aoi_ll, year, ...) scene_db[[as.character(year)]]
dem_fake <- function(aoi_ll, template) {
  el <- template; v <- crds(el, na.rm = FALSE); values(el) <- 300 + 0.02 * (v[, 2] - min(v[, 2]))
  out <- c(el, terrain(el, "slope", unit = "degrees")); names(out) <- c("elevation", "slope"); out
}

res <- build_dataset(file.path(out, "perim.gpkg"), FIRE, res = 20, n_points = 500,
                     out_dir = file.path(out, "data"), search_fun = search_fake, dem_fun = dem_fake)

pts <- read.csv(file.path(out, "data", "sample_points.csv"))
st <- rast(file.path(out, "data", "raster_stack.tif"))
stopifnot(identical(names(pts), c("dNBR","NBR_pre","NDVI_pre","NDVI_1","NDVI_2","NDVI_3","elevation","slope","lon","lat")),
          nrow(pts) == 500, !anyNA(pts),
          all(pts$lon > -79.11 & pts$lon < -78.99), all(pts$lat > 48.19 & pts$lat < 48.27))
# 1) offset removed and clouds masked: pre-fire NDVI must be ~ (0.35-0.04)/(0.35+0.04) = 0.795
stopifnot(abs(median(pts$NDVI_pre) - 0.795) < 0.02)
# 2) severity gradient is recovered: dNBR increases eastwards
print(cor(pts$dNBR, pts$lon, method = "spearman"))
stopifnot(cor(pts$dNBR, pts$lon, method = "spearman") > 0.9)
# 3) recovery over time where it burned
hi <- pts$dNBR > quantile(pts$dNBR, 0.8)
stopifnot(mean(pts$NDVI_1[hi]) < mean(pts$NDVI_2[hi]), mean(pts$NDVI_2[hi]) < mean(pts$NDVI_3[hi]))
# 4) the repo's loader reads the output
source("R/functions.R")
d <- load_data(file.path(out, "data"))
cat("n points after filter:", nrow(d$pts), " | raster layers:", paste(names(d$raster), collapse = ", "), "\n")
cat("ALL ACQUISITION TESTS PASSED\n")

# Sentinel-2 acquisition and compositing in R (STAC + terra).
# No account, no Python, no JavaScript. Scenes are read directly from public
# cloud-optimized GeoTIFFs, so only the pixels that overlap the fire are downloaded.
#
# Main entry point: build_dataset(). The network access is isolated in
# search_scenes() and read_dem_stac(), which can be replaced (see the `search_fun`
# and `dem_fun` arguments) to run the whole pipeline offline, as tests/ does.

`%||%` <- function(a, b) if (is.null(a)) b else a

STAC_URL <- "https://earth-search.aws.element84.com/v1"

# SCL classes removed from the composites:
# 0 no data, 3 cloud shadow, 8 and 9 clouds, 10 thin cirrus, 11 snow or ice.
SCL_BAD <- c(0, 3, 8, 9, 10, 11)

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

# UTM EPSG code of the centroid of an area of interest in lon/lat
utm_epsg <- function(aoi_ll) {
  xy <- terra::crds(terra::centroids(aoi_ll))
  zone <- floor((xy[1, 1] + 180) / 6) + 1
  epsg <- (if (xy[1, 2] >= 0) 32600 else 32700) + zone
  paste0("EPSG:", epsg)
}

# Empty raster covering the area of interest on a regular grid (metres)
make_template <- function(aoi_utm, res = 20) {
  e <- terra::align(terra::ext(aoi_utm), res)
  terra::rast(e, resolution = res, crs = terra::crs(aoi_utm))
}

# ---------------------------------------------------------------------------
# Scene search (network)
# ---------------------------------------------------------------------------

# Searches Sentinel-2 L2A scenes of one summer and returns a table with one row per
# scene: id, datetime, cloud cover, processing baseline and the asset URLs.
search_scenes <- function(aoi_ll, year, start_md = "06-15", end_md = "09-15",
                          cloud_max = 60, max_scenes = 40, stac_url = STAC_URL) {
  bb <- as.vector(terra::ext(aoi_ll))
  bbox <- c(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])

  items <- rstac::stac(stac_url) |>
    rstac::stac_search(
      collections = "sentinel-2-l2a",
      bbox = bbox,
      datetime = sprintf("%d-%sT00:00:00Z/%d-%sT23:59:59Z", year, start_md, year, end_md),
      limit = 200
    ) |>
    rstac::post_request() |>
    rstac::items_fetch()

  feats <- items$features
  if (length(feats) == 0) stop(sprintf("No Sentinel-2 scene found for %d.", year))

  tab <- dplyr::tibble(
    id       = vapply(feats, function(f) f$id, character(1)),
    datetime = vapply(feats, function(f) f$properties$datetime %||% NA_character_, character(1)),
    cloud    = vapply(feats, function(f) as.numeric(f$properties[["eo:cloud_cover"]] %||% NA), numeric(1)),
    baseline = vapply(feats, function(f) as.character(f$properties[["s2:processing_baseline"]] %||% NA_character_), character(1)),
    red      = vapply(feats, function(f) f$assets$red$href %||% NA_character_, character(1)),
    nir      = vapply(feats, function(f) f$assets$nir$href %||% NA_character_, character(1)),
    swir22   = vapply(feats, function(f) f$assets$swir22$href %||% NA_character_, character(1)),
    scl      = vapply(feats, function(f) f$assets$scl$href %||% NA_character_, character(1))
  )

  tab <- tab[stats::complete.cases(tab[, c("red", "nir", "swir22", "scl")]), ]
  tab <- tab[is.na(tab$cloud) | tab$cloud <= cloud_max, ]
  tab <- tab[order(tab$cloud), ]
  if (nrow(tab) == 0) stop(sprintf("No usable scene for %d (cloud cover <= %d%%).", year, cloud_max))
  utils::head(tab, max_scenes)
}

# ---------------------------------------------------------------------------
# Reading and compositing
# ---------------------------------------------------------------------------

# Digital-number offset of Sentinel-2 L2A. Scenes processed with baseline 04.00 or
# later (January 2022 onward) carry a +1000 offset that must be removed so that
# reflectance is comparable across years (this is the "harmonization" GEE applies).
offset_dn <- function(baseline, datetime) {
  if (is.na(baseline) || !nzchar(baseline)) {
    return(if (!is.na(datetime) && substr(datetime, 1, 10) >= "2022-01-25") 1000 else 0)
  }
  if (baseline >= "04.00") 1000 else 0
}

# Reads one asset (URL or local path), crops it to the area of interest in the
# scene's own CRS and resamples it to the template grid.
read_band <- function(href, template, method = "bilinear", zero_is_na = TRUE) {
  path <- if (grepl("^https?://", href)) paste0("/vsicurl/", href) else href
  src <- terra::rast(path)
  window <- terra::project(
    terra::as.polygons(terra::ext(template), crs = terra::crs(template)),
    terra::crs(src)
  )
  src <- terra::crop(src, terra::ext(window), snap = "out")
  if (zero_is_na) src <- terra::ifel(src == 0, NA, src)
  terra::project(src, template, method = method)
}

# Median composite (red, nir, swir22 reflectance) of a set of scenes
build_composite <- function(scenes, template, verbose = TRUE) {
  bands <- c("red", "nir", "swir22")
  acc <- stats::setNames(vector("list", length(bands)), bands)
  tmp <- function() tempfile(fileext = ".tif")

  for (i in seq_len(nrow(scenes))) {
    s <- scenes[i, ]
    if (verbose) message(sprintf("  scene %d/%d (%s)", i, nrow(scenes), s$id))
    layers <- tryCatch({
      scl <- read_band(s$scl, template, method = "near", zero_is_na = FALSE)
      bad <- scl %in% SCL_BAD
      off <- offset_dn(s$baseline, s$datetime)
      lapply(bands, function(b) {
        x <- read_band(s[[b]], template, method = "bilinear")
        x <- terra::clamp(x - off, lower = 0, values = TRUE) / 10000
        terra::mask(x, bad, maskvalues = TRUE, filename = tmp(),
                    overwrite = TRUE, datatype = "FLT4S")
      })
    }, error = function(e) {
      message("    skipped: ", conditionMessage(e))
      NULL
    })
    if (is.null(layers)) next
    for (k in seq_along(bands)) {
      acc[[k]][[length(acc[[k]]) + 1]] <- layers[[k]]
    }
  }

  if (length(acc$red) == 0) stop("No usable scene for this composite.")
  out <- lapply(acc, function(l) terra::app(terra::rast(l), "median", na.rm = TRUE))
  out <- terra::rast(out)
  names(out) <- bands
  out
}

# Composite of one summer, cached on disk so that it is only built once
get_composite <- function(year, aoi_ll, template, cache_dir, search_fun, ...) {
  f <- file.path(cache_dir, sprintf("composite_%d.tif", year))
  if (file.exists(f)) {
    message("Using cached composite: ", f)
    return(terra::rast(f))
  }
  message("Building the ", year, " composite")
  scenes <- search_fun(aoi_ll, year, ...)
  comp <- build_composite(scenes, template)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  terra::writeRaster(comp, f, overwrite = TRUE)
  comp
}

# ---------------------------------------------------------------------------
# Terrain (network)
# ---------------------------------------------------------------------------

# Elevation and slope from the Copernicus 30 m DEM, on the template grid
read_dem_stac <- function(aoi_ll, template, stac_url = STAC_URL) {
  bb <- as.vector(terra::ext(aoi_ll))
  bbox <- c(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])

  items <- rstac::stac(stac_url) |>
    rstac::stac_search(collections = "cop-dem-glo-30", bbox = bbox, limit = 50) |>
    rstac::post_request() |>
    rstac::items_fetch()

  hrefs <- vapply(items$features, function(f) f$assets$data$href %||% NA_character_, character(1))
  hrefs <- hrefs[!is.na(hrefs)]
  if (length(hrefs) == 0) stop("No DEM tile found for the area of interest.")

  tiles <- lapply(hrefs, function(h) {
    tryCatch(read_band(h, template, method = "bilinear", zero_is_na = FALSE),
             error = function(e) NULL)
  })
  tiles <- Filter(Negate(is.null), tiles)
  if (length(tiles) == 0) stop("DEM tiles could not be read.")

  elev <- if (length(tiles) == 1) tiles[[1]] else terra::merge(terra::sprc(tiles))
  out <- c(elev, terra::terrain(elev, v = "slope", unit = "degrees"))
  names(out) <- c("elevation", "slope")
  out
}

# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

# Builds the two files used by analysis.Rmd:
#   <out_dir>/raster_stack.tif   multiband raster (dNBR, NBR_pre, NDVI_pre, NDVI_1..3,
#                                elevation, slope), masked to the fire perimeter
#   <out_dir>/sample_points.csv  random sample of pixels with the same variables + lon/lat
build_dataset <- function(perimeter_path, fire_year,
                          start_md = "06-15", end_md = "09-15",
                          cloud_max = 60, max_scenes = 40,
                          res = 20, n_points = 3000, seed = 42,
                          out_dir = "data",
                          search_fun = search_scenes, dem_fun = read_dem_stac) {
  aoi_ll <- terra::aggregate(terra::project(terra::vect(perimeter_path), "EPSG:4326"))
  aoi <- terra::project(aoi_ll, utm_epsg(aoi_ll))
  template <- make_template(aoi, res)
  cache_dir <- file.path(out_dir, "composites")

  composite <- function(year) {
    get_composite(year, aoi_ll, template, cache_dir, search_fun,
                  start_md = start_md, end_md = end_md,
                  cloud_max = cloud_max, max_scenes = max_scenes)
  }
  pre <- composite(fire_year - 1)
  y1 <- composite(fire_year + 1)
  y2 <- composite(fire_year + 2)
  y3 <- composite(fire_year + 3)

  nbr  <- function(x) (x[["nir"]] - x[["swir22"]]) / (x[["nir"]] + x[["swir22"]])
  ndvi <- function(x) (x[["nir"]] - x[["red"]]) / (x[["nir"]] + x[["red"]])

  message("Reading terrain")
  dem <- dem_fun(aoi_ll, template)

  stack <- c(nbr(pre) - nbr(y1), nbr(pre), ndvi(pre), ndvi(y1), ndvi(y2), ndvi(y3),
             dem[["elevation"]], dem[["slope"]])
  names(stack) <- c("dNBR", "NBR_pre", "NDVI_pre", "NDVI_1", "NDVI_2", "NDVI_3",
                    "elevation", "slope")
  stack <- terra::mask(stack, aoi)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  terra::writeRaster(stack, file.path(out_dir, "raster_stack.tif"), overwrite = TRUE)

  set.seed(seed)
  pts <- terra::spatSample(stack, size = n_points, method = "random",
                           na.rm = TRUE, as.points = TRUE)
  tab <- as.data.frame(pts)
  xy <- terra::crds(terra::project(pts, "EPSG:4326"))
  tab$lon <- xy[, 1]
  tab$lat <- xy[, 2]
  utils::write.csv(tab, file.path(out_dir, "sample_points.csv"), row.names = FALSE)

  message("Done: ", nrow(tab), " sampled pixels written to ", out_dir)
  invisible(list(stack = stack, points = tab))
}

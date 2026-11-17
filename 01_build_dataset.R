# Step 1: builds the dataset from Sentinel-2 scenes, entirely in R.
#
# What it does:
#   1. reads the fire perimeter and sets up a 20 m grid in the right UTM zone;
#   2. for the summer before the fire and for years +1, +2 and +3, searches the public
#      Sentinel-2 L2A catalogue (STAC), masks clouds and shadows with the SCL band,
#      and builds a median composite;
#   3. computes dNBR, NBR and NDVI, adds elevation and slope (Copernicus 30 m DEM);
#   4. writes data/raster_stack.tif and data/sample_points.csv.
#
# Requirements: internet access and the packages listed in README.md. No account.
# The first run reads several dozen scenes per summer and can take a while.
# Each yearly composite is cached in data/composites/, so a second run is fast.

source("R/acquisition.R")

# ---- Configuration ---------------------------------------------------------

perimeter_path <- "data/fire_perimeter.gpkg"  # your fire perimeter (.gpkg, .shp, ...)
fire_year      <- 2023                        # the summer of year +3 must be over

start_md   <- "06-15"   # summer window used for the composites (month-day)
end_md     <- "09-15"
cloud_max  <- 60        # maximum scene cloud cover (%)
max_scenes <- 40        # scenes kept per summer (the least cloudy ones)
res        <- 20        # grid resolution (m)
n_points   <- 3000      # number of sampled pixels
seed       <- 42

# ---- Run -------------------------------------------------------------------

build_dataset(
  perimeter_path = perimeter_path,
  fire_year = fire_year,
  start_md = start_md, end_md = end_md,
  cloud_max = cloud_max, max_scenes = max_scenes,
  res = res, n_points = n_points, seed = seed,
  out_dir = "data"
)

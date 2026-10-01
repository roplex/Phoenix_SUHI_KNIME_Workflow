rm(list = ls())
library(terra)
library(tigris)
options(tigris_use_cache = TRUE)

# 1. Define Paths
worldcover_raw_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/WorldCover/ESA_WorldCover_10m_2021_v200_N33W114_Map.tif"
lst_wide_path        <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/LST/LST_composite_superrelaxed_tsharp.tif" # Update to your wide LST file path

output_dir          <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Preprocessed"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Fetch & Save Phoenix Urbanized Area Boundary (.gpkg)
cat("Fetching Phoenix Urban Area boundary...\n")
ua <- urban_areas(year = 2020, progress_bar = TRUE)

# Filter using the 2010 schema name column (NAME10)
phx_ua <- ua[grepl("Phoenix", ua$NAME10, ignore.case = TRUE) &
               grepl(", AZ$", ua$NAME10), ]

if (nrow(phx_ua) == 0) {
  stop("Phoenix urban area not found in ua$NAME10.")
}

# Convert sf to SpatVector
phx_ua_vect <- vect(phx_ua)

cat("Selected UA:", phx_ua_vect$NAME10, "| area (km^2):", expanse(phx_ua_vect, unit = "km"), "\n")

# Keep the primary Phoenix urban footprint if multiple features match
if (length(phx_ua_vect) > 1) {
  phx_ua_vect <- phx_ua_vect[which.max(expanse(phx_ua_vect)), ]
}

writeVector(phx_ua_vect, file.path(output_dir, "phoenix_urban_area.gpkg"), overwrite = TRUE)
cat("Successfully saved Phoenix Urban Area boundary to phoenix_urban_area.gpkg!\n")

# 3. Resample & Crop ESA WorldCover to Match Wide LST Surface (.tif)
cat("Cropping and resampling WorldCover raster...\n")
lst_wide   <- rast(lst_wide_path)[[1]]
worldcover <- rast(worldcover_raw_path)

# Reproject Phoenix UA to LST CRS for spatial alignment
phx_ua_proj <- project(phx_ua_vect, crs(lst_wide))

# Convert LST extent to a spatial polygon and reproject to WorldCover's CRS
lst_poly_wc_crs <- project(as.polygons(ext(lst_wide), crs = crs(lst_wide)), crs(worldcover))

# Crop WorldCover to match wide LST extent, then project & align to LST raster grid
wc_cropped <- crop(worldcover, lst_poly_wc_crs)
wc_aligned <- project(wc_cropped, lst_wide, method = "near")

# Classify non-built-up classes (20: Shrub, 30: Grass, 40: Crop, 60: Bare)
non_urban_classes <- c(20, 30, 40, 60)
wc_non_urban <- classify(
  wc_aligned,
  cbind(non_urban_classes, rep(1, length(non_urban_classes))),
  others = NA
)

# Save lightweight non-urban mask raster
writeRaster(wc_non_urban, file.path(output_dir, "worldcover_non_urban_aligned.tif"), overwrite = TRUE)
cat("Pre-processing complete! Files written to:", output_dir, "\n")
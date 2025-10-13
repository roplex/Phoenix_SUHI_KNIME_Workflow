library(terra)

# --- Inputs from KNIME ---
lst_path  <- as.character(knime.in[["lst_clipped_path"]])   # clipped/resampled LST
ndvi_path <- as.character(knime.in[["ndvi_clipped_path"]]) # clipped NDVI
city_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData/Boundary/City_Limit_Light_Outline.geojson"

stopifnot(file.exists(lst_path), file.exists(ndvi_path), file.exists(city_path))

# --- Load rasters and polygon ---
lst  <- rast(lst_path)[[1]]
ndvi <- rast(ndvi_path)[[1]]
city <- vect(city_path)

# --- Rasterize city boundary to match LST grid ---
urban_rast <- rasterize(city, lst, field = 1, background = 0)

# --- Save output mask ---
out_mask <- file.path(dirname(lst_path), "urban_mask_Phoenix.tif")
writeRaster(urban_rast, out_mask, overwrite=TRUE)

# --- Return paths downstream ---
knime.out <- data.frame(
  urban_mask_path = out_mask,
  lst_path  = lst_path,
  ndvi_path = ndvi_path
)

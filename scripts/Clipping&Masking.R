library(terra)

# ---- Inputs from KNIME ----
lst_resampled_path <- as.character(knime.in[["resampled_lst_path"]])  # LST composite resampled to NDVI grid
ndvi_aligned_path  <- as.character(knime.in[["ndvi_aligned_path"]])   # aligned NDVI composite
phoenix_path       <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData/Boundary/City_Limit_Light_Outline.geojson"

# ---- Load rasters and AOI ----
lst_resampled <- rast(lst_resampled_path)
ndvi_aligned  <- rast(ndvi_aligned_path)
phoenix       <- vect(phoenix_path)

# ---- Clip + mask LST ----
lst_clipped <- mask(crop(lst_resampled, phoenix), phoenix)

# ---- Also clip NDVI for consistency ----
ndvi_clipped <- mask(crop(ndvi_aligned, phoenix), phoenix)

# ---- Save outputs ----
out_lst  <- file.path(dirname(lst_resampled_path), "LST_resampled_clipped_Phoenix.tif")
out_ndvi <- file.path(dirname(ndvi_aligned_path),  "NDVI_resampled_clipped_Phoenix.tif")

writeRaster(lst_clipped,  out_lst,  overwrite=TRUE)
writeRaster(ndvi_clipped, out_ndvi, overwrite=TRUE)

# ---- Return outputs to KNIME ----
knime.out <- data.frame(
  lst_clipped_path  = out_lst,
  ndvi_clipped_path = out_ndvi
)

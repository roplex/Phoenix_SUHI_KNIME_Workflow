library(terra)

# ---- Inputs from KNIME ----
lst_path  <- as.character(knime.in[["lst_composite_path"]])   # e.g., superrelaxed LST composite
ndvi_path <- as.character(knime.in[["ndvi_composite_path"]])  # strict NDVI composite (reference grid)

# ---- Load rasters ----
lst  <- rast(lst_path)
ndvi <- rast(ndvi_path)

# ---- Resample LST to NDVI resolution ----
lst_resampled <- resample(lst, ndvi, method="bilinear")

# ---- Save resampled LST ----
out_lst  <- gsub(".tif", "_resampledToNDVI.tif", lst_path)

# ---- Save NDVI aligned copy in same folder ----
out_ndvi <- file.path(dirname(out_lst), "NDVI_composite_aligned.tif")
writeRaster(ndvi, out_ndvi, overwrite=TRUE)

# ---- Save resampled LST ----
writeRaster(lst_resampled, out_lst, overwrite=TRUE)

# ---- Return results to KNIME ----
knime.out <- data.frame(
  resampled_lst_path = out_lst,
  ndvi_aligned_path  = out_ndvi
)

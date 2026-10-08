# ============================================================================
# KNIME R Snippet -- Node: "Clipping & Masking"
# Unchanged logic for the urban-side clip; adds pass-through of the wide,
# pre-clip raster so the Urban vs Rural Mask node has genuine countryside
# to work with (the old rural mask collapsed because it rasterized the
# Phoenix polygon onto a raster already clipped to that same polygon).
#
# Inputs (unchanged):
#   knime.in[["resampled_lst_path"]]  -- TsHARP-sharpened LST, wide extent
#   knime.in[["ndvi_aligned_path"]]   -- aligned NDVI, wide extent
#
# Outputs:
#   lst_clipped_path, ndvi_clipped_path  -- unchanged, city-clipped rasters
#   lst_wide_path, ndvi_wide_path        -- NEW: the wide rasters, untouched,
#                                            for the redefined rural mask
# ============================================================================

library(terra)

phx_root <- Sys.getenv("PHX_ROOT", unset = "PhoenixData2")  # set PHX_ROOT to your local data folder (see README)
# ---- Inputs from KNIME ----
lst_resampled_path <- as.character(knime.in[["resampled_lst_path"]])  # wide, TsHARP-sharpened
ndvi_aligned_path  <- as.character(knime.in[["ndvi_aligned_path"]])   # wide, aligned NDVI
phoenix_path       <- file.path(phx_root, "Boundary", "City_Limit_Light_Outline.geojson")

# ---- Load rasters and AOI ----
lst_resampled <- rast(lst_resampled_path)
ndvi_aligned  <- rast(ndvi_aligned_path)
phoenix       <- vect(phoenix_path)

# ---- Clip + mask LST (urban side, unchanged) ----
lst_clipped <- mask(crop(lst_resampled, phoenix), phoenix)

# ---- Also clip NDVI for consistency ----
ndvi_clipped <- mask(crop(ndvi_aligned, phoenix), phoenix)

# ---- Save clipped (urban) outputs ----
out_lst  <- file.path(dirname(lst_resampled_path), "LST_resampled_clipped_Phoenix.tif")
out_ndvi <- file.path(dirname(ndvi_aligned_path),  "NDVI_resampled_clipped_Phoenix.tif")

writeRaster(lst_clipped,  out_lst,  overwrite = TRUE)
writeRaster(ndvi_clipped, out_ndvi, overwrite = TRUE)

# ---- Return outputs to KNIME ----
# lst_wide_path / ndvi_wide_path are simply the input paths passed through
# unmodified -- the Urban vs Rural Mask node needs the full, un-clipped
# extent (approx. 12,000 sq km, Phoenix roughly centered) to find genuine
# non-urban reference pixels; the city-clipped rasters no longer contain any.
knime.out <- data.frame(
  lst_clipped_path  = out_lst,
  ndvi_clipped_path = out_ndvi,
  lst_wide_path      = lst_resampled_path,
  ndvi_wide_path     = ndvi_aligned_path
)

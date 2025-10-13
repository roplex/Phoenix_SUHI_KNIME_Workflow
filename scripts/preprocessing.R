library(terra)

# ---- Get file path from KNIME ----
lst_path <- knime.flow.in[["file_path"]]

# Construct QC path
qc_path <- gsub("LST_Day_1km", "QC_Day", lst_path)
stopifnot(file.exists(lst_path), file.exists(qc_path))

# ---- Load rasters ----
lst <- rast(lst_path)  # Already scaled to Kelvin (GDAL applied scale_factor=0.02)
qc  <- rast(qc_path)

# ---- Define QC rules ----
mask_strict       <- qc %in% c(0, 65)        # strictest
mask_relaxed      <- qc %in% c(0, 17, 65)    # moderate
mask_superrelaxed <- qc != 81                # everything except clouds

# ---- Apply masks in Kelvin ----
lst_strict_K       <- mask(lst, mask_strict, maskvalues=FALSE)
lst_relaxed_K      <- mask(lst, mask_relaxed, maskvalues=FALSE)
lst_superrelaxed_K <- mask(lst, mask_superrelaxed, maskvalues=FALSE)

# ---- Convert to Celsius (Kelvin → °C) ----
lst_strict_C       <- lst_strict_K  - 273.15
lst_relaxed_C      <- lst_relaxed_K - 273.15
lst_superrelaxed_C <- lst_superrelaxed_K - 273.15

# ---- Sanity mask (remove impossible values) ----
# Keep only temperatures within realistic bounds
lst_strict_C[lst_strict_C < -100 | lst_strict_C > 70] <- NA
lst_relaxed_C[lst_relaxed_C < -100 | lst_relaxed_C > 70] <- NA
lst_superrelaxed_C[lst_superrelaxed_C < -100 | lst_superrelaxed_C > 70] <- NA

# ---- Coverage statistics ----
valid_strict   <- global(!is.na(lst_strict_C), "sum", na.rm=TRUE)[[1]]
valid_relaxed  <- global(!is.na(lst_relaxed_C), "sum", na.rm=TRUE)[[1]]
valid_super    <- global(!is.na(lst_superrelaxed_C), "sum", na.rm=TRUE)[[1]]
total_pixels   <- ncell(lst)

coverage_strict   <- round(100 * valid_strict  / total_pixels, 2)
coverage_relaxed  <- round(100 * valid_relaxed / total_pixels, 2)
coverage_super    <- round(100 * valid_super   / total_pixels, 2)

# ---- Save outputs ----
out_strict       <- gsub(".tif", "_clean_strict.tif",       lst_path)
out_relaxed      <- gsub(".tif", "_clean_relaxed.tif",      lst_path)
out_superrelaxed <- gsub(".tif", "_clean_superrelaxed.tif", lst_path)

writeRaster(lst_strict_C,       out_strict,       overwrite=TRUE)
writeRaster(lst_relaxed_C,      out_relaxed,      overwrite=TRUE)
writeRaster(lst_superrelaxed_C, out_superrelaxed, overwrite=TRUE)

# ---- Return results to KNIME ----
knime.out <- data.frame(
  strict_lst_path       = out_strict,
  relaxed_lst_path      = out_relaxed,
  superrelaxed_lst_path = out_superrelaxed,
  coverage_strict       = coverage_strict,
  coverage_relaxed      = coverage_relaxed,
  coverage_superrelaxed = coverage_super
)

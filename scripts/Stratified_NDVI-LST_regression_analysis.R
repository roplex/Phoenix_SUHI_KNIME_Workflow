# ============================================================================
# KNIME R Snippet -- Node: "Stratified NDVI-LST Regression"
# Chapter 7 (SUHI) -- corrected version
#
# Fits one NDVI~LST regression per Phoenix urban village (as opposed to the
# Global node, which fits one regression across the whole city).
#
# Inputs:
#   knime.in[["lst_path"]]   -- city-clipped LST (TsHARP-sharpened)
#   knime.in[["ndvi_path"]]  -- city-clipped NDVI, passed through
#   knime.in[["suhi_path"]]  -- SUHI raster path; used here to flag which
#                                villages contain any of the previously
#                                identified extreme-SUHI pixels
#
# Outputs: one row per village --
#   village, slope, intercept, r_squared, pearson_r, n_points,
#   n_extreme_pixels, r_squared_excl_extreme, suhi_path, lst_path, ndvi_path
# ============================================================================

library(terra)
library(dplyr)

# ---- Inputs from KNIME ----
lst_path  <- as.character(knime.in[["lst_path"]])
ndvi_path <- as.character(knime.in[["ndvi_path"]])
suhi_path <- as.character(knime.in[["suhi_path"]])

stopifnot(file.exists(lst_path), file.exists(ndvi_path), file.exists(suhi_path))

# ---- Villages path ----
# CORRECTED from .../PhoenixData/Boundary/Villages.geojson to
# .../PhoenixData2/Boundary/Villages.geojson -- every other path in this
# pipeline (DEM, WorldCover, city-limit boundary, LST/NDVI composites) lives
# under PhoenixData2; the original PhoenixData (no "2") path is almost
# certainly stale from before the data reorganization.
# VERIFY this file actually exists at this path before trusting results --
# it was not independently confirmed as part of this correction pass.
villages_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Boundary/Villages.geojson"
stopifnot(file.exists(villages_path))

# ---- Load rasters and vector ----
lst      <- rast(lst_path)[[1]]
ndvi     <- rast(ndvi_path)[[1]]
suhi     <- rast(suhi_path)[[1]]
villages <- vect(villages_path)

# ---- Defensive check: confirm the village-name field is actually called
# NAME before relying on it. A tigris field-naming mismatch (NAME10 vs
# NAME20) bit this pipeline earlier, so this kind of check is now standard
# practice rather than an assumption. ----
if (!"NAME" %in% names(villages)) {
  stop(sprintf("Expected a 'NAME' field in villages_path but found: %s -- update village_names below to the correct field.",
               paste(names(villages), collapse = ", ")))
}

# ---- Align CRS between rasters and villages ----
if (!identical(crs(lst), crs(ndvi)))     ndvi     <- project(ndvi, crs(lst))
if (!identical(crs(lst), crs(suhi)))     suhi     <- project(suhi, crs(lst))
if (!identical(crs(lst), crs(villages))) villages <- project(villages, crs(lst))

# ---- Align extent/resolution ----
# compareGeom() replaces all.equal(ext(...)) -- see the Global Regression
# node's header for why the original all.equal()-based check is unreliable
# (it can silently fail to trigger a needed resample).
if (!compareGeom(lst, ndvi, stopOnError = FALSE)) ndvi <- resample(ndvi, lst, method = "bilinear")
if (!compareGeom(lst, suhi, stopOnError = FALSE)) suhi <- resample(suhi, lst, method = "bilinear")

# ---- Stack rasters for extraction ----
stacked <- c(lst, ndvi, suhi)
names(stacked) <- c("LST", "NDVI", "SUHI")

# ---- Extract data and run regressions by village ----
village_names <- villages$NAME
results_list <- list()

for (i in seq_along(village_names)) {
  v <- villages[i]
  df_v <- terra::extract(stacked, v, na.rm = TRUE)

  na_row <- function(n_pts) {
    data.frame(
      village = village_names[i],
      slope = NA_real_, intercept = NA_real_, r_squared = NA_real_, pearson_r = NA_real_,
      n_points = n_pts, n_extreme_pixels = NA_integer_, r_squared_excl_extreme = NA_real_,
      suhi_path = suhi_path, lst_path = lst_path, ndvi_path = ndvi_path
    )
  }

  if (is.null(df_v) || nrow(df_v) == 0) {
    message(paste("Too few valid pixels for:", village_names[i], "-> Filling with NA"))
    results_list[[length(results_list) + 1]] <- na_row(0)
    next
  }

  df_v <- df_v[complete.cases(df_v[, c("LST", "NDVI", "SUHI")]), ]
  if (nrow(df_v) < 5) {
    # FIX: the original script's second length check used `next` here with
    # no row appended, which silently dropped the village from the output
    # entirely -- a village could vanish from results with nothing to flag
    # it. Every village now gets a row (NA-filled if excluded) so nothing
    # goes missing without a trace.
    message(paste("Too few valid pixels for:", village_names[i], "-> Filling with NA"))
    results_list[[length(results_list) + 1]] <- na_row(nrow(df_v))
    next
  }

  fit <- lm(LST ~ NDVI, data = df_v)
  r2 <- summary(fit)$r.squared
  pearson_r <- suppressWarnings(cor(df_v$NDVI, df_v$LST))

  # ---- Per-village extreme-pixel check: smaller per-village sample sizes
  # make an individual leverage point more influential than in the Global
  # regression, so this is worth checking here even where it wasn't very
  # consequential city-wide. Reuses the same SUHI < -8 / > 4 degC threshold
  # as the SUHI Computation node and the Global Regression node. ----
  extreme_idx_v <- which(df_v$SUHI < -8 | df_v$SUHI > 4)
  n_extreme_v <- length(extreme_idx_v)
  if (n_extreme_v > 0 && (nrow(df_v) - n_extreme_v) >= 5) {
    fit_clean_v <- lm(LST ~ NDVI, data = df_v[-extreme_idx_v, ])
    r2_excl_v <- summary(fit_clean_v)$r.squared
  } else {
    r2_excl_v <- r2
  }

  res <- data.frame(
    village = village_names[i],
    slope = coef(fit)[2], intercept = coef(fit)[1], r_squared = r2, pearson_r = pearson_r,
    n_points = nrow(df_v), n_extreme_pixels = n_extreme_v, r_squared_excl_extreme = r2_excl_v,
    suhi_path = suhi_path, lst_path = lst_path, ndvi_path = ndvi_path
  )
  results_list[[length(results_list) + 1]] <- res
}

# ---- Combine all results safely for KNIME output ----
knime.out <- as.data.frame(bind_rows(results_list))

# ---- Flag villages where the extreme-pixel check meaningfully shifts R2 ----
flagged <- knime.out[!is.na(knime.out$n_extreme_pixels) & knime.out$n_extreme_pixels > 0, ]
if (nrow(flagged) > 0) {
  cat("Villages containing previously-flagged extreme-SUHI pixels:\n")
  print(flagged[, c("village", "n_points", "n_extreme_pixels", "r_squared", "r_squared_excl_extreme")])
} else {
  cat("No villages contain any of the previously-flagged extreme-SUHI pixels.\n")
}

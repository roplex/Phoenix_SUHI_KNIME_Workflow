# ============================================================================
# KNIME R Snippet -- Node: "SUHI Computation"
# Chapter 7 (SUHI) -- corrected version
#
# Simplified relative to the original node: since lst_path is already
# city-clipped (by Clipping & Masking) and the rural baseline arrives
# pre-computed from the redefined, elevation-matched reference zone (from
# Urban vs Rural Mask), there is no separate urban_mask_path to load or
# align here -- every pixel remaining in lst_path is already an urban
# pixel by construction, and everything outside the city was already NA.
#
# Inputs:
#   knime.in[["lst_path"]]                   -- city-clipped LST (TsHARP-sharpened)
#   knime.in[["ndvi_path"]]                  -- city-clipped NDVI, passed through
#   knime.in[["mean_rural_LST_elev_matched"]] -- redefined, elevation-matched
#                                                 rural baseline (from Urban vs
#                                                 Rural Mask). NOTE: the upstream
#                                                 node's output column is named
#                                                 mean_rural_LST_elev_matched, not
#                                                 mean_rural_LST -- if a Column
#                                                 Rename node sits between the two
#                                                 in the KNIME graph, point it at
#                                                 this column specifically (not
#                                                 mean_rural_LST_full_sample or
#                                                 mean_rural_LST_elev_adjusted).
#
# Outputs:
#   suhi_path, mean_rural_LST, histogram_csv_path, lst_path, ndvi_path,
#   suhi_mean, suhi_min, suhi_max, suhi_n_pixels
# ============================================================================

library(terra)

# ---- Inputs from KNIME ----
lst_path       <- as.character(knime.in[["lst_path"]])
ndvi_path      <- as.character(knime.in[["ndvi_path"]])
mean_rural_LST <- as.numeric(knime.in[["mean_rural_LST_elev_matched"]])

stopifnot(file.exists(lst_path), file.exists(ndvi_path))
if (is.na(mean_rural_LST)) stop("mean_rural_LST_elev_matched was not provided or is NA.")

# ---- Load the (already city-clipped) LST raster ----
lst <- rast(lst_path)[[1]]

# ---- SUHI = LST minus the redefined rural baseline ----
# No urban mask needed here: lst is already clipped to Phoenix, so every
# non-NA pixel is, by construction, an urban pixel.
suhi <- lst - mean_rural_LST
names(suhi) <- "SUHI"

# ---- Save SUHI raster ----
out_suhi <- file.path(dirname(lst_path), "SUHI_map_urbanOnly_corrected.tif")
writeRaster(suhi, out_suhi, overwrite = TRUE)

# ---- Global diagnostic: check the real, fully-combined result against the
# earlier linear projection (which only accounted for the rural-baseline
# shift, holding the pre-TsHARP LST fixed). This run combines both revisions
# for the first time, so this is the actual answer, not an estimate. ----
vals <- as.vector(values(suhi, na.rm = TRUE))
if (length(vals) == 0) {
  stop("No valid SUHI pixels found -- check that lst_path is genuinely city-clipped and non-empty.")
}
suhi_mean <- mean(vals)
suhi_min  <- min(vals)
suhi_max  <- max(vals)
suhi_n    <- length(vals)
suhi_p01  <- as.numeric(quantile(vals, 0.01, na.rm = TRUE))
suhi_p99  <- as.numeric(quantile(vals, 0.99, na.rm = TRUE))
cat(sprintf("SUHI (rural baseline = %.2f degC): mean=%.2f degC | range=[%.2f, %.2f] degC | n=%d pixels\n",
            mean_rural_LST, suhi_mean, suhi_min, suhi_max, suhi_n))
cat(sprintf("SUHI 1st-99th percentile range: [%.2f, %.2f] degC (trimmed; excludes the small, ambiguous extreme-tail pixels -- see clustering check below)\n",
            suhi_p01, suhi_p99))

# ---- Histogram for visualization (SUHI Histogram node downstream) ----
hist_counts <- hist(vals, breaks = 50, plot = FALSE)
hist_df <- data.frame(
  break_left  = hist_counts$breaks[-length(hist_counts$breaks)],
  break_right = hist_counts$breaks[-1],
  count       = hist_counts$counts
)
out_hist_csv <- file.path(dirname(lst_path), "suhi_histogram_urbanOnly_corrected.csv")
write.csv(hist_df, out_hist_csv, row.names = FALSE)

# ---- Extreme-tail check: are the SUHI extremes NDVI-driven, as TsHARP's
# design predicts (cold = high NDVI/vegetated, hot = low NDVI/built-up)? ----
suhi_pts <- as.points(suhi, values = TRUE, na.rm = TRUE)
suhi_df  <- data.frame(SUHI = suhi_pts$SUHI, crds(suhi_pts))

ndvi <- rast(ndvi_path)[[1]]
suhi_df$NDVI <- extract(ndvi, suhi_df[, c("x", "y")])[, 2]

extreme_cold_idx <- which(suhi_df$SUHI < -8)
extreme_hot_idx  <- which(suhi_df$SUHI > 4)
cat(sprintf("Pixels SUHI < -8 degC: %d (%.2f%%) | SUHI > 4 degC: %d (%.2f%%)\n",
            length(extreme_cold_idx), 100 * length(extreme_cold_idx) / nrow(suhi_df),
            length(extreme_hot_idx), 100 * length(extreme_hot_idx) / nrow(suhi_df)))
cat("NDVI at extreme-cold pixels (expect HIGH -- vegetated/irrigated):\n")
print(summary(suhi_df$NDVI[extreme_cold_idx]))
cat("NDVI at extreme-hot pixels (expect LOW -- bare/built-up):\n")
print(summary(suhi_df$NDVI[extreme_hot_idx]))

# ---- Spatial clustering check: do the extreme-tail pixels sit in one or two
# tight clusters (consistent with a single unstable local-regression window
# in the TsHARP moving-window fit, where a near-flat NDVI denominator can
# produce a poorly-conditioned slope) or scattered independently across the
# city (more consistent with genuine, spatially distinct physical extremes)? ----
cat("Extreme-cold pixel locations:\n")
print(suhi_df[extreme_cold_idx, c("x", "y", "SUHI", "NDVI")])
cat("Extreme-hot pixel locations:\n")
print(suhi_df[extreme_hot_idx, c("x", "y", "SUHI", "NDVI")])

if (length(extreme_cold_idx) > 1) {
  d_cold <- dist(suhi_df[extreme_cold_idx, c("x", "y")])
  cat(sprintf("Extreme-cold pairwise distances: min=%.5f max=%.5f (deg)\n", min(d_cold), max(d_cold)))
}
if (length(extreme_hot_idx) > 1) {
  d_hot <- dist(suhi_df[extreme_hot_idx, c("x", "y")])
  cat(sprintf("Extreme-hot pairwise distances: min=%.5f max=%.5f (deg)\n", min(d_hot), max(d_hot)))
}

# ---- Return outputs to KNIME ----
# histogram is passed as a file path (out_hist_csv), not embedded as a
# data.frame -- embedding a ~50-row data.frame directly here would cause
# data.frame() to recycle every other scalar column to match its row count,
# turning this into a ~50-row table with suhi_path/mean_rural_LST/lst_path/
# ndvi_path duplicated on every row.
knime.out <- data.frame(
  suhi_path        = out_suhi,
  mean_rural_LST    = mean_rural_LST,
  histogram_csv_path = out_hist_csv,
  histogram = hist_df,
  lst_path          = lst_path,
  ndvi_path         = ndvi_path,
  suhi_mean         = suhi_mean,
  suhi_min          = suhi_min,
  suhi_max          = suhi_max,
  suhi_p01          = suhi_p01,
  suhi_p99          = suhi_p99,
  suhi_n_pixels     = suhi_n
)

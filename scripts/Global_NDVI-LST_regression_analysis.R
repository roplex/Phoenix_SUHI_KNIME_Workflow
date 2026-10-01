# ============================================================================
# KNIME R Snippet -- Node: "Global NDVI-LST Regression"
# Chapter 7 (SUHI) -- corrected version
#
# Fits one NDVI~LST regression across the entire city-clipped raster (as
# opposed to the Stratified node, which fits one regression per village).
# Feeds Table 7.2's "effect on regression R2" comparison once this same
# node is also run against the bilinear-resampled LST -- see NOTE below.
#
# Inputs:
#   knime.in[["lst_path"]]   -- city-clipped LST (TsHARP-sharpened)
#   knime.in[["ndvi_path"]]  -- city-clipped NDVI, passed through
#   knime.in[["suhi_path"]]  -- SUHI raster path; used here to reuse the
#                                SUHI Computation node's own extreme-pixel
#                                threshold for the robustness check below
#
# Outputs (one row per sampled scatter point, stats replicated across rows
# so a downstream Scatter Plot node has slope/intercept alongside each
# point for a regression-line overlay -- deliberate, unlike the earlier
# SUHI histogram issue where row-recycling was accidental):
#   slope, intercept, r_squared, pearson_r, r_squared_excl_extreme,
#   n_excluded_extreme, n_points, suhi_path, lst_path, ndvi_path, NDVI, LST
#
# NOTE ON TABLE 7.2: this node reports R2 for whichever lst_path is fed in
# -- currently the TsHARP output. To fill in Table 7.2's bilinear-vs-TsHARP
# R2 comparison, this same node also needs to run against the ORIGINAL
# bilinear-resampled LST, routed through the identical downstream chain
# (Clipping & Masking -> Urban vs Rural Mask -> SUHI Computation -> here) so
# only the resampling method differs. The LST Resampling node's bilinear
# comparison raster (lst_bilinear) was computed in-memory for its own RMSE/
# bias diagnostic but was never written to disk -- it will need a
# writeRaster() added there and to be pushed through this same chain before
# Table 7.2 can be completed.
# ============================================================================

library(terra)

# ---- Inputs from KNIME ----
lst_path  <- as.character(knime.in[["lst_path"]])
ndvi_path <- as.character(knime.in[["ndvi_path"]])
suhi_path <- as.character(knime.in[["suhi_path"]])

stopifnot(file.exists(lst_path), file.exists(ndvi_path), file.exists(suhi_path))

# ---- Load rasters ----
lst  <- rast(lst_path)[[1]]
ndvi <- rast(ndvi_path)[[1]]
suhi <- rast(suhi_path)[[1]]

# ---- Align CRS ----
if (!identical(crs(lst), crs(ndvi))) ndvi <- project(ndvi, crs(lst))
if (!identical(crs(lst), crs(suhi))) suhi <- project(suhi, crs(lst))

# ---- Align extent/resolution ----
# compareGeom() replaces the original all.equal(ext(lst), ext(ndvi)) check:
# all.equal() on SpatExtent objects returns a descriptive STRING when the
# extents differ, not FALSE -- !all.equal(...) then coerces that string via
# as.logical(), which returns NA, not TRUE, so the resample step could be
# silently skipped (or the node could error with "missing value where
# TRUE/FALSE needed") exactly when it's needed most. compareGeom() is the
# correct terra-native check and is used consistently elsewhere in this
# pipeline (Urban vs Rural Mask, DEM alignment).
if (!compareGeom(lst, ndvi, stopOnError = FALSE)) ndvi <- resample(ndvi, lst, method = "bilinear")
if (!compareGeom(lst, suhi, stopOnError = FALSE)) suhi <- resample(suhi, lst, method = "bilinear")

# ---- Extract paired values, filtering complete cases across ALL THREE
# together so LST/NDVI/SUHI stay aligned to the same pixels throughout --
# the same discipline established in the Urban vs Rural node after an
# earlier index-misalignment bug there (lm()'s internal na.omit silently
# shifting positions between related vectors). ----
v_lst  <- as.numeric(values(lst))
v_ndvi <- as.numeric(values(ndvi))
v_suhi <- as.numeric(values(suhi))

df <- data.frame(LST = v_lst, NDVI = v_ndvi, SUHI = v_suhi)
df <- df[complete.cases(df), ]

if (nrow(df) < 5) {
  stop("Too few valid pixel pairs for regression")
}

# ---- Regression + correlation (full sample) ----
fit <- lm(LST ~ NDVI, data = df)
r2  <- summary(fit)$r.squared
pearson_r <- suppressWarnings(cor(df$NDVI, df$LST))
cat(sprintf("Global regression (n=%d): slope=%.4f intercept=%.2f R2=%.4f Pearson r=%.4f\n",
            nrow(df), coef(fit)[2], coef(fit)[1], r2, pearson_r))

# ---- Robustness check: exclude the previously-flagged ambiguous
# extreme-SUHI pixels (SUHI < -8 or > 4 degC -- the same threshold used in
# the SUHI Computation node's clustering diagnostic, reused here rather than
# a fresh generic outlier rule so this stays tied to the same 18 pixels
# already identified) and see whether R2 shifts meaningfully. ----
extreme_idx <- which(df$SUHI < -8 | df$SUHI > 4)
n_excluded_extreme <- length(extreme_idx)
if (n_excluded_extreme > 0) {
  fit_clean <- lm(LST ~ NDVI, data = df[-extreme_idx, ])
  r2_excl_extreme <- summary(fit_clean)$r.squared
  cat(sprintf("Excluding %d extreme-SUHI pixels (%.2f%%): R2=%.4f (full-sample R2=%.4f, shift=%.4f)\n",
              n_excluded_extreme, 100 * n_excluded_extreme / nrow(df),
              r2_excl_extreme, r2, r2_excl_extreme - r2))
} else {
  r2_excl_extreme <- r2
  cat("No pixels matched the SUHI < -8 or > 4 degC threshold in this sample; R2 unaffected.\n")
}

# ---- Sample data for scatterplot ----
sample_size <- min(5000, nrow(df))
set.seed(42)

df_sample <- df[sample.int(nrow(df), size = sample_size), ]

# Use [1] to ensure lst_path is a single scalar path, not a multi-element vector
target_dir <- dirname(as.character(lst_path)[1])
out_csv    <- file.path(target_dir, "regression_sample_points.csv")

write.csv(df_sample, out_csv, row.names = FALSE)

# ---- Replicate stats across rows of df_sample ----
n <- nrow(df_sample)
knime.out <- data.frame(
  slope                  = rep(as.numeric(coef(fit)[2]), n),
  intercept              = rep(as.numeric(coef(fit)[1]), n),
  r_squared              = rep(as.numeric(r2), n),
  pearson_r              = rep(as.numeric(pearson_r), n),
  r_squared_excl_extreme = rep(as.numeric(r2_excl_extreme), n),
  n_excluded_extreme     = rep(as.integer(n_excluded_extreme), n),
  n_points               = rep(as.integer(nrow(df)), n),
  suhi_path              = rep(suhi_path, n),
  lst_path               = rep(lst_path, n),
  ndvi_path              = rep(ndvi_path, n),
  NDVI                   = df_sample$NDVI,
  LST                    = df_sample$LST
)

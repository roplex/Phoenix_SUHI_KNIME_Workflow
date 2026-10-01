# ============================================================================
# KNIME R Snippet -- Node: "LST Resampling"
# Replaces bilinear interpolation with TsHARP thermal sharpening.
#
# Method: Agam, N., Kustas, W. P., Anderson, M. C., Li, F., & Neale, C. M. U.
# (2007). A vegetation index based technique for spatial sharpening of
# thermal imagery. Remote Sensing of Environment, 107(4), 545-558.
#
# Inputs (unchanged from the original node):
#   knime.in[["lst_composite_path"]]   -- e.g., superrelaxed LST composite
#   knime.in[["ndvi_composite_path"]]  -- strict NDVI composite (reference grid)
#
# Outputs (same column names as the original bilinear node, so every
# downstream node -- Clipping & Masking, Urban vs Rural Mask, SUHI
# Computation -- needs no changes at all):
#   resampled_lst_path, ndvi_aligned_path
# Plus new diagnostic columns feeding Table 7.2 directly:
#   bilinear_rmse, bilinear_bias, tsharp_rmse, tsharp_bias
# ============================================================================

library(terra)

# ---- Inputs from KNIME ----
lst_path  <- as.character(knime.in[["lst_composite_path"]])
ndvi_path <- as.character(knime.in[["ndvi_composite_path"]])

# ---- Load rasters ----
lst  <- rast(lst_path)
ndvi <- rast(ndvi_path)

# ----------------------------------------------------------------------------
# TsHARP disaggregation function
# ----------------------------------------------------------------------------
tsharp_disaggregate <- function(lst_coarse, ndvi_fine, window_size = 9) {

  agg_factor <- round(res(lst_coarse)[1] / res(ndvi_fine)[1])

  # Step 1: aggregate NDVI to the coarse LST grid so the regression is fit
  # at the same spatial support as the original LST observations.
  ndvi_coarse <- aggregate(ndvi_fine, fact = agg_factor, fun = "mean", na.rm = TRUE)
  ndvi_coarse <- resample(ndvi_coarse, lst_coarse, method = "near")

  fit_global <- NULL

  if (is.null(window_size)) {
    # --- TsHARP-global: single regression across the full study area ---
    df <- data.frame(lst = values(lst_coarse)[, 1], ndvi = values(ndvi_coarse)[, 1])
    df <- df[complete.cases(df), ]
    fit_global <- lm(lst ~ ndvi, data = df)
    a <- coef(fit_global)[1]; b <- coef(fit_global)[2]

    lst_predicted_fine   <- a + b * ndvi_fine
    lst_predicted_coarse <- a + b * ndvi_coarse

  } else {
    # --- TsHARP-local: moving-window regression, one (a, b) pair per window ---
    fw <- matrix(1, window_size, window_size)

    # NOTE: terra's focalPairs() forwards na.rm (and possibly other internal
    # arguments) straight through to `fun` itself, not just to the window
    # extraction step -- so `fun` must accept them even if it ignores them.
    # `...` absorbs whatever terra passes; NA-pairing is then handled
    # explicitly below rather than relied upon implicitly.
    #
    # IMPORTANT: focalPairs() passes window values to `fun` positionally, in
    # the SAME ORDER the layers were stacked into x = c(lst_coarse, ndvi_coarse)
    # above -- i.e., first argument = LST values, second argument = NDVI
    # values, regardless of what you name those parameters. To make this
    # impossible to silently invert (as happened once already), the
    # parameters below are named for what they physically are rather than
    # generic x/y -- do not rename them to x/y or swap their order.
    coeffs <- focalPairs(
      x = c(lst_coarse, ndvi_coarse), w = fw, na.rm = TRUE,
      fun = function(lst_vals, ndvi_vals, ...) {
        ok <- stats::complete.cases(lst_vals, ndvi_vals)
        lst_vals <- lst_vals[ok]; ndvi_vals <- ndvi_vals[ok]
        # Require a handful of valid, non-identical pairs before trusting a
        # local fit -- a 9x9 window has at most 81 cells, but edges, masked
        # pixels, and QC gaps can leave far fewer valid pairs than that.
        if (length(ndvi_vals) < 5 || length(unique(ndvi_vals)) < 2) return(c(NA_real_, NA_real_))
        # LST as the response, NDVI as the predictor -- must match the
        # global branch's lm(lst ~ ndvi) above for the two branches to be
        # comparable.
        fit <- lm(lst_vals ~ ndvi_vals)
        c(coef(fit)[1], coef(fit)[2])
      }
    )
    a_coarse <- coeffs[[1]]; b_coarse <- coeffs[[2]]
    a_fine <- resample(a_coarse, ndvi_fine, method = "bilinear")
    b_fine <- resample(b_coarse, ndvi_fine, method = "bilinear")

    lst_predicted_fine   <- a_fine + b_fine * ndvi_fine
    lst_predicted_coarse <- a_coarse + b_coarse * ndvi_coarse
  }

  # Step 2: residual at coarse scale = observed - regression-predicted.
  residual_coarse <- lst_coarse - lst_predicted_coarse

  # Step 3: smooth the residual to fine resolution and add it back, so the
  # sharpened output still aggregates to match the original coarse LST
  # observation (energy conservation) -- this is what makes it a genuine
  # disaggregation rather than a smoothing operation.
  residual_fine <- resample(residual_coarse, ndvi_fine, method = "bilinear")
  lst_sharpened_fine <- lst_predicted_fine + residual_fine
  names(lst_sharpened_fine) <- "LST_sharpened_250m"

  list(
    sharpened  = lst_sharpened_fine,
    predicted  = lst_predicted_fine,
    residual   = residual_fine,
    regression = fit_global
  )
}

# ----------------------------------------------------------------------------
# Sharpening-quality diagnostic: re-aggregate back to coarse resolution and
# compare against the original observed coarse LST. Feeds Table 7.2 directly.
# ----------------------------------------------------------------------------
diagnose_sharpening_quality <- function(lst_coarse_observed, sharpened_fine, method_label) {
  fact <- round(res(lst_coarse_observed)[1] / res(sharpened_fine)[1])
  reaggregated <- aggregate(sharpened_fine, fact = fact, fun = "mean", na.rm = TRUE)
  reaggregated <- resample(reaggregated, lst_coarse_observed, method = "near")
  resid <- values(lst_coarse_observed)[, 1] - values(reaggregated)[, 1]
  resid <- resid[!is.na(resid)]
  rmse <- sqrt(mean(resid^2))
  bias <- mean(resid)
  cat(sprintf("[%s] RMSE vs. original coarse LST: %.4f degC | Mean bias: %.4f degC | n = %d\n",
              method_label, rmse, bias, length(resid)))
  list(rmse = rmse, bias = bias, n = length(resid))
}

# ---- Run TsHARP (replaces bilinear resample) ----
tsharp_result <- tsharp_disaggregate(lst, ndvi, window_size = 9)
lst_resampled <- tsharp_result$sharpened

# ---- Diagnostic comparison: bilinear (original) vs. TsHARP (revised) ----
lst_bilinear  <- resample(lst, ndvi, method = "bilinear")
diag_bilinear <- diagnose_sharpening_quality(lst, lst_bilinear, "Bilinear (original)")
diag_tsharp   <- diagnose_sharpening_quality(lst, lst_resampled, "TsHARP (revised)")

# ---- Save resampled LST ----
out_lst  <- gsub(".tif", "_tsharp.tif", lst_path)
out_lst_bilinear <- gsub(".tif", "_bilinear.tif", lst_path)

# ---- Save NDVI aligned copy in same folder ----
out_ndvi <- file.path(dirname(out_lst), "NDVI_composite_aligned.tif")
writeRaster(ndvi, out_ndvi, overwrite = TRUE)
writeRaster(lst_resampled, out_lst, overwrite = TRUE)
writeRaster(lst_bilinear, out_lst_bilinear, overwrite = TRUE)

# ---- Return results to KNIME ----
knime.out <- data.frame(
  resampled_lst_path = out_lst,
  ndvi_aligned_path  = out_ndvi,
  bilinear_rmse = diag_bilinear$rmse,
  bilinear_bias = diag_bilinear$bias,
  tsharp_rmse   = diag_tsharp$rmse,
  tsharp_bias   = diag_tsharp$bias
)

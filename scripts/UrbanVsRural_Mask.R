# ============================================================================
# KNIME R Snippet -- Node: "Urban vs Rural Mask using LST" (Pre-processed)
# Chapter 7 (SUHI) -- consolidated, corrected version
# ============================================================================
library(terra)

# ---- Inputs from KNIME ----
lst_clipped_path  <- as.character(knime.in[["lst_clipped_path"]])
ndvi_clipped_path <- as.character(knime.in[["ndvi_clipped_path"]])
lst_wide_path     <- as.character(knime.in[["lst_wide_path"]])

# Local pre-processed assets
phx_ua_path  <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Preprocessed/phoenix_urban_area.gpkg"
wc_mask_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Preprocessed/worldcover_non_urban_aligned.tif"
dem_path     <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Preprocessed/dem_aligned.tif"

stopifnot(
  file.exists(lst_clipped_path),
  file.exists(lst_wide_path),
  file.exists(phx_ua_path),
  file.exists(wc_mask_path),
  file.exists(dem_path)
)

# ---- Load Wide LST & Pre-processed Inputs ----
lst_wide <- rast(lst_wide_path)[[1]]
phx_ua   <- project(vect(phx_ua_path), crs(lst_wide))
wc_mask  <- rast(wc_mask_path)
dem_wide <- rast(dem_path)

if (!compareGeom(lst_wide, wc_mask, stopOnError = FALSE)) {
  stop("lst_wide_path grid does not match the pre-processed WorldCover mask grid -- rebuild wc_mask against the current lst_wide before proceeding.")
}
if (!compareGeom(lst_wide, dem_wide, stopOnError = FALSE)) {
  stop("dem_aligned.tif grid does not match lst_wide grid -- rebuild dem_aligned.tif against the current lst_wide before proceeding.")
}

# ---- Extent-Aware Radius Sizing ----
grid_ext <- ext(lst_wide)
ua_ext   <- ext(phx_ua)
lat0     <- (as.numeric(grid_ext$ymin) + as.numeric(grid_ext$ymax)) / 2
km_per_deg_lon <- 111.32 * cos(lat0 * pi / 180)
km_per_deg_lat <- 111.0

margin_km <- c(
  west  = (as.numeric(ua_ext$xmin) - as.numeric(grid_ext$xmin)) * km_per_deg_lon,
  east  = (as.numeric(grid_ext$xmax) - as.numeric(ua_ext$xmax)) * km_per_deg_lon,
  south = (as.numeric(ua_ext$ymin) - as.numeric(grid_ext$ymin)) * km_per_deg_lat,
  north = (as.numeric(grid_ext$ymax) - as.numeric(ua_ext$ymax)) * km_per_deg_lat
)

# ---- Search Extent: buffer at the intended radius, let crop() truncate per-direction ----
# (Do NOT scale the buffer width by the tightest margin -- that forces every side down to
# match the worst side, even ones with real room. Buffer at the full design radius and let
# the raster's own extent truncate wherever it falls short.)
search_radius_km <- 60
search_buffer <- buffer(phx_ua, width = search_radius_km * 1000)
search_extent <- crop(search_buffer, grid_ext)  # exact geometric clip to the raster extent

cat(sprintf("Requested radius: %d km | margin (km): west=%.1f east=%.1f south=%.1f north=%.1f\n",
            search_radius_km, margin_km["west"], margin_km["east"], margin_km["south"], margin_km["north"]))

# ---- Rural Reference Mask Extraction ----
rural_mask <- mask(wc_mask, phx_ua, inverse = TRUE)
# mask = TRUE restricts to the *actual* buffered polygon shape, not just its bounding box --
# otherwise pixels in the far corners of the bounding rectangle (beyond the true ring) leak in.
rural_mask <- crop(rural_mask, search_extent, mask = TRUE)
names(rural_mask) <- "rural_reference_mask"

# Align lst_wide extent with rural_mask extent before masking
lst_wide_cropped <- crop(lst_wide, rural_mask)

# ---- Diagnostics & Mean Computation ----
n_valid        <- global(!is.na(rural_mask), "sum", na.rm = TRUE)[[1]]
mean_rural_LST <- global(mask(lst_wide_cropped, rural_mask), "mean", na.rm = TRUE)[[1]]
if (is.na(mean_rural_LST)) {
  stop("No valid rural LST pixels found under redefined mask.")
}

# ---- Directional breakdown diagnostic ----
ua_centroid <- centroids(phx_ua)
cx <- crds(ua_centroid)[1, 1]; cy <- crds(ua_centroid)[1, 2]
rural_pts <- as.points(rural_mask, values = FALSE, na.rm = TRUE)
pt_crds   <- crds(rural_pts)
bearing   <- (atan2(pt_crds[, 1] - cx, pt_crds[, 2] - cy) * 180 / pi) %% 360
sector    <- cut(bearing, breaks = c(-1, 45, 135, 225, 315, 361),
                  labels = c("N", "E", "S", "W", "N"))
lst_vals  <- extract(lst_wide_cropped, pt_crds)[, 1]

cat("Mean LST by sector (deg C):\n")
print(tapply(lst_vals, sector, mean, na.rm = TRUE))
cat("Pixel count by sector:\n")
print(tapply(lst_vals, sector, length))

# ---- Elevation Covariate Check ----
dem_wide_cropped <- crop(dem_wide, rural_mask)
elev_vals <- extract(dem_wide_cropped, pt_crds)[, 1]

# ---- Filter to complete cases ONCE, keeping the SAME variable names ----
# lm() silently drops NA rows internally; filtering here -- in place, no
# renaming -- keeps every downstream block (Option A, Option B, tolerance
# sensitivity, linearity check, outlier lookup) aligned automatically, with
# no risk of a name or length mismatch anywhere further down.
complete_idx <- which(complete.cases(lst_vals, elev_vals))
n_dropped <- length(lst_vals) - length(complete_idx)
if (n_dropped > 0) {
  cat(sprintf("Dropping %d incomplete pixels (NA in LST or elevation) before modeling.\n", n_dropped))
  lst_vals  <- lst_vals[complete_idx]
  elev_vals <- elev_vals[complete_idx]
  sector    <- sector[complete_idx]
  pt_crds   <- pt_crds[complete_idx, , drop = FALSE]
}

cat("Mean elevation by sector (m):\n")
print(tapply(elev_vals, sector, mean, na.rm = TRUE))
elev_lst_cor <- cor(elev_vals, lst_vals, use = "complete.obs")
cat(sprintf("Correlation between elevation and rural LST: %.3f\n", elev_lst_cor))

# Original municipal-boundary baseline for Table 7.3 comparison
mean_rural_LST_original <- 44.8

# ---- Urban Core Elevation (comparison target) ----
urban_elev_vals <- unlist(extract(dem_wide_cropped, phx_ua)[, -1])
urban_elev_mean <- mean(urban_elev_vals, na.rm = TRUE)
cat(sprintf("Urban core mean elevation: %.1f m (n=%d valid pixels)\n",
            urban_elev_mean, sum(!is.na(urban_elev_vals))))

# ---- Option A: Elevation-Matched Rural Sample ----
elev_tolerance <- 100  # meters; widen if too few pixels survive
elev_matched_idx <- which(abs(elev_vals - urban_elev_mean) <= elev_tolerance)
n_elev_matched <- length(elev_matched_idx)

if (n_elev_matched < 500) {
  warning(sprintf("Only %d pixels survive elevation matching at +/-%dm -- consider widening elev_tolerance.",
                   n_elev_matched, elev_tolerance))
}

mean_rural_LST_elev_matched <- if (n_elev_matched > 0) {
  mean(lst_vals[elev_matched_idx], na.rm = TRUE)
} else {
  NA_real_
}

cat(sprintf("Option A -- Elevation-matched rural LST (+/-%dm of urban %.1fm): %.2f degC (n=%d)\n",
            elev_tolerance, urban_elev_mean, mean_rural_LST_elev_matched, n_elev_matched))
cat("Sector composition of elevation-matched sample:\n")
print(table(sector[elev_matched_idx]))

# ---- Option B: Elevation-Adjusted (Detrended) Rural LST ----
elev_lst_model <- lm(lst_vals ~ elev_vals)
model_summary  <- summary(elev_lst_model)
cat(sprintf("Elevation-LST regression: slope=%.4f degC/m, R2=%.3f\n",
            coef(elev_lst_model)[2], model_summary$r.squared))

mean_rural_LST_elev_adjusted <- as.numeric(predict(
  elev_lst_model, newdata = data.frame(elev_vals = urban_elev_mean)
))
cat(sprintf("Option B -- Elevation-adjusted rural LST (predicted at urban elevation %.1fm): %.2f degC\n",
            urban_elev_mean, mean_rural_LST_elev_adjusted))

# ---- Tolerance Sensitivity Check (Option A) ----
tolerance_grid <- c(50, 75, 100, 150, 200)
sensitivity_results <- data.frame(
  tolerance_m    = tolerance_grid,
  mean_rural_LST = NA_real_,
  n_pixels       = NA_integer_
)
for (i in seq_along(tolerance_grid)) {
  idx <- which(abs(elev_vals - urban_elev_mean) <= tolerance_grid[i])
  sensitivity_results$mean_rural_LST[i] <- mean(lst_vals[idx], na.rm = TRUE)
  sensitivity_results$n_pixels[i]       <- length(idx)
}
cat("Tolerance sensitivity (Option A):\n")
print(sensitivity_results)

# ---- Linearity Check (Option B) ----
elev_lst_model_quad <- lm(lst_vals ~ elev_vals + I(elev_vals^2))
anova_result <- anova(elev_lst_model, elev_lst_model_quad)
cat("Linear vs quadratic model comparison:\n")
print(anova_result)

cat("Linear model residual summary (deg C):\n")
print(summary(residuals(elev_lst_model)))

pred_linear <- as.numeric(predict(elev_lst_model,      newdata = data.frame(elev_vals = urban_elev_mean)))
pred_quad   <- as.numeric(predict(elev_lst_model_quad, newdata = data.frame(elev_vals = urban_elev_mean)))
cat(sprintf("Linear: %.2f degC | Quadratic: %.2f degC | Difference: %.2f degC\n",
            pred_linear, pred_quad, pred_quad - pred_linear))

# If running interactively in RStudio (not headless KNIME), this renders a
# visual check -- linear fit vs. a flexible local (lowess) fit:
# plot(elev_vals, lst_vals, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.1),
#      xlab = "Elevation (m)", ylab = "Rural LST (deg C)")
# abline(elev_lst_model, col = "red", lwd = 2)
# lines(lowess(elev_vals, lst_vals), col = "blue", lwd = 2)
# legend("topright", c("Linear fit", "Lowess fit"), col = c("red", "blue"), lwd = 2)

# ---- Outlier / Cold-Tail Diagnostics (aligned -- see complete-case filter above) ----
resids <- residuals(elev_lst_model)
outlier_idx <- which.min(resids)
cat(sprintf("Worst residual: %.2f degC | elevation %.1fm | observed LST %.2f degC | sector %s\n",
            resids[outlier_idx], elev_vals[outlier_idx], lst_vals[outlier_idx],
            as.character(sector[outlier_idx])))

cold_tail_idx <- which(resids < -5)
cat(sprintf("Pixels with residual < -5 degC: %d (%.2f%% of sample)\n",
            length(cold_tail_idx), 100 * length(cold_tail_idx) / length(resids)))
cat("Sector composition of cold-tail pixels:\n")
print(table(sector[cold_tail_idx]))
cat(sprintf("Cold-tail elevation range: %.1f - %.1f m (mean %.1f m)\n",
            min(elev_vals[cold_tail_idx]), max(elev_vals[cold_tail_idx]), mean(elev_vals[cold_tail_idx])))

# ---- Robustness Check: Option A excluding cold-tail pixels ----
# The cold tail's mean elevation (390.3m in the prior run) sits inside the
# elevation-matching window around the urban core, so some of these pixels
# are likely already inside the Option A sample, quietly pulling it down.
elev_matched_clean_idx <- setdiff(elev_matched_idx, cold_tail_idx)
n_elev_matched_clean <- length(elev_matched_clean_idx)
mean_rural_LST_elev_matched_clean <- if (n_elev_matched_clean > 0) {
  mean(lst_vals[elev_matched_clean_idx], na.rm = TRUE)
} else {
  NA_real_
}
n_cold_tail_in_matched <- n_elev_matched - n_elev_matched_clean
cat(sprintf("Cold-tail pixels inside the elevation-matched sample: %d (%.2f%% of matched sample)\n",
            n_cold_tail_in_matched, 100 * n_cold_tail_in_matched / n_elev_matched))
cat(sprintf("Option A excluding cold-tail: %.2f degC (n=%d) vs %.2f degC unfiltered (n=%d) | shift: %.2f degC\n",
            mean_rural_LST_elev_matched_clean, n_elev_matched_clean,
            mean_rural_LST_elev_matched, n_elev_matched,
            mean_rural_LST_elev_matched_clean - mean_rural_LST_elev_matched))

# ---- What land cover are the cold-tail pixels actually? ----
# wc_mask only stores a binary non-urban flag (classes 20/30/40/60 already
# collapsed to 1); to see WHICH of those classes the cold-tail pixels are,
# read the original multi-class WorldCover raster directly.
worldcover_raw_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/WorldCover/ESA_WorldCover_10m_2021_v200_N33W114_Map.tif"
stopifnot(file.exists(worldcover_raw_path))
worldcover_raw <- rast(worldcover_raw_path)

cold_tail_pts_wc <- project(vect(pt_crds[cold_tail_idx, , drop = FALSE], crs = crs(lst_wide)),
                             crs(worldcover_raw))
cold_tail_wc_class <- extract(worldcover_raw, cold_tail_pts_wc)[, 2]

cat("WorldCover class composition of cold-tail pixels (10=Tree 20=Shrub 30=Grass 40=Cropland 60=Bare):\n")
print(table(cold_tail_wc_class))
cat("WorldCover class by sector, cold-tail pixels only:\n")
print(table(sector[cold_tail_idx], cold_tail_wc_class))

# ---- Pass Outputs to Downstream KNIME Nodes ----
knime.out <- data.frame(
  lst_path                          = lst_clipped_path,
  ndvi_path                         = ndvi_clipped_path,
  mean_rural_LST_original           = mean_rural_LST_original,
  mean_rural_LST_full_sample        = mean_rural_LST,
  rural_pixels_found                = n_valid,
  rural_search_radius_km_requested  = search_radius_km,
  margin_west_km                    = margin_km["west"],
  margin_east_km                    = margin_km["east"],
  margin_south_km                   = margin_km["south"],
  margin_north_km                   = margin_km["north"],
  elev_lst_correlation              = elev_lst_cor,
  urban_core_mean_elev_m            = urban_elev_mean,
  mean_rural_LST_elev_matched       = mean_rural_LST_elev_matched,
  n_elev_matched_pixels             = n_elev_matched,
  mean_rural_LST_elev_matched_clean = mean_rural_LST_elev_matched_clean,
  n_elev_matched_clean_pixels       = n_elev_matched_clean,
  elev_tolerance_m                  = elev_tolerance,
  mean_rural_LST_elev_adjusted      = mean_rural_LST_elev_adjusted,
  elev_lst_regression_r2            = model_summary$r.squared
)

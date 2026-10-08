# ============================================================================
# Chapter 7 (SUHI) -- Sensitivity analysis of the rural reference-zone correction
# Standalone RStudio script (runs OUTSIDE the KNIME workflow)
#
# WHY THIS EXISTS
# The redefinition of the rural reference zone moved the rural-mean LST from
# 44.8 degC (municipal-boundary mask) to 47.89 degC (Census Urbanized Area +
# ESA WorldCover + elevation-matched to the urban core). A shift of that size
# is only credible if it is robust to the reasonable analyst choices made along
# the way. This script re-computes the baseline under systematic alternatives
# and reports how far the baseline -- and the SUHI conclusions that depend on
# it -- move.
#
# FIVE FAMILIES OF ALTERNATIVES (one per point raised in review)
#   A. Elevation handling      : matching tolerance, one-sided windows, absolute
#                                elevation bands, regression adjustment (linear
#                                and quadratic), histogram re-weighting
#   B. Land-cover exclusions   : alternative WorldCover class sets and purity /
#                                built-up-fraction thresholds
#   C. Spatial buffers         : inner exclusion buffer around the Urbanized
#                                Area x outer search radius
#   D. Retained-pixel sample   : number and spatial (N/E/S/W) distribution of
#                                retained pixels, leave-one-sector-out,
#                                sector-balanced mean, spatial block bootstrap
#   E. Temporal sampling       : per-date and per-month baselines from the
#                                individual 8-day MOD11A2 scenes versus the
#                                median composite used in the chapter
#
# OUTPUTS (written to <root>/Sensitivity/)
#   sensitivity_summary.csv      every scenario, one row each
#   sensitivity_decision.txt     plain-language verdict + numbers to quote
#   sensitivity_forest_plot.png  baseline by scenario against the primary value
#   sensitivity_temporal.png     per-date baseline and city-wide SUHI
#   sessionInfo.txt              package versions for reproducibility
#
# HOW TO RUN (RStudio)
#   1. install.packages("terra")            # the only required package
#   2. Check the PATHS block below (defaults follow your PhoenixData2 layout).
#   3. Source the whole file. The first run builds a WorldCover class-fraction
#      raster from the 10 m tile (a few minutes, a few GB of temporary disk) and
#      caches it in Preprocessed/; later runs reuse the cache.
#   4. Read the console log; it first checks that the primary scenario
#      reproduces the chapter's 47.89 degC before reporting anything else.
# ============================================================================

suppressPackageStartupMessages(library(terra))
set.seed(42)

# ------------------------------- PATHS --------------------------------------
root <- Sys.getenv("PHX_ROOT", "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2")

P <- list(
  # 1 km LST composite used to compute the rural baseline (the "lst_wide" input
  # of the KNIME "Urban vs Rural Mask" node). Change if your flow variable
  # points at a different file; the baseline check below will tell you.
  lst_wide   = file.path(root, "LST/LST_composite_superrelaxed.tif"),
  # Sharpened LST used for the urban side of SUHI (optional; falls back to lst_wide)
  lst_urban  = file.path(root, "LST/LST_composite_superrelaxed_tsharp.tif"),
  dem        = file.path(root, "Preprocessed/dem_aligned.tif"),
  wc_binary  = file.path(root, "Preprocessed/worldcover_non_urban_aligned.tif"),
  ua         = file.path(root, "Preprocessed/phoenix_urban_area.gpkg"),
  city       = file.path(root, "Boundary/City_Limit_Light_Outline.geojson"),
  villages   = file.path(root, "Boundary/Villages.geojson"),
  wc_raw     = file.path(root, "WorldCover/ESA_WorldCover_10m_2021_v200_N33W114_Map.tif"),
  wc_frac    = file.path(root, "Preprocessed/worldcover_class_fractions.tif"),  # cache
  lst_dir    = file.path(root, "LST"),
  lst_scene_pattern = "^MOD11A2\\.061_LST_Day_1km_doy2024\\d{3}000000_aid0001_clean_superrelaxed\\.tif$",
  out_dir    = file.path(root, "Sensitivity")
)

CFG <- list(
  primary_tol_m      = 100,    # chapter's primary elevation-matching tolerance
  expected_baseline  = 47.89,  # value reported in Table 7.3
  baseline_tolerance = 0.15,   # how close the reproduction must be (degC)
  search_radius_km   = 60,     # chapter's outer search radius
  min_pixels         = 200,    # scenarios with fewer retained pixels are flagged
  wc_purity          = 0.50,   # default share of an LST cell that must be an allowed class
  block_km           = 10,     # block size for the spatial bootstrap
  n_boot             = 2000,
  robust_threshold_c = 1.0     # spread (degC) below which the baseline is called robust
)

dir.create(P$out_dir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ------------------------------- LOAD ---------------------------------------
need <- c("lst_wide", "dem", "wc_binary", "ua", "city", "villages")
miss <- need[!file.exists(unlist(P[need]))]
if (length(miss)) stop("Missing input file(s): ", paste(unlist(P[miss]), collapse = "\n  "))

lst <- rast(P$lst_wide)[[1]]
crs_lst <- crs(lst)
ua      <- project(vect(P$ua),       crs_lst)
city    <- project(vect(P$city),     crs_lst)
vill    <- project(vect(P$villages), crs_lst)

align_to_lst <- function(r, method) {
  if (compareGeom(lst, r, stopOnError = FALSE)) r else resample(r, lst, method = method)
}
dem <- align_to_lst(rast(P$dem)[[1]], "bilinear")
wcb <- align_to_lst(rast(P$wc_binary)[[1]], "near")

urban_lst <- if (file.exists(P$lst_urban)) rast(P$lst_urban)[[1]] else lst
if (!file.exists(P$lst_urban)) log_msg("NOTE: sharpened LST not found; city-wide SUHI uses lst_wide.")

# ------------------------- CELL TABLE (one row per LST cell) -----------------
# Cells touched by the Urbanized Area polygon are excluded, exactly as the KNIME
# node does with mask(..., inverse = TRUE) (terra's default touches = TRUE).
ua_r   <- rasterize(ua, lst, field = 1, touches = TRUE)
dist_r <- distance(ua_r)                      # metres to nearest UA cell (0 inside)
in_ua  <- !is.na(values(ua_r)[, 1])

xy   <- xyFromCell(lst, seq_len(ncell(lst)))
ctr  <- crds(centroids(ua))[1, ]
bearing <- (atan2(xy[, 1] - ctr[1], xy[, 2] - ctr[2]) * 180 / pi) %% 360
sector  <- c("N", "E", "S", "W", "N")[findInterval(bearing, c(0, 45, 135, 225, 315, 360))]

lat0   <- mean(range(xy[, 2]))
x_km   <- (xy[, 1] - min(xy[, 1])) * 111.32 * cos(lat0 * pi / 180)
y_km   <- (xy[, 2] - min(xy[, 2])) * 111.0

px <- data.frame(
  cell   = seq_len(ncell(lst)),
  lst    = values(lst)[, 1],
  elev   = values(dem)[, 1],
  wcb    = values(wcb)[, 1],
  in_ua  = in_ua,
  dist_km = values(dist_r)[, 1] / 1000,
  sector = sector,
  bx = floor(x_km / CFG$block_km), by = floor(y_km / CFG$block_km)
)
px$block <- paste(px$bx, px$by)
px <- px[is.finite(px$lst) & is.finite(px$elev), ]

# ------------------------- WORLDCOVER CLASS FRACTIONS ------------------------
# Fraction of each LST cell covered by each WorldCover class
# (10 tree, 20 shrub, 30 grass, 40 crop, 50 built-up, 60 bare/sparse, 80 water, 90 wetland)
wc_classes <- c(10, 20, 30, 40, 50, 60, 80, 90)
have_frac <- FALSE
if (file.exists(P$wc_frac)) {
  fr <- align_to_lst(rast(P$wc_frac), "bilinear"); have_frac <- TRUE
  log_msg("Loaded cached WorldCover class fractions: %s", P$wc_frac)
} else if (file.exists(P$wc_raw)) {
  log_msg("Building WorldCover class fractions from the 10 m tile (one-off, a few minutes)...")
  wc_raw <- rast(P$wc_raw)
  e <- as.polygons(ext(lst), crs = crs_lst) |> project(crs(wc_raw)) |> ext()
  wc_raw <- crop(wc_raw, e)
  seg <- segregate(wc_raw, classes = wc_classes, other = 0)
  f <- max(1, floor(min(res(lst) / res(wc_raw)) ))
  if (f > 1) seg <- aggregate(seg, fact = f, fun = "mean", na.rm = TRUE)
  fr <- resample(seg, lst, method = "average")
  names(fr) <- paste0("fr_", wc_classes)
  writeRaster(fr, P$wc_frac, overwrite = TRUE)
  have_frac <- TRUE
} else {
  log_msg("NOTE: no WorldCover raw tile or cache found -- land-cover family (B) will be limited.")
}
if (have_frac) {
  names(fr) <- paste0("fr_", wc_classes)
  fv <- values(fr)
  px <- cbind(px, fv[px$cell, , drop = FALSE])
}

# ------------------------------ CORE HELPERS ---------------------------------
urban_elev <- mean(unlist(extract(dem, ua)[, -1]), na.rm = TRUE)
urban_elev_vals <- unlist(extract(dem, ua)[, -1]); urban_elev_vals <- urban_elev_vals[is.finite(urban_elev_vals)]
log_msg("Urban core mean elevation: %.1f m (%d cells)", urban_elev, length(urban_elev_vals))

# Urban side of SUHI: city-wide mean LST and per-village mean LST
urban_mean <- global(mask(urban_lst, city), "mean", na.rm = TRUE)[1, 1]
vill_mean  <- extract(urban_lst, vill, fun = mean, na.rm = TRUE)[, 2]
log_msg("City-wide mean urban LST (sharpened): %.2f degC | villages with data: %d",
        urban_mean, sum(is.finite(vill_mean)))

base_pool <- function(r_in = 0, r_out = CFG$search_radius_km) {
  px$dist_km > 0 & !px$in_ua & px$dist_km > r_in & px$dist_km <= r_out
}

# Baseline statistic on a set of retained cells d (data.frame rows of px)
stat_matched <- function(d, tol = CFG$primary_tol_m, lo = NULL, hi = NULL) {
  lo <- if (is.null(lo)) urban_elev - tol else lo
  hi <- if (is.null(hi)) urban_elev + tol else hi
  k <- d$elev >= lo & d$elev <= hi
  list(value = if (any(k)) mean(d$lst[k]) else NA_real_, n = sum(k), idx = which(k))
}
stat_regadj <- function(d, quad = FALSE) {
  if (nrow(d) < 30) return(list(value = NA_real_, n = nrow(d)))
  m <- if (quad) lm(lst ~ elev + I(elev^2), d) else lm(lst ~ elev, d)
  list(value = as.numeric(predict(m, data.frame(elev = urban_elev))), n = nrow(d))
}
stat_histw <- function(d, bin = 25) {
  b_r <- floor(d$elev / bin); b_u <- floor(urban_elev_vals / bin)
  pu <- table(b_u) / length(b_u); pr <- table(b_r) / length(b_r)
  common <- intersect(names(pu), names(pr))
  if (!length(common)) return(list(value = NA_real_, n = 0))
  w_bin <- pu[common] / pr[common]
  w <- w_bin[as.character(b_r)]; w[is.na(w)] <- 0
  ess <- sum(w)^2 / sum(w^2)
  list(value = sum(w * d$lst) / sum(w), n = round(ess))
}

results <- list()
add <- function(family, scenario, d, st, note = "", suhi = NULL) {
  if (is.null(st$n)) st$n <- nrow(d)
  sec <- if (!is.null(st$idx)) d$sector[st$idx] else d$sector
  tb  <- table(factor(sec, levels = c("N", "E", "S", "W")))
  results[[length(results) + 1]] <<- data.frame(
    family = family, scenario = scenario, n_pixels = st$n,
    baseline_C = round(st$value, 3),
    n_N = tb[["N"]], n_E = tb[["E"]], n_S = tb[["S"]], n_W = tb[["W"]],
    max_sector_share = round(max(tb) / max(1, sum(tb)), 3),
    citywide_SUHI_C = round(if (is.null(suhi)) urban_mean - st$value else suhi, 3),
    villages_positive = if (is.null(suhi)) sum(vill_mean - st$value > 0, na.rm = TRUE) else NA_integer_,
    flag = if (!is.na(st$n) && st$n < CFG$min_pixels) "LOW_N" else "",
    note = note, stringsAsFactors = FALSE)
}

# ------------------------------ PRIMARY + CHECK ------------------------------
pool0 <- px[base_pool(), ]
pool_nonurb <- pool0[pool0$wcb %in% 1, ]
prim <- stat_matched(pool_nonurb)
log_msg("\n=== PRIMARY SCENARIO REPRODUCTION ===")
log_msg("Rural cells (UA excluded, WorldCover non-urban, <= %d km): %d", CFG$search_radius_km, nrow(pool_nonurb))
log_msg("Elevation-matched (+/-%d m) rural LST: %.2f degC (n = %d) | chapter value: %.2f",
        CFG$primary_tol_m, prim$value, prim$n, CFG$expected_baseline)
repro_ok <- is.finite(prim$value) && abs(prim$value - CFG$expected_baseline) <= CFG$baseline_tolerance
if (!repro_ok) {
  warning(sprintf(paste0("Primary scenario (%.2f) does not reproduce the chapter value (%.2f) within %.2f degC. ",
                         "Check P$lst_wide points at the same LST raster the KNIME node used before trusting ",
                         "any numbers below."), prim$value, CFG$expected_baseline, CFG$baseline_tolerance),
          call. = FALSE)
}
add("Primary", "Chapter definition (UA excl. + WC non-urban + 60 km + +/-100 m)", pool_nonurb, prim)
add("Primary", "Same, no elevation matching (full rural sample)", pool_nonurb,
    list(value = mean(pool_nonurb$lst), n = nrow(pool_nonurb)))
add("Primary", "Original municipal-boundary mask (value reported in Table 7.3)",
    pool_nonurb[0, ], list(value = 44.8, n = NA_integer_), "Reported value, not recomputed here")

# ------------------------------ A. ELEVATION ---------------------------------
for (tol in c(25, 50, 75, 100, 150, 200, 300, 500)) {
  add("A_elevation", sprintf("Matched +/-%d m", tol), pool_nonurb, stat_matched(pool_nonurb, tol))
}
add("A_elevation", "Lower half-window [urban-100, urban]", pool_nonurb,
    stat_matched(pool_nonurb, lo = urban_elev - 100, hi = urban_elev))
add("A_elevation", "Upper half-window [urban, urban+100]", pool_nonurb,
    stat_matched(pool_nonurb, lo = urban_elev, hi = urban_elev + 100))
add("A_elevation", "Regression-adjusted, linear (predict at urban elevation)", pool_nonurb, stat_regadj(pool_nonurb))
add("A_elevation", "Regression-adjusted, quadratic", pool_nonurb, stat_regadj(pool_nonurb, TRUE))
add("A_elevation", "Histogram re-weighted to urban elevation distribution (25 m bins)", pool_nonurb, stat_histw(pool_nonurb),
    "n = effective sample size (Kish)")
for (b in list(c(250, 350), c(350, 450), c(450, 550), c(550, 700))) {
  s <- stat_matched(pool_nonurb, lo = b[1], hi = b[2])
  add("A_elevation_gradient", sprintf("Absolute band %d-%d m (not matched to urban)", b[1], b[2]), pool_nonurb, s,
      "Shows the LST-elevation gradient, not a candidate baseline")
}

# ------------------------------ B. LAND COVER --------------------------------
lc_rule <- function(pool, classes, min_frac) {
  cols <- paste0("fr_", classes)
  if (!all(cols %in% names(pool))) return(rep(FALSE, nrow(pool)))
  rowSums(pool[, cols, drop = FALSE], na.rm = TRUE) >= min_frac
}
pool_all <- px[base_pool(), ]          # UA excluded, 60 km, no land-cover filter yet
add("B_landcover", "Existing aligned WorldCover mask (chapter)", pool_nonurb, prim)
if (have_frac) {
  sets <- list(
    "Shrub+grass+crop+bare, >=50% of cell (reconstruction of chapter mask)" = list(c(20, 30, 40, 60), 0.5),
    "Shrub+grass+crop+bare, >=80% of cell (stricter purity)"               = list(c(20, 30, 40, 60), 0.8),
    "Excluding cropland: shrub+grass+bare, >=50%"                          = list(c(20, 30, 60), 0.5),
    "Native desert only: shrub+grass, >=50%"                               = list(c(20, 30), 0.5),
    "Bare/sparse only, >=50%"                                              = list(c(60), 0.5),
    "Cropland only, >=50%"                                                 = list(c(40), 0.5),
    "All non-built, non-water classes, >=50%"                              = list(c(10, 20, 30, 40, 60, 90), 0.5)
  )
  for (nm in names(sets)) {
    k <- lc_rule(pool_all, sets[[nm]][[1]], sets[[nm]][[2]])
    d <- pool_all[k, ]
    add("B_landcover", nm, d, stat_matched(d))
  }
  for (bt in c(0.02, 0.05, 0.10, 0.20)) {
    d <- pool_all[pool_all$fr_50 <= bt & pool_all$fr_80 <= 0.05, ]
    add("B_landcover", sprintf("Any land cover with built-up fraction <= %d%%", round(bt * 100)), d, stat_matched(d))
  }
  agree <- mean(lc_rule(pool_all, c(20, 30, 40, 60), 0.5) == (pool_all$wcb %in% 1))
  log_msg("Agreement between reconstructed and existing WorldCover mask: %.1f%% of cells", 100 * agree)
}

# ------------------------------ C. BUFFERS -----------------------------------
for (r_in in c(0, 2.5, 5, 10, 15, 20)) {
  for (r_out in c(20, 30, 40, 60, 80)) {
    if (r_out <= r_in + 5) next
    d <- pool_nonurb[pool_nonurb$dist_km > r_in & pool_nonurb$dist_km <= r_out, ]
    add("C_buffer", sprintf("Inner exclusion %g km, outer radius %g km", r_in, r_out), d, stat_matched(d))
  }
}

# ------------------------------ D. SAMPLE & DISTRIBUTION ---------------------
sel <- prim$idx; dm <- pool_nonurb[sel, ]
secs <- c("N", "E", "S", "W")
for (s in secs) {
  d <- dm[dm$sector != s, ]
  add("D_distribution", paste("Leave out", s, "sector"), d, list(value = mean(d$lst), n = nrow(d)))
}
sec_means <- tapply(dm$lst, factor(dm$sector, levels = secs), mean)
add("D_distribution", "Sector-balanced mean (equal weight N/E/S/W)", dm,
    list(value = mean(sec_means, na.rm = TRUE), n = nrow(dm)), "Mean of the four sector means")
# Block bootstrap
bsum <- tapply(dm$lst, dm$block, sum); bn <- tapply(dm$lst, dm$block, length)
boot <- replicate(CFG$n_boot, { i <- sample(seq_along(bn), replace = TRUE); sum(bsum[i]) / sum(bn[i]) })
ci <- quantile(boot, c(0.025, 0.975))
add("D_distribution", sprintf("Spatial block bootstrap (%d km blocks, %d reps): median", CFG$block_km, CFG$n_boot),
    dm, list(value = median(boot), n = nrow(dm)),
    sprintf("95%% CI %.2f to %.2f degC; %d blocks", ci[1], ci[2], length(bn)))
# Random subsampling
for (frac in c(0.10, 0.25, 0.50)) {
  m <- replicate(500, mean(sample(dm$lst, max(2, round(frac * nrow(dm))))))
  add("D_distribution", sprintf("Random %d%% subsample (500 reps): mean", round(frac * 100)),
      dm, list(value = mean(m), n = round(frac * nrow(dm))),
      sprintf("SD across reps %.3f degC", sd(m)))
}

# ------------------------------ E. TEMPORAL ----------------------------------
scene_files <- list.files(P$lst_dir, pattern = P$lst_scene_pattern, full.names = TRUE)
temporal <- NULL
if (length(scene_files) >= 2) {
  doy <- as.integer(sub(".*doy2024(\\d{3}).*", "\\1", basename(scene_files)))
  ord <- order(doy); scene_files <- scene_files[ord]; doy <- doy[ord]
  scenes <- lapply(scene_files, function(f) align_to_lst(rast(f)[[1]], "bilinear"))
  date <- as.Date(doy - 1, origin = "2024-01-01")
  mask_idx <- pool_nonurb$cell[sel]
  rows <- lapply(seq_along(scenes), function(i) {
    v <- values(scenes[[i]])[, 1]
    rl <- v[mask_idx]; rl <- rl[is.finite(rl)]
    um <- global(mask(scenes[[i]], city), "mean", na.rm = TRUE)[1, 1]
    data.frame(date = date[i], doy = doy[i], month = format(date[i], "%Y-%m"),
               n_rural = length(rl), rural_C = if (length(rl)) mean(rl) else NA_real_, urban_C = um)
  })
  temporal <- do.call(rbind, rows)
  temporal$SUHI_C <- temporal$urban_C - temporal$rural_C
  write.csv(temporal, file.path(P$out_dir, "sensitivity_temporal_by_scene.csv"), row.names = FALSE)
  ok <- temporal[temporal$n_rural >= CFG$min_pixels, ]
  for (m in unique(ok$month)) {
    s <- ok[ok$month == m, ]
    add("E_temporal", paste("Month", m, sprintf("(%d scenes, mean of per-scene baselines)", nrow(s))),
        dm[0, ], list(value = mean(s$rural_C), n = round(mean(s$n_rural))),
        "Seasonal baseline: compare SUHI (urban - rural of the SAME scenes), not the absolute baseline",
        suhi = mean(s$SUHI_C))
  }
  add("E_temporal", sprintf("All %d scenes, mean of per-scene baselines", nrow(ok)), dm[0, ],
      list(value = mean(ok$rural_C), n = round(mean(ok$n_rural))),
      sprintf("per-scene baseline range %.2f-%.2f degC", min(ok$rural_C), max(ok$rural_C)),
      suhi = mean(ok$SUHI_C))
  stk <- rast(scenes)
  for (fn in c("median", "mean")) {
    comp <- app(stk, fn, na.rm = TRUE)
    v <- values(comp)[, 1][mask_idx]; v <- v[is.finite(v)]
    add("E_temporal", sprintf("Re-built %s composite of the scenes", fn), dm[0, ],
        list(value = mean(v), n = length(v)),
        "Elevation-matched primary mask; composite recomputed from the per-scene files")
  }
} else {
  log_msg("NOTE: fewer than two per-scene LST files matched P$lst_scene_pattern -- family E skipped.")
}

# ------------------------------ ASSEMBLE & DECIDE ----------------------------
res <- do.call(rbind, results)
write.csv(res, file.path(P$out_dir, "sensitivity_summary.csv"), row.names = FALSE)

defensible <- res[res$family %in% c("A_elevation", "B_landcover", "C_buffer", "D_distribution") &
                    res$flag == "" & is.finite(res$baseline_C) &
                    !grepl("Cropland only|Bare/sparse only|Absolute band|Random|Leave out|Existing|Lower half|Upper half|Matched \\+/-(300|500) m", res$scenario), ]
rng  <- range(defensible$baseline_C)
sprd <- diff(rng)
by_fam <- aggregate(baseline_C ~ family, defensible, function(x) round(diff(range(x)), 2))
prim_v <- prim$value
dec <- c(
  "SENSITIVITY OF THE RURAL REFERENCE-ZONE BASELINE (Chapter 7)",
  sprintf("Primary (chapter) baseline, reproduced here : %.2f degC (n = %d) [chapter: %.2f]", prim_v, prim$n, CFG$expected_baseline),
  sprintf("Reproduction within %.2f degC                : %s", CFG$baseline_tolerance, if (repro_ok) "YES" else "NO -- check P$lst_wide"),
  sprintf("Scenarios retained as 'defensible' (n >= %d)  : %d", CFG$min_pixels, nrow(defensible)),
  sprintf("Range of baseline across them               : %.2f to %.2f degC (spread %.2f degC)", rng[1], rng[2], sprd),
  sprintf("Spread by family (degC)                     : %s", paste(by_fam$family, by_fam$baseline_C, sep = "=", collapse = "; ")),
  sprintf("Bootstrap 95%% CI of primary baseline        : %.2f to %.2f degC", ci[1], ci[2]),
  sprintf("City-wide mean urban LST used for SUHI      : %.2f degC", urban_mean),
  sprintf("Implied city-wide SUHI across scenarios     : %.2f to %.2f degC (primary %.2f)",
          urban_mean - rng[2], urban_mean - rng[1], urban_mean - prim_v),
  sprintf("Sign of city-wide SUHI stable?              : %s",
          if (sign(urban_mean - rng[2]) == sign(urban_mean - rng[1])) "YES" else "NO -- sign flips inside the range"),
  if (!is.null(temporal)) sprintf("Temporal (family E, kept separate): per-scene city-wide SUHI mean %.2f, range %.2f to %.2f degC across %d scenes; seasonal change in baseline is mirrored by the urban side, so SUHI is the quantity to compare",
          mean(temporal$SUHI_C, na.rm = TRUE), min(temporal$SUHI_C, na.rm = TRUE), max(temporal$SUHI_C, na.rm = TRUE), nrow(temporal)) else "Temporal (family E): skipped",
  sprintf("Villages with SUHI > 0 across scenarios     : %d to %d of %d",
          min(defensible$villages_positive), max(defensible$villages_positive), sum(is.finite(vill_mean))),
  "",
  if (sprd <= CFG$robust_threshold_c) {
    sprintf("VERDICT: ROBUST. The baseline moves by <= %.1f degC across all defensible alternatives. The sensitivity table can be reported as supporting evidence.", CFG$robust_threshold_c)
  } else {
    sprintf("VERDICT: SENSITIVE. The baseline moves by %.2f degC (> %.1f degC) across defensible alternatives. Report the range as an explicit uncertainty on the SUHI magnitude and name the families driving it (see 'Spread by family').", sprd, CFG$robust_threshold_c)
  },
  "",
  "NOTE: the baseline is a single constant subtracted from every pixel, so these choices shift the MAGNITUDE (and possibly the sign) of SUHI",
  "everywhere by the same amount; they do not change the ranking of villages or pixels.",
  "Scenarios flagged LOW_N, single-class land-cover tests, absolute elevation bands, matching tolerances of +/-300 and +/-500 m (which no longer match the urban core),",
  "leave-one-out and random subsamples are in sensitivity_summary.csv but excluded from the range above because they are diagnostics, not candidate definitions."
)
writeLines(dec, file.path(P$out_dir, "sensitivity_decision.txt"))
cat("\n", paste(dec, collapse = "\n"), "\n", sep = "")

# ------------------------------ PLOTS ----------------------------------------
fam_col <- c(Primary = "black", A_elevation = "#1b9e77", A_elevation_gradient = "#bbbbbb", B_landcover = "#d95f02",
             C_buffer = "#7570b3", D_distribution = "#e7298a", E_temporal = "#66a61e")
pl <- res[is.finite(res$baseline_C) & res$family != "A_elevation_gradient" & res$flag == "", ]
pl <- pl[nrow(pl):1, ]
png(file.path(P$out_dir, "sensitivity_forest_plot.png"), width = 2000, height = max(1400, 28 * nrow(pl) + 300), res = 200)
par(mar = c(4.5, 2, 2.5, 1), mfrow = c(1, 1))
plot(pl$baseline_C, seq_len(nrow(pl)), pch = 16, col = fam_col[pl$family], yaxt = "n", xlab = "Rural baseline LST (degC)",
     ylab = "", main = "Rural reference-zone baseline under alternative definitions",
     xlim = range(c(pl$baseline_C, prim_v - 0.5, prim_v + 0.5)))
abline(v = prim_v, lty = 2); rect(prim_v - 0.5, 0, prim_v + 0.5, nrow(pl) + 1, col = rgb(0, 0, 0, 0.06), border = NA)
text(par("usr")[1], seq_len(nrow(pl)), substr(pl$scenario, 1, 70), pos = 4, cex = 0.45, xpd = TRUE)
legend("bottomright", legend = names(fam_col)[names(fam_col) %in% pl$family], col = fam_col[names(fam_col) %in% pl$family],
       pch = 16, cex = 0.6, bty = "n")
dev.off()

if (!is.null(temporal)) {
  png(file.path(P$out_dir, "sensitivity_temporal.png"), width = 2000, height = 1000, res = 200)
  par(mar = c(4.5, 4.5, 2.5, 4.5))
  plot(temporal$date, temporal$rural_C, type = "b", pch = 16, col = "#d95f02", ylab = "Rural baseline LST (degC)",
       xlab = "Scene date (8-day MOD11A2)", main = "Per-scene rural baseline and city-wide SUHI")
  abline(h = prim_v, lty = 2)
  par(new = TRUE)
  plot(temporal$date, temporal$SUHI_C, type = "b", pch = 17, col = "#1b9e77", axes = FALSE, xlab = "", ylab = "")
  axis(4); mtext("City-wide SUHI (degC)", side = 4, line = 3)
  legend("topleft", c("Rural baseline (left)", "City-wide SUHI (right)", "Composite primary baseline"),
         col = c("#d95f02", "#1b9e77", "black"), pch = c(16, 17, NA), lty = c(1, 1, 2), cex = 0.7, bty = "n")
  dev.off()
}
capture.output(sessionInfo(), file = file.path(P$out_dir, "sessionInfo.txt"))
log_msg("\nDone. Results in: %s", P$out_dir)

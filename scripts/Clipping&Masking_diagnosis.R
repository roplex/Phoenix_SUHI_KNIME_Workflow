library(terra)

# ---- Inputs from KNIME ----
lst_clipped_path  <- as.character(knime.in[["lst_clipped_path"]])   # Already clipped LST
ndvi_clipped_path <- as.character(knime.in[["ndvi_clipped_path"]])  # Already clipped NDVI

# ---- Safety checks ----
if (!file.exists(lst_clipped_path))  stop(paste("Clipped LST not found:", lst_clipped_path))
if (!file.exists(ndvi_clipped_path)) stop(paste("Clipped NDVI not found:", ndvi_clipped_path))

# ---- Load clipped rasters ----
lst_clipped  <- rast(lst_clipped_path)
ndvi_clipped <- rast(ndvi_clipped_path)

# ---- Diagnostics function ----
get_stats <- function(r) {
  v <- values(r, na.rm=TRUE)
  n_total <- ncell(r)
  n_valid <- length(v)
  data.frame(
    valid_pixels = n_valid,
    total_pixels = n_total,
    valid_pct    = round(100 * n_valid / n_total, 2),
    min  = ifelse(n_valid > 0, min(v), NA),
    max  = ifelse(n_valid > 0, max(v), NA),
    mean = ifelse(n_valid > 0, mean(v), NA),
    sd   = ifelse(n_valid > 0, sd(v), NA),
    p05  = ifelse(n_valid > 0, quantile(v, 0.05), NA),
    p10  = ifelse(n_valid > 0, quantile(v, 0.10), NA),
    p25  = ifelse(n_valid > 0, quantile(v, 0.25), NA),
    p50  = ifelse(n_valid > 0, quantile(v, 0.50), NA),
    p75  = ifelse(n_valid > 0, quantile(v, 0.75), NA),
    p90  = ifelse(n_valid > 0, quantile(v, 0.90), NA),
    p95  = ifelse(n_valid > 0, quantile(v, 0.95), NA)
  )
}

# ---- Run diagnostics ----
lst_stats  <- get_stats(lst_clipped)
ndvi_stats <- get_stats(ndvi_clipped)

# ---- Optionally: save diagnostics to CSVs ----
out_dir <- dirname(lst_clipped_path)
lst_csv  <- file.path(out_dir, "LST_clipped_diagnostics.csv")
ndvi_csv <- file.path(out_dir, "NDVI_clipped_diagnostics.csv")

write.csv(lst_stats,  lst_csv,  row.names=FALSE)
write.csv(ndvi_stats, ndvi_csv, row.names=FALSE)

# ---- Return outputs to KNIME ----
knime.out <- list(
  lst_diagnostics  = lst_stats,
  ndvi_diagnostics = ndvi_stats,
  csv_paths = data.frame(
    lst_csv  = lst_csv,
    ndvi_csv = ndvi_csv
  )
)

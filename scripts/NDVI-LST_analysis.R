# R Snippet: ndvi_lst_analysis_export_flat_fixed.R
library(terra)

# --- Inputs from KNIME ---
lst_path  <- as.character(knime.in[["lst_path"]])
ndvi_path <- as.character(knime.in[["ndvi_path"]])
suhi_path <- as.character(knime.in[["suhi_path"]])

# --- Load rasters ---
lst  <- rast(lst_path)[[1]]
ndvi <- rast(ndvi_path)[[1]]

# Align CRS
if (!identical(crs(lst), crs(ndvi))) {
  ndvi <- project(ndvi, crs(lst))
}

# Align extent/resolution
if (!all.equal(ext(lst), ext(ndvi)) || ncell(lst) != ncell(ndvi)) {
  ndvi <- resample(ndvi, lst, method="bilinear")
}

# --- Extract paired values ---
v_lst  <- as.numeric(values(lst))
v_ndvi <- as.numeric(values(ndvi))

df <- data.frame(LST = v_lst, NDVI = v_ndvi)
df <- df[complete.cases(df), ]

# Guard: stop if no valid data
if (nrow(df) < 5) {
  stop("Too few valid pixel pairs for regression")
}

# --- Regression + correlation ---
fit <- lm(LST ~ NDVI, data=df)
r2  <- summary(fit)$r.squared
pearson_r <- suppressWarnings(cor(df$NDVI, df$LST))

# --- Sample data for scatterplot ---
sample_size <- min(5000, nrow(df))
set.seed(42)
df_sample <- df[sample(nrow(df), sample_size), ]

# --- Replicate stats across rows of df_sample ---
n <- nrow(df_sample)
knime.out <- data.frame(
  slope     = rep(as.numeric(coef(fit)[2]), n),
  intercept = rep(as.numeric(coef(fit)[1]), n),
  r_squared = rep(as.numeric(r2), n),
  pearson_r = rep(as.numeric(pearson_r), n),
  n_points  = rep(as.integer(nrow(df)), n),
  suhi_path = rep(suhi_path, n),
  lst_path  = rep(lst_path, n),
  ndvi_path = rep(ndvi_path, n),
  NDVI      = df_sample$NDVI,
  LST       = df_sample$LST
)

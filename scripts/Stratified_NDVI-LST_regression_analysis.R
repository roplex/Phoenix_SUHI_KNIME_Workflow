# --- Libraries ---
library(terra)
library(dplyr)

# --- Inputs from KNIME ---
lst_path  <- as.character(knime.in[["lst_path"]])
ndvi_path <- as.character(knime.in[["ndvi_path"]])
suhi_path <- as.character(knime.in[["suhi_path"]])

# --- Fixed village path ---
villages_path <- "/Users/roplex/Desktop/EO_Harmonization/PhoenixData/Boundary/Villages.geojson"

# --- Load rasters and vector ---
lst  <- rast(lst_path)[[1]]
ndvi <- rast(ndvi_path)[[1]]
villages <- vect(villages_path)

# Align CRS between rasters and villages
if (!identical(crs(lst), crs(ndvi))) {
  ndvi <- project(ndvi, crs(lst))
}
if (!identical(crs(lst), crs(villages))) {
  villages <- project(villages, crs(lst))
}

# Align extent/resolution
if (!all.equal(ext(lst), ext(ndvi)) || ncell(lst) != ncell(ndvi)) {
  ndvi <- resample(ndvi, lst, method = "bilinear")
}

# --- Stack rasters for extraction ---
stacked <- c(lst, ndvi)
names(stacked) <- c("LST", "NDVI")

# --- Extract data and run regressions by village ---
village_names <- villages$NAME
results_list <- list()

for (i in seq_along(village_names)) {
  v <- villages[i]
  df_v <- terra::extract(stacked, v, na.rm = TRUE)
  
  # Handle cases with missing or too few pixels
  if (is.null(df_v) || nrow(df_v) < 5) {
    message(paste("Too few valid pixels for:", village_names[i], "-> Filling with NA"))
    res <- data.frame(
      village = village_names[i],
      slope = NA,
      intercept = NA,
      r_squared = NA,
      pearson_r = NA,
      n_points = ifelse(is.null(df_v), 0, nrow(df_v)),
      suhi_path = suhi_path,
      lst_path = lst_path,
      ndvi_path = ndvi_path
    )
    results_list[[length(results_list) + 1]] <- res
    next
  }
  
  # Clean and model
  df_v <- df_v[complete.cases(df_v), ]
  if (nrow(df_v) < 5) next
  
  fit <- lm(LST ~ NDVI, data = df_v)
  r2 <- summary(fit)$r.squared
  pearson_r <- suppressWarnings(cor(df_v$NDVI, df_v$LST))
  
  # Store valid results
  res <- data.frame(
    village = village_names[i],
    slope = coef(fit)[2],
    intercept = coef(fit)[1],
    r_squared = r2,
    pearson_r = pearson_r,
    n_points = nrow(df_v),
    suhi_path = suhi_path,
    lst_path = lst_path,
    ndvi_path = ndvi_path
  )
  
  results_list[[length(results_list) + 1]] <- res
}

# --- Combine all results safely for KNIME output ---
knime.out <- as.data.frame(bind_rows(results_list))

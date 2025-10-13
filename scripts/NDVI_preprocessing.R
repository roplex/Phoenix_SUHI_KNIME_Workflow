# R Snippet for KNIME: NDVI QA decoding + masks (works with SpatRaster bit ops)
library(terra)

# --- Get KNIME input (robust to knime.flow.in vs knime.in) ---
kn <- if (exists("knime.flow.in")) knime.flow.in else if (exists("knime.in")) knime.in else stop("KNIME input object not found")
ndvi_path <- as.character(kn[["ndvi_path"]])
if (is.na(ndvi_path) || ndvi_path == "") stop("ndvi_path not provided in KNIME input")

# --- Build QC path (adjust pattern if your filenames differ) ---
qc_path <- gsub("NDVI", "VI_Quality", ndvi_path)
stopifnot(file.exists(ndvi_path), file.exists(qc_path))

# --- Load rasters ---
ndvi <- rast(ndvi_path)
qc   <- rast(qc_path)

# Debug prints (optional)
print(paste("NDVI raster:", ndvi_path))
print(paste("QC raster :", qc_path))
print(paste("NDVI dims:", paste(dim(ndvi), collapse = " x ")))
print(paste("QC  dims:", paste(dim(qc), collapse = " x ")))

# --- Decode QC bits (use app() so results stay SpatRaster) ---
# Bits per MOD13: bits 0-1 = MODLAND QA, bits 2-5 = VI usefulness (shifted right by 2)
modland_qc <- app(qc, fun = function(x) bitwAnd(as.integer(x), 3))                 # 0..3
vi_use_r   <- app(qc, fun = function(x) floor(bitwAnd(as.integer(x), 60) / 4))    # 0..15

# Quick check (print a small frequency sample)
print("Sample unique QC values (qc):")
print(sort(unique(head(values(qc, mat=FALSE), 50))))
print("Sample decoded modland_qc unique (first 30):")
print(sort(unique(head(values(modland_qc, mat=FALSE), 30))))
print("Sample decoded vi_use unique (first 30):")
print(sort(unique(head(values(vi_use_r, mat=FALSE), 30))))

# --- Pre-scale NDVI (we will keep scaled values where QA passes) ---
ndvi_scaled <- ndvi * 0.0001

# --- Masks using ifel so alignment/type preserved ---
# Strict: modland_qc == 0 and vi_use <= 1
ndvi_strict  <- ifel((modland_qc == 0) & (vi_use_r <= 1), ndvi_scaled, NA)

# Relaxed: modland_qc <= 1 and vi_use <= 2
ndvi_relaxed <- ifel((modland_qc <= 1) & (vi_use_r <= 2), ndvi_scaled, NA)

# Super-relaxed: vi_use <= 3
ndvi_super   <- ifel(vi_use_r <= 3, ndvi_scaled, NA)

# --- Save outputs ---
out_dir <- dirname(ndvi_path)
out_strict  <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(ndvi_path)), "_clean_strict.tif"))
out_relaxed <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(ndvi_path)), "_clean_relaxed.tif"))
out_super   <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(ndvi_path)), "_clean_superrelaxed.tif"))

writeRaster(ndvi_strict,  out_strict,  overwrite=TRUE)
writeRaster(ndvi_relaxed, out_relaxed, overwrite=TRUE)
writeRaster(ndvi_super,   out_super,   overwrite=TRUE)

# --- Coverage stats (percent) ---
tot_pix <- ncell(ndvi)
cov_strict  <- 100 * sum(!is.na(values(ndvi_strict)))  / tot_pix
cov_relaxed <- 100 * sum(!is.na(values(ndvi_relaxed))) / tot_pix
cov_super   <- 100 * sum(!is.na(values(ndvi_super)))   / tot_pix

# Small diagnostics (min/max/mean on strict mask)
safe_min <- function(r) tryCatch(global(r, "min", na.rm=TRUE)[[1]], error = function(e) NA)
safe_max <- function(r) tryCatch(global(r, "max", na.rm=TRUE)[[1]], error = function(e) NA)
safe_mean<- function(r) tryCatch(global(r, "mean", na.rm=TRUE)[[1]], error = function(e) NA)

diag_strict <- c(min = safe_min(ndvi_strict), max = safe_max(ndvi_strict), mean = safe_mean(ndvi_strict))
print("NDVI strict diagnostics (min,max,mean):")
print(diag_strict)

# --- Return outputs to KNIME ---
knime.out <- data.frame(
  strict_path  = out_strict,
  relaxed_path = out_relaxed,
  super_path   = out_super,
  strict_cov   = cov_strict,
  relaxed_cov  = cov_relaxed,
  super_cov    = cov_super
)

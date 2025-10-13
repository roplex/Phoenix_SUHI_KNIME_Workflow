# R Snippet: compute_suhi_urban_only.R
library(terra)

lst_path  <- as.character(knime.in[["lst_path"]])
ndvi_path  <- as.character(knime.in[["ndvi_path"]])
mask_path <- as.character(knime.in[["urban_mask_path"]])

lst  <- rast(lst_path)[[1]]
mask <- rast(mask_path)[[1]]

# ensure alignment
if (!identical(crs(lst), crs(mask))) mask <- project(mask, crs(lst))
if (!all.equal(ext(lst), ext(mask)) || ncell(lst) != ncell(mask)) {
  mask <- resample(mask, lst, method="near")
}

# rural mean (mask == 0)
mean_rural <- global(mask(lst, mask, maskvalue=1), "mean", na.rm=TRUE)[[1]]

if (is.na(mean_rural)) stop("No valid rural LST pixels found (mean_rural is NA)")

# SUHI = ΔT but only for urban pixels (mask == 1)
suhi <- ifel(mask == 1, lst - mean_rural, NA)

out_suhi <- file.path(dirname(lst_path), "SUHI_map_urbanOnly.tif")
writeRaster(suhi, out_suhi, overwrite=TRUE)

# quick histogram for SUHI values (urban only)
vals <- values(suhi, na.rm=TRUE)
hist_counts <- hist(vals, breaks=50, plot=FALSE)
hist_df <- data.frame(
  break_left  = hist_counts$breaks[-length(hist_counts$breaks)],
  break_right = hist_counts$breaks[-1],
  count       = hist_counts$counts
)
write.csv(hist_df, file.path(dirname(lst_path), "suhi_histogram_urbanOnly.csv"), row.names=FALSE)

knime.out <- data.frame(
  suhi_path       = out_suhi,
  mean_rural_LST  = mean_rural,
  histogram = hist_df,
  lst_path = lst_path,
  ndvi_path = ndvi_path
)

# Run inside KNIME R Snippet (Table) node
library(terra)

# Input path from KNIME (force character)
lst_path <- as.character(knime.in[["resampled_lst_path"]])  # e.g. LST_composite_superrelaxed_resampledToNDVI.tif

# Load raster
r <- rast(lst_path)

# Inspect metadata
print(crs(r))
print(ext(r))
print(dim(r))

# Extract raw values
vals_all <- values(r, mat=FALSE)

# Frequency table of most common values (can help catch fill values)
freq <- sort(table(vals_all), decreasing=TRUE)
print("Top value frequencies:")
print(head(freq, 20))

# Handle sentinel nodata values if present
possible_sentinels <- c(-9999, -32768, 32767, 9999)
present_sentinels <- intersect(as.numeric(names(freq)), possible_sentinels)
if (length(present_sentinels) > 0) {
  cat("Detected sentinel nodata values:", paste(present_sentinels, collapse=", "), "\n")
  for (s in present_sentinels) {
    r[r == s] <- NA
  }
}

# Compute stats
r_min  <- global(r, "min", na.rm=TRUE)[[1]]
r_max  <- global(r, "max", na.rm=TRUE)[[1]]
r_mean <- global(r, "mean", na.rm=TRUE)[[1]]
r_sd   <- global(r, "sd", na.rm=TRUE)[[1]]
valid_count <- sum(!is.na(values(r)))
total_count <- ncell(r)

cat(sprintf("Valid pixels: %d / %d (%.2f%%)\n", valid_count, total_count, 
            100*valid_count/total_count))
cat(sprintf("LST min=%.2f  max=%.2f  mean=%.2f  sd=%.2f\n", 
            r_min, r_max, r_mean, r_sd))

# Percentiles
vals <- values(r, mat=FALSE)
p <- quantile(vals, 
              probs = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1), 
              na.rm=TRUE)
print("Percentiles:")
print(round(p, 2))

# Histogram (10 bins)
h <- hist(vals, breaks=10, plot=FALSE)
hist_df <- data.frame(
  break_left  = h$breaks[-length(h$breaks)],
  break_right = h$breaks[-1],
  count       = h$counts
)
print("Histogram (10 bins):")
print(hist_df)

# Save diagnostics to CSV
out_csv <- file.path(dirname(lst_path), "lst_resampled_diagnostics.csv")
diag_df <- data.frame(
  stat  = c("valid_pixels","total_pixels","valid_pct","min","max","mean","sd",
            "p05","p10","p25","p50","p75","p90","p95"),
  value = c(valid_count, total_count, 100*valid_count/total_count,
            r_min, r_max, r_mean, r_sd,
            p["5%"], p["10%"], p["25%"], p["50%"], p["75%"], p["90%"], p["95%"])
)
write.csv(diag_df, out_csv, row.names=FALSE)

# Return to KNIME
knime.out <- data.frame(
  lst_diag_csv = out_csv,
  lst_min = r_min,
  lst_max = r_max,
  lst_mean = r_mean,
  lst_sd = r_sd,
  lst_p10 = as.numeric(p["10%"]),
  lst_p25 = as.numeric(p["25%"])
)

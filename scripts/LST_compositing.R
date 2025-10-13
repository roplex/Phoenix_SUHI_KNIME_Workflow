library(terra)

# Collect file paths from KNIME input (force character)
files <- as.character(knime.in[["superrelaxed_lst_path"]])

# Load rasters into a stack
rlist <- lapply(files, rast)
stack <- rast(rlist)

# Median composite (robust to outliers)
lst_composite <- app(stack, median, na.rm=TRUE)

# Save composite
out_path <- file.path(dirname(files[1]), "LST_composite_superrelaxed.tif")
writeRaster(lst_composite, out_path, overwrite=TRUE)

# Return to KNIME
knime.out <- data.frame(lst_composite_path = out_path)

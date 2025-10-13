library(terra)

files <- as.character(knime.in[["super_path"]])

rlist <- lapply(files, rast)
stack <- rast(rlist)

ndvi_composite <- app(stack, median, na.rm=TRUE)

out_path <- file.path(dirname(files[1]), "NDVI_composite_super.tif")
writeRaster(ndvi_composite, out_path, overwrite=TRUE)

knime.out <- data.frame(ndvi_composite_path = out_path)

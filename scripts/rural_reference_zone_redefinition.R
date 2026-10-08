# ============================================================================
# Redefining the SUHI rural reference zone using genuine non-urban land,
# replacing the original "inside the bounding box, outside the city line"
# definition.
#
# WHY THIS CHANGES:
# The original "rural" mask kept any pixel inside the Phoenix city bounding
# box but outside the Phoenix municipal line. Phoenix sits inside a
# contiguous built-up metropolitan area (Scottsdale, Tempe, Mesa, Chandler,
# Glendale, Peoria, and others) with no natural break at the municipal
# boundary, so that mask was mostly still measuring urban or suburban
# surface temperature - contaminating the SUHI baseline with urban heat and
# very likely causing the reported SUHI intensity to be UNDERESTIMATED.
#
# This script replaces the municipal-boundary approach with a two-part
# non-urban mask:
#   1. Outside the US Census Bureau's Phoenix-Mesa-Scottsdale Urbanized Area
#      boundary - a standard, externally defined, defensible measure of the
#      contiguous built-up extent, avoiding an ad hoc list of city names.
#   2. Classified as non-built-up land cover in ESA WorldCover, so isolated
#      exurban development within the search radius is also excluded.
#
# Requires: install.packages(c("terra", "tigris"))
# WorldCover: download the v200 tile(s) covering the Phoenix region from the
# Copernicus Data Space or the WorldCover S3 bucket (this workflow is
# GEE-free, consistent with the rest of the chapter's tooling).
# ============================================================================

library(terra)
library(tigris)   # US Census TIGER/Line boundaries

#' Build a genuine non-urban reference mask for SUHI baseline computation
#'
#' @param lst_grid SpatRaster defining the target analysis grid (e.g., the
#'   TsHARP-sharpened 250 m LST composite) - the mask is aligned to this.
#' @param worldcover_path Local path to the ESA WorldCover v200 tile(s)
#'   covering the search radius around Phoenix.
#' @param search_radius_km How far out from the urbanized-area boundary to
#'   search for valid non-urban reference pixels. 60 km reaches genuine
#'   Sonoran Desert and irrigated agricultural land west and south of the
#'   metro area without crossing into a climatically distinct zone.
#' @return SpatRaster mask (1 = valid rural reference pixel, NA otherwise)
build_rural_reference_mask <- function(lst_grid, worldcover_path, search_radius_km = 60) {

  # Step 1: Census Urbanized Area boundary - the standard alternative to an
  # ad hoc list of municipality names, and what the reviewer's critique
  # effectively asks for: a boundary that actually reflects contiguous
  # built-up extent rather than an administrative line.
  ua <- urban_areas(year = 2020)
  phx_ua <- ua[grepl("Phoenix--Mesa--Scottsdale", ua$NAME20), ]
  if (nrow(phx_ua) == 0) {
    stop("Phoenix-Mesa-Scottsdale Urbanized Area not found - check the ",
         "NAME20 field against the current TIGER/Line vintage.")
  }
  phx_ua <- vect(phx_ua) |> project(crs(lst_grid))

  # Step 2: search extent - a buffer around the urbanized area, so the rural
  # reference is drawn from the surrounding desert/agricultural land rather
  # than from an arbitrarily distant, climatically different region.
  search_extent <- buffer(phx_ua, width = search_radius_km * 1000)

  # Step 3: land-cover mask - keep only non-built-up classes.
  # ESA WorldCover v200 codes: 10 Tree cover, 20 Shrubland, 30 Grassland,
  # 40 Cropland, 50 Built-up (EXCLUDED), 60 Bare/sparse vegetation,
  # 70 Snow/ice, 80 Water, 90 Wetland, 95 Mangroves, 100 Moss/lichen.
  # For the Sonoran Desert setting, the realistic rural classes are
  # Shrubland, Grassland, Cropland, and Bare/sparse vegetation.
  worldcover <- rast(worldcover_path)
  worldcover <- crop(worldcover, project(search_extent, crs(worldcover)))
  worldcover <- project(worldcover, crs(lst_grid), method = "near")

  non_urban_classes <- c(20, 30, 40, 60)
  non_urban_mask <- classify(
    worldcover,
    cbind(non_urban_classes, rep(1, length(non_urban_classes))),
    others = NA
  )
  non_urban_mask <- resample(non_urban_mask, lst_grid, method = "near")

  # Step 4: combine - rural = non-built-up land cover AND outside the
  # urbanized area boundary (belt-and-braces against isolated exurban
  # development just past the UA line).
  rural_mask <- mask(non_urban_mask, project(phx_ua, crs(lst_grid)), inverse = TRUE)
  names(rural_mask) <- "rural_reference_mask"
  rural_mask
}

#' Quantify the contamination in the original definition by comparing rural
#' mean LST under both masks. Report this in the chapter's Quality
#' Assessment / Uncertainty section - it directly substantiates (and
#' quantifies the direction of) the reviewer's underestimation claim.
compare_rural_definitions <- function(lst_grid, old_rural_mask, new_rural_mask) {
  old_mean <- global(mask(lst_grid, old_rural_mask), "mean", na.rm = TRUE)[1, 1]
  new_mean <- global(mask(lst_grid, new_rural_mask), "mean", na.rm = TRUE)[1, 1]
  cat(sprintf("Old (bounding-box) rural mean LST:      %.2f degC\n", old_mean))
  cat(sprintf("New (land-cover-based) rural mean LST:  %.2f degC\n", new_mean))
  cat(sprintf(
    "Difference: %.2f degC (positive => old mask was warmer, i.e., contaminated; ",
    old_mean - new_mean
  ))
  cat("every village's SUHI value will shift upward by roughly this amount.)\n")
  invisible(list(old_mean = old_mean, new_mean = new_mean, delta = old_mean - new_mean))
}

# ============================================================================
# USAGE - drop into the KNIME R Snippet node currently doing the urban/rural
# bounding-box mask (Section 7.5, "Urban vs. Rural Masking")
# ============================================================================
#
# lst_grid <- rast(knime.in[["lst_sharpened_path"]])   # output of TsHARP step
#
# new_rural_mask <- build_rural_reference_mask(
#   lst_grid,
#   worldcover_path = "worldcover_phoenix_region.tif",
#   search_radius_km = 60
# )
#
# # Recompute the SUHI baseline with the corrected mask:
# rural_lst_mean <- global(mask(lst_grid, new_rural_mask), "mean", na.rm = TRUE)[1, 1]
#
# # Recommended: report the before/after comparison directly in the chapter -
# old_rural_mask <- ...  # reconstruct the original bounding-box mask for comparison
# compare_rural_definitions(lst_grid, old_rural_mask, new_rural_mask)
#
# knime.out <- data.frame(status = "rural_mask_redefined")

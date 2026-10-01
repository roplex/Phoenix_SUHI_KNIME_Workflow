# --- KNIME: Multi-Raster Zonal Statistics ---
# Chapter 7 (SUHI) -- corrected version
import os
import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import rasterio
import knime.scripting.io as knio
import warnings

# Silence harmless pandas/geopandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# --- Read KNIME input table ---
input_df = knio.input_tables[0].to_pandas()

# Extract raster paths from input table
lst_path = input_df["lst_path"].iloc[0]
ndvi_path = input_df["ndvi_path"].iloc[0]
suhi_path = input_df["suhi_path"].iloc[0]

# --- Path to vector file (villages) ---
# CORRECTED from .../PhoenixData/Boundary/Villages.geojson to
# .../PhoenixData2/Boundary/Villages.geojson -- every verified run in this
# pipeline (Stratified NDVI-LST Regression, the village-confounding
# diagnostic) has used PhoenixData2; the original PhoenixData (no "2") path
# is stale from before the data reorganization.
villages_path = "/Users/roplex/Desktop/EO_Harmonization/PhoenixData2/Boundary/Villages.geojson"

# --- File-existence checks -- none existed in the original script ---
for label, p in [("lst_path", lst_path), ("ndvi_path", ndvi_path),
                  ("suhi_path", suhi_path), ("villages_path", villages_path)]:
    if not os.path.exists(p):
        raise FileNotFoundError(f"{label} does not exist: {p}")

# --- Load vector data ---
villages = gpd.read_file(villages_path)
if "NAME" not in villages.columns:
    raise ValueError(f"Expected a 'NAME' field in villages_path but found: {list(villages.columns)}")

n_villages_in = len(villages)
print(f"Loaded {n_villages_in} village polygons from {villages_path}")

# --- Known per-village pixel counts from the Stratified NDVI-LST Regression
# node's verified TsHARP run (n_points there = pixels where LST, NDVI, AND
# SUHI are ALL simultaneously valid). LST_count/SUHI_count below should be
# >= these values (a single-raster count is a superset of the three-way
# intersection); NDVI_count may be noticeably higher for every village,
# since raw NDVI has more valid pixels citywide than LST/SUHI (see the
# village-confounding diagnostic: 31,255 raw NDVI vs ~19,194 LST/SUHI). Use
# this as a sanity floor, not an exact-match requirement. ---
known_stratified_n_points = {
    "Ahwatukee Foothills": 1582, "Alhambra": 4, "Camelback East": 338,
    "Central City": 89, "Deer Valley": 1835, "Desert View": 3809,
    "Encanto": 0, "Estrella": 2049, "Laveen": 1472, "Maryvale": 1009,
    "North Gateway": 2570, "North Mountain": 490, "Paradise Valley": 481,
    "Rio Vista": 1619, "South Mountain": 1417,
}

# --- Helper to compute and return zonal stats DataFrame ---
def compute_zonal_stats(raster_path, prefix):
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        raster_nodata = src.nodata
    print(f"[{prefix}] {os.path.basename(raster_path)} -- nodata={raster_nodata}, crs={raster_crs}")
    if raster_nodata is None:
        print(f"  WARNING: {prefix} raster has no NoData value set in its header. "
              f"zonal_stats(nodata=None) will rely on this being genuinely absent -- "
              f"verify this isn't silently including invalid pixels.")

    if villages.crs != raster_crs:
        v_proj = villages.to_crs(raster_crs)
    else:
        v_proj = villages.copy()

    stats = zonal_stats(
        vectors=v_proj,
        raster=raster_path,
        stats=["count", "mean", "median", "std", "min", "max"],
        nodata=None,
        geojson_out=False
    )
    df = pd.DataFrame(stats)
    if len(df) != n_villages_in:
        raise RuntimeError(f"{prefix}: zonal_stats returned {len(df)} rows, expected {n_villages_in}")
    df = df.add_prefix(f"{prefix}_")
    return df

# --- Compute stats for each raster ---
stats_lst  = compute_zonal_stats(lst_path,  "LST")
stats_ndvi = compute_zonal_stats(ndvi_path, "NDVI")
stats_suhi = compute_zonal_stats(suhi_path, "SUHI")

# --- Combine with village attributes ---
villages_df = villages.reset_index(drop=True)
combined_df = pd.concat([villages_df, stats_lst, stats_ndvi, stats_suhi], axis=1)

# --- Cross-check LST_count/SUHI_count against the known Stratified n_points ---
print("\nCross-check against Stratified Regression node's known n_points:")
for _, row in combined_df.iterrows():
    name = row["NAME"]
    expected = known_stratified_n_points.get(name)
    if expected is None:
        print(f"  {name}: not in reference table -- skipping cross-check")
        continue
    lst_count = row["LST_count"]
    suhi_count = row["SUHI_count"]
    flag = ""
    if pd.isna(lst_count) or lst_count < expected or pd.isna(suhi_count) or suhi_count < expected:
        flag = "  <-- LOWER than expected floor, worth investigating"
    print(f"  {name}: LST_count={lst_count}, SUHI_count={suhi_count}, "
          f"expected >= {expected}{flag}")

# --- Flag villages with any missing (all-NaN) stats, e.g. Encanto (0 valid
# pixels, a known, pre-existing limitation reproduced identically in the
# original bilinear analysis -- not introduced by any of this pipeline's
# corrections) ---
missing = combined_df[combined_df[["LST_mean", "NDVI_mean", "SUHI_mean"]].isna().any(axis=1)]
if len(missing) > 0:
    print(f"\n{len(missing)} village(s) with at least one missing raster mean:")
    print(missing[["NAME", "LST_mean", "NDVI_mean", "SUHI_mean"]])
else:
    print("\nNo villages have missing raster means.")

# --- Save outputs ---
# NOTE: filenames don't currently disambiguate which LST resampling method
# (bilinear vs TsHARP) fed this run -- a fixed-filename mixup between the
# two already caused one false alarm earlier in this pipeline (Global/
# Stratified Regression nodes). Worth appending a method tag here too, e.g.
# "villages_zonal_stats_py_tsharp.csv", once you've confirmed which run
# this is from.
out_dir = os.path.dirname(lst_path)
csv_out = os.path.join(out_dir, "villages_zonal_stats_py.csv")
geojson_out = os.path.join(out_dir, "villages_zonal_stats_py.geojson")

combined_df.to_csv(csv_out, index=False)
combined_df.to_file(geojson_out, driver="GeoJSON")

# --- Return outputs to KNIME ---
knio.output_tables[0] = knio.Table.from_pandas(combined_df)
knio.flow_variables["csv_path"] = csv_out
knio.flow_variables["geojson_path"] = geojson_out
knio.flow_variables["n_villages"] = len(combined_df)

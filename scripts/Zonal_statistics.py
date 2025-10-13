# --- KNIME: Multi-Raster Zonal Statistics ---
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

# Path to vector file (villages)
villages_path = "/Users/roplex/Desktop/EO_Harmonization/PhoenixData/Boundary/Villages.geojson"

# --- Load vector data ---
villages = gpd.read_file(villages_path)

# --- Helper to compute and return zonal stats DataFrame ---
def compute_zonal_stats(raster_path, prefix):
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
    if villages.crs != raster_crs:
        v_proj = villages.to_crs(raster_crs)
    else:
        v_proj = villages.copy()
    
    stats = zonal_stats(
        vectors=v_proj,
        raster=raster_path,
        stats=["mean", "median", "std", "min", "max"],
        nodata=None,
        geojson_out=False
    )
    df = pd.DataFrame(stats)
    df = df.add_prefix(f"{prefix}_")
    return df

# --- Compute stats for each raster ---
stats_lst  = compute_zonal_stats(lst_path,  "LST")
stats_ndvi = compute_zonal_stats(ndvi_path, "NDVI")
stats_suhi = compute_zonal_stats(suhi_path, "SUHI")

# --- Combine with village attributes ---
villages_df = villages.reset_index(drop=True)
combined_df = pd.concat([villages_df, stats_lst, stats_ndvi, stats_suhi], axis=1)

# --- Save outputs ---
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

import cdsapi
import os

client = cdsapi.Client()

# Absolute path to your Data folder
data_dir = os.path.expanduser("~/Downloads/Stats Modeling project/Data")
target   = os.path.join(data_dir, "era5_monthly.nc")

client.retrieve(
    "reanalysis-era5-single-levels-monthly-means",
    {
        "product_type": "monthly_averaged_reanalysis",
        "variable": ["2m_temperature", "10m_wind_speed"],
        "year": [str(y) for y in range(2015, 2021)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": "00:00",
        "area": [80, 180, 40, -180],   # N, W, S, E  — covers both sides of the dateline
        "format": "netcdf"
    },
    target
)

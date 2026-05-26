# ============================================================
# ERA5 data extraction for the Pacific Arctic hex grid
# ============================================================
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(ggplot2)

data_path <- "Downloads/Stats Modeling project/Data"

# ---- Load data (once) ----
era5     <- rast(file.path(data_path, "era5_monthly.nc"))
hex_grid <- readRDS(file.path(data_path, "hex_grid_template.rds"))

# Separate temperature and wind layers
temp_layers <- era5[[1:72]]
wind_layers <- era5[[73:144]]

# ------------------------------------------------------------
# Chunk 1 – Single‑month extraction (January 2015) for quick check
# ------------------------------------------------------------
# Reproject the FIRST temperature layer (January 2015)
temp_jan <- project(temp_layers[[1]], "EPSG:3338", method = "bilinear")

# Extract area‑weighted mean per hex for this single month
vals <- exact_extract(temp_jan, hex_grid, fun = "mean", progress = FALSE)

jan2015_temp <- data.frame(
  hexID    = hex_grid$hexID,
  air_temp = vals - 273.15   # Kelvin → Celsius
)

# Quick checks
cat("Number of hex cells:", nrow(jan2015_temp), "\n")
cat("Number of NA values:", sum(is.na(jan2015_temp$air_temp)), "\n")
cat("Summary of air temperature (Celsius):\n")
summary(jan2015_temp$air_temp)

# Map the result
hex_map <- left_join(hex_grid, jan2015_temp, by = "hexID")
ggplot(hex_map) +
  geom_sf(aes(fill = air_temp), colour = NA, size = 0.05) +
  scale_fill_viridis_c(option = "plasma", name = "air temp (°C)") +
  coord_sf(crs = st_crs(hex_map)) +
  theme_void() +
  ggtitle("ERA5 2m air temperature – January 2015 (hex average)")

# ------------------------------------------------------------
# Chunk 2 – Extract temperature & wind for ALL 72 months
# ------------------------------------------------------------
# Use the already loaded era5 and hex_grid (no reload needed)

temp_3338 <- project(temp_layers, "EPSG:3338", method = "bilinear")
wind_3338 <- project(wind_layers, "EPSG:3338", method = "bilinear")

# Correct date extraction from NetCDF time units
secs <- as.numeric(depth(temp_layers))                 # seconds since 1970-01-01
month_dates <- as.Date(as.POSIXct(secs, origin = "1970-01-01", tz = "UTC"))

temp_df <- NULL
wind_df <- NULL

for (m in seq_along(month_dates)) {
  val_temp <- exact_extract(temp_3338[[m]], hex_grid, fun = "mean", progress = FALSE)
  val_wind <- exact_extract(wind_3338[[m]], hex_grid, fun = "mean", progress = FALSE)
  
  temp_df <- bind_rows(temp_df, data.frame(
    hexID    = hex_grid$hexID,
    date     = month_dates[m],
    air_temp = val_temp - 273.15
  ))
  wind_df <- bind_rows(wind_df, data.frame(
    hexID      = hex_grid$hexID,
    date       = month_dates[m],
    wind_speed = val_wind
  ))
  cat("Month", m, "of", length(month_dates), "done\n")
}

# Save the tables
saveRDS(temp_df, file.path(data_path, "air_temp.rds"))
saveRDS(wind_df, file.path(data_path, "wind_speed.rds"))
cat("Saved air_temp.rds and wind_speed.rds with correct dates\n")

# Quick check on dates (no need to re‑read; temp_df is already here)
cat("First 6 dates:\n")
print(head(temp_df$date))
cat("Unique dates:", length(unique(temp_df$date)), "\n")

# ------------------------------------------------------------
# Faceted temperature maps – winter vs summer for 3 years
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(sf)

# Load the saved files (if this script is run in a fresh session)
temp_df <- readRDS(file.path(data_path, "air_temp.rds"))
hex_grid <- readRDS(file.path(data_path, "hex_grid_template.rds"))

selected <- expand.grid(
  year   = c(2015, 2018, 2020),
  month  = c("01", "07")
) %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

months_plot <- temp_df %>%
  filter(date %in% selected$date) %>%
  left_join(hex_grid, by = "hexID") %>%
  st_as_sf() %>%
  left_join(selected, by = "date") %>%
  mutate(season = ifelse(month == "01", "Winter (Jan)", "Summer (Jul)"))

ggplot(months_plot) +
  geom_sf(aes(fill = air_temp), colour = NA, size = 0.05) +
  scale_fill_viridis_c(option = "plasma", name = "air temp (°C)") +
  facet_grid(rows = vars(season), cols = vars(year)) +
  coord_sf(crs = st_crs(months_plot)) +
  theme_void() +
  labs(title = "ERA5 2m air temperature – Winter vs Summer")

# ------------------------------------------------------------
# Faceted WIND speed maps – winter vs summer for 3 years
# ------------------------------------------------------------
wind_df <- readRDS(file.path(data_path, "wind_speed.rds"))
hex_grid <- readRDS(file.path(data_path, "hex_grid_template.rds"))

selected <- expand.grid(
  year   = c(2015, 2018, 2020),
  month  = c("01", "07")
) %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

months_plot <- wind_df %>%
  filter(date %in% selected$date) %>%
  left_join(hex_grid, by = "hexID") %>%
  st_as_sf() %>%
  left_join(selected, by = "date") %>%
  mutate(season = ifelse(month == "01", "Winter (Jan)", "Summer (Jul)"))

ggplot(months_plot) +
  geom_sf(aes(fill = wind_speed), colour = NA, size = 0.05) +
  scale_fill_viridis_c(option = "plasma", name = "wind speed (m/s)") +
  facet_grid(rows = vars(season), cols = vars(year)) +
  coord_sf(crs = st_crs(months_plot)) +
  theme_void() +
  labs(title = "ERA5 10m wind speed – Winter vs Summer")

# ------------------------------------------------------------
# Merge Datasets
# ------------------------------------------------------------
ship_ice <- readRDS(file.path(data_path, "merged_shipping_ice.rds"))
temp     <- readRDS(file.path(data_path, "air_temp.rds"))
wind     <- readRDS(file.path(data_path, "wind_speed.rds"))

# Join temperature and wind
df <- ship_ice %>%
  left_join(temp, by = c("hexID", "date")) %>%
  left_join(wind, by = c("hexID", "date"))

# Save the complete analytic dataset (dates are already correct)
saveRDS(df, file.path(data_path, "analytic_dataset.rds"))

cat("Merged dimensions:", nrow(df), "rows,", ncol(df), "columns\n")
cat("Date range:", min(df$date, na.rm = TRUE), "to", max(df$date, na.rm = TRUE), "\n")

# Summary of the four key variables
df %>% select(nOprD_A, ice_conc, air_temp, wind_speed) %>% summary()
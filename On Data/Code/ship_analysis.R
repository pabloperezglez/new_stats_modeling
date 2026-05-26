library(sf)
library(dplyr)
library(purrr)
library(lubridate)
library(tidyr)
library(ggplot2)


# Inspect One Month 

# Path to the dataset (adjust if needed)
base_path <- "Downloads/Stats_Modeling_Data/North Pacific and Arctic Marine Vessel Traffic Dataset"


# There is an unzipped example folder "SpeedHex_2015_01"
one_shp <- file.path(base_path, "SpeedHex_2015_01", "SpeedHex_2015-01.shp")
gdf_jan15 <- st_read(one_shp, quiet = TRUE)

cat("Dimensions:", nrow(gdf_jan15), "rows,", ncol(gdf_jan15), "columns\n")
cat("CRS:", st_crs(gdf_jan15)$input, "\n\n")
print(head(gdf_jan15))


# Load and Stack all 72 months

zip_files <- list.files(base_path, pattern = "SpeedHex_.*\\.zip$", full.names = TRUE)
cat("Found", length(zip_files), "ZIP files.\n")

read_speedhex <- function(zip_path) {
  # Unzip to a temp directory
  tmp <- tempfile()
  dir.create(tmp)
  unzip(zip_path, exdir = tmp)
  shp <- list.files(tmp, pattern = "\\.shp$", full.names = TRUE)[1]
  gdf <- st_read(shp, quiet = TRUE)
  unlink(tmp, recursive = TRUE)
  
  gdf <- gdf %>%
    mutate(
      date  = ymd(paste(year, month, "01", sep = "-")),
      year  = as.integer(year),
      month = as.integer(month)
    )
  return(gdf)
}

# stacking (may take a minute)
options(readr.show_col_types = FALSE) # suppress readr messages
shipping_sf <- map_dfr(zip_files, read_speedhex)

cat("Stacked dataset:", nrow(shipping_sf), "rows,", ncol(shipping_sf), "columns\n")
cat("Time span:", format(min(shipping_sf$date)), "to", format(max(shipping_sf$date)), "\n")



# Recode Missing Values Count variables: NA → 0 (no vessel). Speed variables (SOGs) remain NA

count_vars <- grep("^(nMMSI|nOprD)", names(shipping_sf), value = TRUE)
shipping_sf <- shipping_sf %>%
  mutate(across(all_of(count_vars), ~ replace_na(.x, 0)))

# quick check: number of zeros in nOprD_A
cat("Zeros in nOprD_A:", sum(shipping_sf$nOprD_A == 0),
    "/", nrow(shipping_sf), "\n")
cat("Proportion zeros:", round(mean(shipping_sf$nOprD_A == 0), 3), "\n")




# Basic Temporal view-total month traffic

monthly_traffic <- shipping_sf %>%
  st_drop_geometry() %>%
  group_by(date) %>%
  summarise(total_traffic = sum(nOprD_A, na.rm = TRUE), .groups = "drop")

ggplot(monthly_traffic, aes(x = date, y = total_traffic)) +
  geom_line() + geom_point() +
  labs(title = "Monthly total shipping intensity (nOprD_A)",
       x = "Date", y = "Sum of nOprD_A") +
  theme_minimal()



# MAPS

# select four months
months_to_plot <- c("2018-01-01", "2018-07-01", "2020-01-01", "2020-07-01")
map_data <- shipping_sf %>%
  filter(as.character(date) %in% months_to_plot)

# use a log scale for better visibility
map_data <- map_data %>%
  mutate(traffic_cat = cut(nOprD_A,
                           breaks = c(0, 1, 10, 50, 100, 500, Inf),
                           labels = c("1","2-10","11-50","51-100","101-500",">500"),
                           include.lowest = TRUE))

ggplot(map_data) +
  geom_sf(aes(fill = traffic_cat), color = NA, size = 0.05) +
  facet_wrap(~ date, ncol = 2) +
  scale_fill_viridis_d(option = "plasma", na.value = "grey80",
                       name = "nOprD_A") +
  theme_void() +
  labs(title = "Shipping intensity (nOprD_A) across selected months")
       

# For better visualization
ggplot(map_data) +
  geom_sf(aes(fill = log1p(nOprD_A)), color = NA, size = 0.05) +
  facet_wrap(~ date, ncol = 2) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey80",
                       name = "log(1 + nOprD_A)") +
  theme_void() +
  labs(title = "Shipping intensity (log scale)") 




# How many areas recieve 0 traffic by month/year?


zero_prop <- shipping_sf %>%
  st_drop_geometry() %>%
  group_by(year, month) %>%
  summarise(zero_frac = mean(nOprD_A == 0), .groups = "drop")

ggplot(zero_prop, aes(x = month, y = zero_frac, color = factor(year))) +
  geom_line() +
  labs(y = "Proportion of cells with zero traffic",
       title = "Zero inflation by month and year") +
  theme_minimal()



## 8. Save the cleaned dataset (two versions) into the Data folder

data_path <- "Downloads/Stats Modeling project/Data"

# 8a) Full spatial stack
saveRDS(shipping_sf,                                    file = file.path(data_path, "shipping_sf_full.rds"))

# 8b) Attribute-only table
saveRDS(sf::st_drop_geometry(shipping_sf),             file = file.path(data_path, "cleaned_shipping.rds"))

# 8c) Hex grid template (6 553 polygons)
saveRDS(dplyr::filter(shipping_sf, date == lubridate::ymd("2015-01-01")) %>% select(hexID, geometry),
        file = file.path(data_path, "hex_grid_template.rds"))

cat("Shipping files written to", data_path, "\n")
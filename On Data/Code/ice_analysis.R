# ============================================================
# ice_analysis.R  –  Cleaned, Working Version
# ============================================================

library(icecdr)
library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(dplyr)
library(exactextractr)

# ---- Where your project files live ----
data_path <- "Downloads/Stats Modeling project/Data"

# ---- 1. Download monthly Arctic ice concentration ----
cat("===== Downloading monthly Arctic ice data =====\n")
ice_file <- cdr_arctic_monthly(
  date_range = c("2015-01-01", "2020-12-01"),
  variables  = "aice",
  version    = 5,
  dir        = "ice_data"
)
cat("File downloaded to:", ice_file, "\n")

# ---- 2. Load raster stack & assign CRS ----
ice_stack <- rast(ice_file)
crs(ice_stack) <- "EPSG:3413"      # NSIDC North Polar Stereographic

cat("\n===== Ice Stack Summary =====\n")
print(ice_stack)
cat("Number of layers (months):", nlyr(ice_stack), "\n")
cat("Resolution (x, y):", res(ice_stack), "\n")
cat("Extent:\n")
print(ext(ice_stack))

# ---- 3. Inspect January 2015 (layer 1) with proper min/max ----
jan2015 <- ice_stack[[1]]
names(jan2015) <- "2015-01"

mm <- global(jan2015, fun = range, na.rm = TRUE)
cat("Min ice conc (%):", mm[1,1], "  Max ice conc (%):", mm[1,2], "\n")
cat("Number of NA cells:", sum(is.na(values(jan2015))), "\n")

# ---- 4. Six‑panel map (Jan/Jul 2015, 2018, 2020) ----
idx <- c(1, 7, 37, 43, 61, 67)
months_to_plot <- c("2015-01", "2015-07", "2018-01", "2018-07",
                    "2020-01", "2020-07")
plot(ice_stack[[idx]], col = rev(viridis(100)), main = months_to_plot)

# ---- 5. Hex grid overlay on January 2015 ----
hex_grid <- readRDS(file.path(data_path, "hex_grid_template.rds"))
hex_ice  <- st_transform(hex_grid, crs = "EPSG:3413")

jan_df <- as.data.frame(jan2015, xy = TRUE)
colnames(jan_df)[3] <- "ice_conc"

ggplot() +
  geom_raster(data = jan_df, aes(x = x, y = y, fill = ice_conc)) +
  scale_fill_viridis_c(name = "ice conc (%)", na.value = "transparent") +
  geom_sf(data = hex_ice, fill = NA, color = "white", linewidth = 0.1) +
  coord_sf(crs = "EPSG:3413") +
  ggtitle("Hex grid over January 2015 ice concentration") +
  theme_minimal()

# ---- 6. Centroid extraction check ----
hex_centroids <- st_centroid(hex_ice)
ice_at_centroids <- extract(jan2015, hex_centroids, fun = mean, na.rm = TRUE)
hex_ice$ice_jan2015 <- ice_at_centroids[,2]
summary(hex_ice$ice_jan2015)

# ============================================================
# FULL EXTRACTION – Ice concentration for every hex × month
# ============================================================

# Reproject ice stack to match hex grid CRS (EPSG:3338)
ice_3338 <- project(ice_stack, "EPSG:3338", method = "bilinear")

n_months <- nlyr(ice_3338)
ice_df   <- NULL

cat("Extracting ice for", n_months, "months ...\n")
for (m in 1:n_months) {
  vals <- exact_extract(ice_3338[[m]], hex_grid, fun = "mean", progress = FALSE)
  month_date <- time(ice_3338)[m]
  df_month <- data.frame(
    hexID    = hex_grid$hexID,
    date     = month_date,
    ice_conc = vals
  )
  ice_df <- bind_rows(ice_df, df_month)
  cat("Month", m, "done\n")
}

# Save the ice‑concentration table
saveRDS(ice_df, file = file.path(data_path, "ice_concentration.rds"))
cat("Ice concentration table (", nrow(ice_df), " rows) saved\n")

summary(ice_df$ice_conc)

# ============================================================
# VISUALISATIONS – Hex‑level ice maps for selected months
# ============================================================

hex_sf <- hex_grid

map_months <- c("2015-01", "2015-07", "2018-01", "2018-07", "2020-01", "2020-07")
map_dates  <- as.Date(paste0(map_months, "-01"))

# Individual hex maps
for (i in seq_along(map_months)) {
  month_str <- map_months[i]
  d <- as.Date(paste0(month_str, "-01"))
  month_ice <- ice_df %>% filter(date == d)
  hex_map <- left_join(hex_sf, month_ice, by = "hexID")
  
  p <- ggplot(hex_map) +
    geom_sf(aes(fill = ice_conc), color = NA, size = 0.05) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey80",
                         name = "ice conc (%)", limits = c(0, 1)) +
    coord_sf(crs = st_crs(hex_map)) +
    ggtitle(paste("Ice concentration on hex grid –", month_str)) +
    theme_void()
  
  file_name <- paste0("ice_hex_", gsub("-", "_", month_str), ".png")
  ggsave(file.path(data_path, file_name), plot = p,
         width = 8, height = 8, dpi = 150)
  cat("Saved map:", file_name, "\n")
}

# Faceted map
hex_maps <- lapply(map_dates, function(d) {
  month_ice <- ice_df %>% filter(date == d)
  left_join(hex_sf, month_ice, by = "hexID") %>%
    mutate(month = as.character(d))
}) %>% bind_rows()

p_facet <- ggplot(hex_maps) +
  geom_sf(aes(fill = ice_conc), color = NA, size = 0.05) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey80",
                       name = "ice conc (%)", limits = c(0, 1)) +
  facet_wrap(~ month, ncol = 3) +
  theme_void() +
  labs(title = "Hex‑level ice concentration (selected months)")

p_facet

ggsave(file.path(data_path, "ice_hex_faceted.png"), plot = p_facet,
       width = 14, height = 10, dpi = 150)
cat("Faceted map saved.\n")


# ---------------------------------------------------------------
# Zoomed hex map – April 2019 (2× zoom, shifted top‑right)
# ---------------------------------------------------------------

# ---- 1. Load data ----
hex_grid <- readRDS(file.path(data_path, "hex_grid_template.rds"))
ice_df   <- readRDS(file.path(data_path, "ice_concentration.rds"))

# ---- 2. April 2019 ----
april2019 <- ice_df %>% filter(date == as.Date("2019-04-01"))
hex_april <- left_join(hex_grid, april2019, by = "hexID")

# ---- 3. Define zoom window (2× magnification, top‑right shift) ----
# Original extent (from shipping shapefile inspection):
# x: -2565510 ... 559770   -> width ≈ 3 115 280 m
# y:   213432 ... 2738670   -> height ≈ 2 525 238 m
# Centroid: x = −1 002 870, y = 1 476 051

zoom_factor <- 2
shift_x     <- 0.15 * (559770 - (-2565510))   # shift to the right
shift_y     <- 0.15 * (2738670 - 213432)       # shift toward the top

new_center_x <- -1002870 + shift_x
new_center_y <-  1476051 + shift_y

half_w <- (559770 - (-2565510)) / zoom_factor / 2
half_h <- (2738670 - 213432) / zoom_factor / 2

xlim <- c(new_center_x - half_w, new_center_x + half_w)
ylim <- c(new_center_y - half_h, new_center_y + half_h)

# ---- 4. Plot ----
ggplot(hex_april) +
  geom_sf(aes(fill = ice_conc), colour = NA, size = 0.1) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey80",
                       name   = "ice concentration\n(April 2019)",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent) +
  coord_sf(crs   = st_crs(hex_april),
           xlim  = xlim, ylim = ylim,
           expand = FALSE) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Sea‑ice concentration – April 2019 (zoomed)")


# ---------------------------------------------------------------
# Merge ice with shipping, then save the combined dataset
# ---------------------------------------------------------------

data_path <- "Downloads/Stats Modeling project/Data"

ice  <- readRDS(file.path(data_path, "ice_concentration.rds"))
ship <- readRDS(file.path(data_path, "cleaned_shipping.rds"))

df <- left_join(ship, ice, by = c("hexID", "date"))

# Quick check
cat("Merged dimensions:", nrow(df), "rows,", ncol(df), "columns\n")
head(df)

# Save the merged data frame
saveRDS(df, file = file.path(data_path, "merged_shipping_ice.rds"))
cat("Merged shipping + ice saved as 'merged_shipping_ice.rds'\n")



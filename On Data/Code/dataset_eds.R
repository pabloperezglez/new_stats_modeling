library(dplyr)

df <- readRDS("Downloads/Stats Modeling project/Data/analytic_dataset.rds")

glimpse(df)
# Check that these exist (they should based on your three pipelines):
#   hexID, date, year, month
#   nOprD_A, nOprD_C, nOprD_T, nOprD_O, nOprD_F, nOprD_L
#   nMMSI_A, nMMSI_C, ... (optional for later)
#   SOG_A, SOG_C, ... (optional for later)
#   ice_conc, air_temp, wind_speed
#   geometry (sf polygon, may or may not be loaded)

# Quick summary of missingness
df %>% summarise(across(everything(), ~ mean(is.na(.x)))) %>% glimpse()


# convert NA in the count columns to 0
count_cols <- grep("n(MMSI|OprD)_", names(df), value = TRUE)
df <- df %>% mutate(across(all_of(count_cols), ~ ifelse(is.na(.x), 0, .x)))

library(ggplot2)
library(mgcv)

# a) Raw histogram
ggplot(df, aes(nOprD_A)) +
  geom_histogram(binwidth = 1, fill = "steelblue") +
  scale_y_log10() +
  labs(title = "Vessel‑days per hex‑month (all vessels)", y = "Count (log10)")

# b) Same, but only when traffic is > 0
df_traffic <- df %>% filter(nOprD_A > 0)
ggplot(df_traffic, aes(nOprD_A)) +
  geom_histogram(binwidth = 2, fill = "steelblue") +
  scale_x_log10() +
  labs(title = "Vessel‑days when traffic present", x = "nOprD_A (log scale)")

# some scatter plot

ggplot(df, aes(ice_conc, nOprD_A)) +
  geom_point(alpha = 0.03, size = 0.1) +          # many points → tiny alpha
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 15),
              method.args = list(family = tw(link = "log")),
              colour = "firebrick", fill = "firebrick") +
  labs(title = "Traffic intensity vs. sea‑ice concentration (all months)")

# the same but only when ice is present:
df_ice <- df %>% filter(ice_conc >= 0.15)

ggplot(df_ice, aes(ice_conc, nOprD_A)) +
  geom_point(alpha = 0.04, size = 0.1) +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, bs = "ps", k = 10),
              method.args = list(family = tw(link = "log")),
              colour = "darkorange", fill = "darkorange") +
  coord_cartesian(ylim = c(0, 50)) +   # clip, but the smooth still uses all data
  labs(title = "Traffic intensity vs. sea‑ice (ice ≥ 15 %)")

# Monthly Totals for all Vessels

monthly <- df %>%
  group_by(year, month) %>%
  summarise(total_traffic = sum(nOprD_A), .groups = "drop")

ggplot(monthly, aes(as.Date(paste(year, month, "01", sep = "-")), total_traffic)) +
  geom_line(colour = "steelblue") +
  labs(x = "Date", y = "Total vessel‑days per month")

# Seasonal cycle by vessel type

library(tidyr)

monthly_type <- df %>%
  group_by(year, month) %>%
  summarise(
    Cargo = sum(nOprD_C),
    Tanker = sum(nOprD_T),
    Other = sum(nOprD_O),
    Fishing = sum(nOprD_F),
    .groups = "drop"
  ) %>%
  pivot_longer(Cargo:Fishing, names_to = "vessel_type", values_to = "total_traffic")

ggplot(monthly_type, aes(as.Date(paste(year, month, "01", sep = "-")), 
                         total_traffic, colour = vessel_type)) +
  geom_line() +
  facet_wrap(~ vessel_type, scales = "free_y") +
  labs(x = "Date", y = "Monthly vessel‑days")


# Vessel Type composition: 

type_totals <- df %>%
  summarise(
    Cargo   = sum(nOprD_C),
    Tanker  = sum(nOprD_T),
    Other   = sum(nOprD_O),
    Fishing = sum(nOprD_F)
  ) %>%
  pivot_longer(everything(), names_to = "type", values_to = "total")

# Bar chart instead of pie (easier to read)
ggplot(type_totals, aes(type, total / 1e6)) +
  geom_col(fill = "steelblue") +
  labs(y = "Total vessel‑days (millions)") +
  theme_minimal()

# Pie Chart
ggplot(type_totals, aes("", total, fill = type)) +
  geom_col() +
  coord_polar("y") +
  labs(fill = "Vessel type")


# joint distribution of ice and temperature:
ggplot(df, aes(air_temp, ice_conc)) +
  geom_bin2d(bins = 70) +
  scale_fill_viridis_c(option = "magma", trans = "log10") +
  labs(title = "Ice concentration vs. 2m air temperature")


# Wind speed and ice concentration

ggplot(df, aes(ice_conc, wind_speed)) +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c(option = "magma", trans = "log10") +
  labs(title = "Wind speed vs. sea‑ice concentration")

# Mean traffic over all months

library(sf)

# Load hex grid (or use the geometry stored in df if it's still an sf object)
hex_grid <- readRDS("Downloads/Stats Modeling project/Data/hex_grid_template.rds")  # EPSG:3338

# Compute mean traffic per hex
hex_means <- df %>%
  group_by(hexID) %>%
  summarise(mean_traffic = mean(nOprD_A), .groups = "drop") %>%
  left_join(hex_grid, by = "hexID") %>%
  st_as_sf()

ggplot(hex_means) +
  geom_sf(aes(fill = mean_traffic), colour = NA) +
  scale_fill_viridis_c(option = "plasma", trans = "log1p") +
  labs(title = "Mean monthly vessel‑days per hex (2015–2020)")


# Low‑effort linear/GLM view (just for intuition)

# Quick and dirty: aggregate by ice concentration bins
df %>%
  mutate(ice_bin = cut(ice_conc, breaks = seq(0, 1, by = 0.1))) %>%
  group_by(ice_bin) %>%
  summarise(mean_ndays = mean(nOprD_A), n = n()) %>%
  ggplot(aes(ice_bin, mean_ndays)) +
  geom_point() +
  geom_smooth(aes(as.numeric(ice_bin)), method = "lm", se = FALSE) +
  labs(x = "Ice concentration bin", y = "Average vessel‑days")
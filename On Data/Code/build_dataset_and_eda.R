# =============================================================================
#  build_dataset_and_eda.R
#  ---------------------------------------------------------------------------
#  SINGLE, self-contained data-construction + EDA script for the Pacific Arctic
#  shipping study (2015-2020). Consolidates the previously scattered scripts
#  (ship_analysis.R, ice_analysis.R, air_analysis.R, download_era5.py,
#  dataset_eds.R) into ONE coherent, narrated pipeline.
#
#  Deliverable scope: DATA retrieval + transformation + EDA only.
#  NO predictive / regression / GAM analysis is fit here. The ONLY regression
#  that appears (Section 9) is a univariate smoother used purely to CONSTRUCT
#  an orthogonalized-temperature covariate; it is a data transformation, not a
#  model of the response.
#
#  Output: a tidy analytic table (one row per hexID x month) plus a full set of
#  saved EDA figures and verbose, copy-pasteable console output (ledger, NA
#  tallies, port table, row-loss localization, correlations, orthogonalization
#  diagnostics) intended to be pasted into the LaTeX document
#  "Data Handling for Melting Ports".
#
#  Author: <PLACEHOLDER>            Date: Sys.Date() at run time
# =============================================================================

# =============================================================================
## 0. CONFIG  — every path / CRS / date-range / bbox lives here, nothing mid-script
# =============================================================================

set.seed(2025)  # determinism (only the orthogonalization spline is "fitty")

# --- self-locate the project root --------------------------------------------
# This script lives in <project>/On Data/Code/. Data live in ../Data and the
# outputs go to ../Data/output. Paths are resolved relative to this file so the
# project is portable (no hard-coded ~/Downloads).
.this_file <- local({
  for (i in seq_len(sys.nframe())) {
    fr <- sys.frame(i); if (!is.null(fr$ofile)) return(normalizePath(fr$ofile))
  }
  a <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", a, value = TRUE)
  if (length(fa)) return(normalizePath(sub("^--file=", "", fa[1])))
  NA_character_
})
CODE_DIR     <- if (!is.na(.this_file)) dirname(.this_file) else getwd()
PROJECT_ROOT <- normalizePath(file.path(CODE_DIR, ".."))   # the "On Data" folder

CONFIG <- list(
  # --- directories (resolved relative to the project root) -----------------
  data_dir = file.path(PROJECT_ROOT, "Data"),                # tables / csv / rds / nc
  out_dir  = file.path(PROJECT_ROOT, "Data", "output"),      # logs + audit csv
  fig_dir  = file.path(PROJECT_ROOT, "Data", "output", "figures"),

  # --- study definition ----------------------------------------------------
  date_start = as.Date("2015-01-01"),
  date_end   = as.Date("2020-12-01"),
  n_months_expected = 72L,
  n_hex_expected    = 6553L,

  # --- coordinate reference systems ---------------------------------------
  crs_ship = 3338,   # NAD83 / Alaska Albers (equal-area, metres) — analysis CRS
  crs_ice  = 3413,   # NSIDC North Polar Stereographic (native ice CRS)
  crs_geo  = 4326,   # WGS84 lon/lat (for geodesic port distances)

  # --- ERA5 download bounding box (N, W, S, E) used by the cdsapi/ecmwfr call
  era5_bbox = c(N = 80, W = 180, S = 40, E = -180),

  # --- ice product pin (resolve CDR-v5 vs Sea-Ice-Index naming) ------------
  ice_product = "NOAA/NSIDC Climate Data Record (CDR) of passive-microwave sea-ice concentration, Version 5 (variable 'aice'); monthly; 25 km polar-stereographic. Retrieved via the icecdr R package (cdr_arctic_monthly).",

  # --- structural-zero / derived-flag thresholds ---------------------------
  ice_free_threshold = 0.15,         # conventional 'open water' cut
  polar_code_date    = as.Date("2017-01-01")  # IMO Polar Code in force (used only as a pre/post indicator)
)

dir.create(CONFIG$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$fig_dir, recursive = TRUE, showWarnings = FALSE)

# convenience path builders
dpath <- function(...) file.path(CONFIG$data_dir, ...)

# fpath(): figures are written with descriptive, human-readable names. The body
# of the script refers to figures by short codenames (fig01..fig16); this lookup
# translates each codename to its descriptive filename so the saved PNGs are
# self-explanatory.
.FIG_RENAME <- c(
  "fig01_zero_prop_by_ice_bin.png"          = "01_zero-proportion-of-traffic_by_ice-concentration-bin.png",
  "fig02_monthly_traffic_timeseries.png"    = "02_monthly-total-vessel-days_timeseries.png",
  "fig03_seasonal_cycle_by_vessel_type.png" = "03_seasonal-cycle_by-vessel-type.png",
  "fig04_vessel_type_composition.png"       = "04_vessel-type-composition_donut.png",
  "fig05_ice_hex_jan2015.png"               = "05_sea-ice_hex-aggregated_jan2015.png",
  "fig06_faceted_ice_winter_summer.png"     = "06_sea-ice_winter-vs-summer_2015-2018-2020.png",
  "fig07_faceted_temp_winter_summer.png"    = "07_air-temperature_winter-vs-summer_2015-2018-2020.png",
  "fig08_faceted_wind_winter_summer.png"    = "08_wind-speed_winter-vs-summer_2015-2018-2020.png",
  "fig09_hist_nOprD_A_all.png"              = "09_histogram_vessel-days_all-hex-months.png",
  "fig10_hist_nOprD_A_positive.png"         = "10_histogram_vessel-days_positive-only.png",
  "fig11_port_dist_distribution.png"        = "11_distance-to-port_distribution.png",
  "fig12_map_port_dist_with_ports.png"      = "12_map_distance-to-nearest-port_with-ports.png",
  "fig13_corr_matrix.png"                   = "13_correlation-matrix_env-covariates.png",
  "fig14_joint_density_ice_temp.png"        = "14_joint-density_ice-vs-air-temperature.png",
  "fig15_air_temp_resid_vs_ice.png"         = "15_orthogonalized-temperature_check.png",
  "fig16_dropped_hexes_map.png"             = "16_row-loss-audit_dropped-ice-NA-hexes.png"
)
fpath <- function(f) {
  if (f %in% names(.FIG_RENAME)) f <- unname(.FIG_RENAME[[f]])
  file.path(CONFIG$fig_dir, f)
}

# =============================================================================
## 1. PACKAGES
# =============================================================================
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(terra)
  library(exactextractr)
  library(viridis)
  library(scales)
})
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
HAS_RNE       <- requireNamespace("rnaturalearth", quietly = TRUE)
HAS_MGCV      <- requireNamespace("mgcv", quietly = TRUE)

banner <- function(txt) {
  cat("\n", strrep("=", 78), "\n## ", txt, "\n", strrep("=", 78), "\n", sep = "")
}
say <- function(...) cat(sprintf(...), "\n")

# -----------------------------------------------------------------------------
# ledger(): the running ROW / HEX / DATE / NA audit. Call after every step and
# announce what the step ADDS or DELETES and why. Nothing silent.
# -----------------------------------------------------------------------------
ledger <- function(df, step, action = "", key_vars = c("nOprD_A", "ice_conc",
                                                        "air_temp", "wind_speed")) {
  nrw <- nrow(df)
  nhx <- if ("hexID" %in% names(df)) length(unique(df$hexID)) else NA
  ndt <- if ("date"  %in% names(df)) length(unique(df$date))  else NA
  cat(sprintf("\n[LEDGER] %-46s\n", step))
  if (nzchar(action)) cat(sprintf("         action: %s\n", action))
  cat(sprintf("         rows=%s  unique_hexID=%s  unique_dates=%s\n",
              format(nrw, big.mark = ","), format(nhx, big.mark = ","), ndt))
  present <- key_vars[key_vars %in% names(df)]
  if (length(present) > 0) {
    na_str <- paste(sprintf("%s=%s", present,
                            vapply(present, function(v) format(sum(is.na(df[[v]])),
                                                               big.mark = ","),
                                   character(1))),
                    collapse = "  ")
    cat(sprintf("         NA: %s\n", na_str))
  }
  invisible(df)
}

cat("\n",
    "#############################################################\n",
    "#  Data Handling for Melting Ports — build_dataset_and_eda  #\n",
    "#  run date: ", as.character(Sys.Date()), "\n",
    "#############################################################\n", sep = "")
say("data_dir : %s", CONFIG$data_dir)
say("out_dir  : %s", CONFIG$out_dir)

# =============================================================================
## 2. SHIP — STACK + STRUCTURAL-ZERO RECODE
#
#  PROVENANCE
#    Source   : Kapsar, K., Sullender, B., Liu, J., & Poe, A. (2022b).
#               "North Pacific and Arctic Marine Vessel Traffic Dataset
#                (2015-2020)", Data in Brief 44:108531.
#    Form     : 72 monthly SpeedHex shapefiles, 6,553 hex cells x 72 months.
#    Native CRS: EPSG:3338 (NAD83 / Alaska Albers).
#    Key vars : nOprD_* (operating vessel-days by class), nMMSI_* (unique
#               vessels by class), SOG_* / SOGsd_* (mean speed & sd, kn).
#               Class suffixes: A=all, C=cargo, T=tanker, F=fishing, O=other,
#               L=... (note the dataset ships a typo column 'nOPrD_L').
#    Accessed : (fill in your download date).
#
#  This consolidated script does NOT re-stack the 72 shapefiles every run; if
#  the attribute table 'cleaned_shipping.rds' already exists it is loaded. The
#  stacking + recode logic that produced it (from ship_analysis.R) is documented
#  inline so the provenance is auditable. To re-stack from scratch, point
#  raw_speedhex_dir at the unzipped dataset and set FORCE_RESTACK <- TRUE.
# =============================================================================
banner("2. SHIP — stack + structural-zero recode")

FORCE_RESTACK <- FALSE
ship_path <- dpath("cleaned_shipping.rds")

# robust count-column matcher (also catches the 'nOPrD_L' typo)
count_re <- "^n(MMSI|OprD|OPrD)_"
sog_re   <- "^SOG"

if (file.exists(ship_path) && !FORCE_RESTACK) {
  say("Loading existing attribute table: %s", ship_path)
  ship <- readRDS(ship_path)
} else {
  stop("cleaned_shipping.rds not found. Re-stacking requires the raw SpeedHex ",
       "shapefiles (see ship_analysis.R); set raw_speedhex_dir and FORCE_RESTACK.")
}
ship <- ledger(ship, "ship table loaded", "ADD: 72 monthly SpeedHex layers stacked upstream")

# ---- 2a. STRUCTURAL-ZERO VERIFICATION (do not blindly recode) --------------
# The recode rule is: COUNT columns (nMMSI_*, nOprD_*) NA -> 0 because a missing
# count means "no vessel of that class observed in that hex-month" (a true,
# structural zero). SPEED columns (SOG_*, SOGsd_*) are LEFT as NA because the
# mean speed of zero ships is undefined (a recode to 0 would be a fabricated
# value). Before trusting that recode we verify the missingness is structural,
# not a reception gap.
count_cols <- grep(count_re, names(ship), value = TRUE)
sog_cols   <- grep(sog_re,   names(ship), value = TRUE)
say("count columns (%d): %s", length(count_cols), paste(count_cols, collapse = ", "))
say("speed columns (%d): %s", length(sog_cols),   paste(sog_cols,   collapse = ", "))

# (a) is any ENTIRE month missing? (would indicate a data-reception gap, not zeros)
month_completeness <- ship %>%
  group_by(date) %>%
  summarise(n_hex = n(), .groups = "drop")
say("\n(a) Month completeness — every month should have %d hexes:", CONFIG$n_hex_expected)
print(as.data.frame(month_completeness)[c(1:3, (nrow(month_completeness)-2):nrow(month_completeness)), ])
say("    distinct hex-counts across months: %s",
    paste(sort(unique(month_completeness$n_hex)), collapse = ", "))
say("    => %s",
    ifelse(all(month_completeness$n_hex == CONFIG$n_hex_expected),
           "OK: no month is short — no reception gap; missing counts are structural.",
           "WARNING: a month is short — investigate a possible reception gap."))

# (b) current NA state of the loaded table (counts already recoded upstream?)
n_na_counts <- sum(vapply(ship[count_cols], function(x) sum(is.na(x)), integer(1)))
n_na_sog    <- sum(vapply(ship[sog_cols],   function(x) sum(is.na(x)), integer(1)))
say("\n(b) NA in count columns (loaded): %s", format(n_na_counts, big.mark = ","))
say("    NA in speed columns (loaded): %s", format(n_na_sog,    big.mark = ","))

if (n_na_counts > 0) {
  # We still have raw NAs -> show per-month / per-class NA fractions, then recode.
  say("\n    Per-month NA fraction for nOprD_A and nMMSI_A (first/last 3):")
  na_by_month <- ship %>% group_by(date) %>%
    summarise(na_nOprD_A = mean(is.na(nOprD_A)),
              na_nMMSI_A = mean(is.na(nMMSI_A)), .groups = "drop")
  print(as.data.frame(na_by_month)[c(1:3, (nrow(na_by_month)-2):nrow(na_by_month)), ])
  say("    Per-class total NA fractions:")
  print(round(sapply(ship[count_cols], function(x) mean(is.na(x))), 4))

  n_cells_recoded <- n_na_counts
  ship[count_cols] <- lapply(ship[count_cols], function(x) { x[is.na(x)] <- 0; x })
  say("    RECODED NA->0 in %s count cells (structural zeros).",
      format(n_cells_recoded, big.mark = ","))
} else {
  say("\n    Count columns already recoded NA->0 upstream (in ship_analysis.R).")
  say("    Verified here: every month is complete and SOG columns retain their")
  say("    NAs (%s), confirming the recode was applied to COUNTS only.",
      format(n_na_sog, big.mark = ","))
}

# (c) resulting zero proportion of the primary intensity variable
zero_prop_all <- mean(ship$nOprD_A == 0)
say("\n(c) Zero proportion of nOprD_A (all vessels): %.4f", zero_prop_all)
ship <- ledger(ship, "after structural-zero handling",
               "DELETE: none (recode is value-level, not row-level)")

# =============================================================================
## 3. ICE — DOWNLOAD & EXTRACT  (loaded from cached extraction if present)
#
#  PROVENANCE
#    Product  : see CONFIG$ice_product (CDR v5, variable 'aice').
#    Package  : icecdr::cdr_arctic_monthly(date_range, variables="aice",
#               version=5).  Native CRS EPSG:3413, 25 km polar-stereographic.
#    Transform: reproject to EPSG:3338 (bilinear), then area-weighted mean per
#               hex via exactextractr::exact_extract(fun="mean") -> ice_conc in
#               [0,1].  (See ice_analysis.R for the extraction loop.)
#    Accessed : (fill in your download date).
#
#  Re-running the download needs the icecdr package + network; this script loads
#  the cached per-hex-month table 'ice_concentration.rds' if it exists.
# =============================================================================
banner("3. ICE — load cached extraction (ice_concentration.rds)")
say("Ice product pinned: %s", CONFIG$ice_product)

ice_path <- dpath("ice_concentration.rds")
if (!file.exists(ice_path))
  stop("ice_concentration.rds not found; run the icecdr extraction (ice_analysis.R).")
ice <- readRDS(ice_path)
ice$hexID <- as.character(ice$hexID)
ice <- ledger(ice, "ice table loaded",
              "ADD: monthly CDR-v5 sea-ice concentration, area-weighted per hex",
              key_vars = "ice_conc")
say("ice_conc range: [%.3f, %.3f];  NA cells: %s",
    min(ice$ice_conc, na.rm = TRUE), max(ice$ice_conc, na.rm = TRUE),
    format(sum(is.na(ice$ice_conc)), big.mark = ","))

# =============================================================================
## 4. AIR — ERA5 2m temperature + 10m wind  (extract from era5_monthly.nc)
#
#  PROVENANCE
#    Product  : ERA5 monthly averaged reanalysis, single levels:
#               '2m_temperature' (K) and '10m_wind_speed' (m/s).
#    Retrieval: Copernicus CDS. In a one-R-file world this is best done with the
#               `ecmwfr` package (wf_request); a documented Python cdsapi call
#               (download_era5.py, area = c(80,180,40,-180)) is the alternative
#               prerequisite that produced era5_monthly.nc here.
#    Transform: reproject to EPSG:3338 (bilinear); area-weighted mean per hex via
#               exactextractr; Kelvin -> Celsius for air_temp.
#    Accessed : (fill in your download date).
#
#  If the cached per-hex tables (air_temp.rds, wind_speed.rds) exist they are
#  loaded; otherwise they are extracted from era5_monthly.nc on the fly.
# =============================================================================
banner("4. AIR — ERA5 temperature + wind (per-hex extraction)")

hex_grid <- readRDS(dpath("hex_grid_template.rds"))
hex_grid$hexID <- as.character(hex_grid$hexID)

air_path  <- dpath("air_temp.rds")
wind_path <- dpath("wind_speed.rds")

if (file.exists(air_path) && file.exists(wind_path)) {
  say("Loading cached ERA5 extractions: air_temp.rds, wind_speed.rds")
  temp_df <- readRDS(air_path);  temp_df$hexID <- as.character(temp_df$hexID)
  wind_df <- readRDS(wind_path); wind_df$hexID <- as.character(wind_df$hexID)
} else if (file.exists(dpath("era5_monthly.nc"))) {
  say("Extracting ERA5 from era5_monthly.nc (no cached tables found) ...")
  era5 <- terra::rast(dpath("era5_monthly.nc"))
  nl <- terra::nlyr(era5)
  half <- nl / 2
  temp_layers <- era5[[1:half]]
  wind_layers <- era5[[(half + 1):nl]]
  temp_3338 <- terra::project(temp_layers, paste0("EPSG:", CONFIG$crs_ship), method = "bilinear")
  wind_3338 <- terra::project(wind_layers, paste0("EPSG:", CONFIG$crs_ship), method = "bilinear")
  secs <- as.numeric(terra::time(temp_layers))
  month_dates <- as.Date(as.POSIXct(secs, origin = "1970-01-01", tz = "UTC"))
  temp_df <- wind_df <- NULL
  for (m in seq_along(month_dates)) {
    vt <- exactextractr::exact_extract(temp_3338[[m]], hex_grid, fun = "mean", progress = FALSE)
    vw <- exactextractr::exact_extract(wind_3338[[m]], hex_grid, fun = "mean", progress = FALSE)
    temp_df <- rbind(temp_df, data.frame(hexID = hex_grid$hexID, date = month_dates[m], air_temp = vt - 273.15))
    wind_df <- rbind(wind_df, data.frame(hexID = hex_grid$hexID, date = month_dates[m], wind_speed = vw))
  }
} else {
  stop("Neither cached ERA5 tables nor era5_monthly.nc found.")
}
temp_df <- ledger(temp_df, "air_temp table", "ADD: ERA5 2m temperature (Celsius), per hex", "air_temp")
wind_df <- ledger(wind_df, "wind_speed table", "ADD: ERA5 10m wind speed (m/s), per hex", "wind_speed")

# =============================================================================
## 5. MERGE — ship + ice + air on (hexID, date)
# =============================================================================
banner("5. MERGE — ship + ice + air")

ship$hexID <- as.character(ship$hexID)
if (!"date" %in% names(ship))
  ship$date <- as.Date(sprintf("%04d-%02d-01", ship$year, ship$month))

df <- ship %>%
  left_join(ice,     by = c("hexID", "date")) %>%
  left_join(temp_df, by = c("hexID", "date")) %>%
  left_join(wind_df, by = c("hexID", "date"))
# Normalise date class: some upstream tables store it as POSIXct; force Date so
# downstream date filters (EDA facet maps) and the final table are clean.
df$date <- as.Date(df$date)
df <- ledger(df, "after ship+ice+air join", "ADD: ice_conc, air_temp, wind_speed columns")

# =============================================================================
## 6. ROW-LOSS LOCALIZATION — which environmental layer introduces NA hexes?
#
#  The analytic table is full at 471,816 rows, but ~10,368 rows carry an NA in
#  an environmental covariate; any later na.omit() would silently drop them. We
#  localize the cause to a layer (ICE vs AIR), count it, and save the dropped
#  hexIDs so the loss is auditable. We do NOT na.omit() silently.
# =============================================================================
banner("6. ROW-LOSS LOCALIZATION — ice vs air")

na_ice  <- is.na(df$ice_conc)
na_air  <- is.na(df$air_temp)
na_wind <- is.na(df$wind_speed)

say("rows with NA ice_conc   : %s", format(sum(na_ice),  big.mark = ","))
say("rows with NA air_temp   : %s", format(sum(na_air),  big.mark = ","))
say("rows with NA wind_speed : %s", format(sum(na_wind), big.mark = ","))
say("rows NA in ICE only     : %s", format(sum(na_ice & !na_air & !na_wind), big.mark = ","))
say("rows NA in AIR/WIND only: %s", format(sum(!na_ice & (na_air | na_wind)), big.mark = ","))

dropped <- df[na_ice | na_air | na_wind, ]
dropped_hex <- sort(unique(dropped$hexID))
n_drop_rows <- nrow(dropped)
n_drop_hex  <- length(dropped_hex)
driver <- if (sum(na_ice) >= sum(na_air | na_wind)) "ICE" else "AIR"
say("\n=> ROW LOSS DRIVER: %s layer.", driver)
say("   %s rows across %s hexes (= %d hexes x %d months) would be deleted.",
    format(n_drop_rows, big.mark = ","), n_drop_hex,
    n_drop_hex, CONFIG$n_months_expected)

# save the audit list
drop_file <- file.path(CONFIG$out_dir, "dropped_hexes.csv")
write.csv(data.frame(hexID = dropped_hex,
                     reason = sprintf("NA in %s extraction (outside raster coverage)", driver)),
          drop_file, row.names = FALSE)
say("   dropped-hex audit written: %s", drop_file)

# explicit, documented deletion (NOT silent)
df_full <- df                                   # keep the full table for the loss map
df <- df[!(na_ice | na_air | na_wind), ]
df <- ledger(df, "after explicit env-NA deletion",
             sprintf("DELETE: %s rows (%d hexes) — NA from %s layer",
                     format(n_drop_rows, big.mark = ","), n_drop_hex, driver))

# Figure 16: where are the dropped hexes? (audit map)
tryCatch({
  cent_all <- sf::st_centroid(hex_grid)
  xy <- sf::st_coordinates(cent_all)
  audit <- data.frame(hexID = hex_grid$hexID, x = xy[,1], y = xy[,2])
  audit$dropped <- audit$hexID %in% dropped_hex
  p <- ggplot(audit, aes(x, y, colour = dropped)) +
    geom_point(size = 0.4) +
    scale_colour_manual(values = c("FALSE" = "grey80", "TRUE" = "red"),
                        labels = c("kept", "dropped (NA ice)")) +
    coord_fixed() + theme_minimal() +
    labs(title = "Row-loss audit — hexes dropped for NA sea-ice",
         subtitle = sprintf("%d hexes dropped (driver: %s layer)", n_drop_hex, driver),
         x = "centroid_x (m, EPSG:3338)", y = "centroid_y (m)", colour = NULL)
  ggsave(fpath("fig16_dropped_hexes_map.png"), p, width = 8, height = 7, dpi = 110)
  say("   saved fig16_dropped_hexes_map.png")
}, error = function(e) say("   (fig16 skipped: %s)", conditionMessage(e)))

# =============================================================================
## 7. DISTANCE-TO-PORT — documented source, explicit rule, geodesic metric
#
#  Source : Natural Earth 'ports' (1:10m, cultural), via rnaturalearth::ne_download.
#           Public, documented, versioned. (World Port Index / NGA is the
#           alternative; we use Natural Earth for one-package reproducibility.)
#  Rule   : keep all ports whose lon/lat fall inside the STUDY bounding box
#           (derived from the hex grid extent, transformed to WGS84). No manual
#           hand-picking.
#  Metric : GEODESIC great-circle distance via sf::st_distance on geographic
#           coordinates (s2). Hex centroids are transformed to WGS84 and the
#           nearest-port distance is taken, then converted to km.
#  Caveat : ENDOGENEITY — ports sit where traffic already concentrates, so
#           port_dist partly encodes the response. Flagged for the reader; we do
#           NOT model it here.
# =============================================================================
banner("7. DISTANCE-TO-PORT — Natural Earth ports, geodesic")

ports_cache <- dpath("ne_ports_raw.rds")
ports_raw <- if (file.exists(ports_cache)) {
  say("Loading cached Natural Earth ports: %s", ports_cache); readRDS(ports_cache)
} else if (HAS_RNE) {
  say("Downloading Natural Earth ports (1:10m) ...")
  rnaturalearth::ne_download(scale = 10, type = "ports",
                             category = "cultural", returnclass = "sf")
} else stop("rnaturalearth not available and no cached ne_ports_raw.rds.")
say("Natural Earth ports available globally: %d", nrow(ports_raw))

# study bounding box from the hex grid, in WGS84
hex_bb_geo <- hex_grid %>% sf::st_transform(CONFIG$crs_geo) %>% sf::st_bbox()
say("Study bbox (WGS84): lon [%.1f, %.1f]  lat [%.1f, %.1f]",
    hex_bb_geo["xmin"], hex_bb_geo["xmax"], hex_bb_geo["ymin"], hex_bb_geo["ymax"])

ports_geo <- sf::st_transform(ports_raw, CONFIG$crs_geo)
pc <- sf::st_coordinates(ports_geo)
ports_geo$lon <- pc[, 1]; ports_geo$lat <- pc[, 2]

# SELECTION RULE: inside the study bbox. The grid crosses the dateline, so the
# bbox in lon can be misleading; we therefore select on latitude band AND the
# union of the two longitude wings (Western Hemisphere wing and the
# Eastern-Hemisphere >150E wing) that the Pacific-Arctic grid actually spans.
lat_ok <- ports_geo$lat >= 50 & ports_geo$lat <= 80
lon_ok <- (ports_geo$lon >= 150 & ports_geo$lon <= 180) |
          (ports_geo$lon >= -180 & ports_geo$lon <= -120)
ports_sel <- ports_geo[lat_ok & lon_ok, ]
say("Ports selected by rule (lat 50-80N AND Pacific-Arctic lon wings): %d",
    nrow(ports_sel))
say("\nFULL TABLE OF PORTS USED:")
print(as.data.frame(ports_sel)[, c("name", "lon", "lat")], row.names = FALSE)

# hex centroids in WGS84 for geodesic distance
cent_geo <- hex_grid %>% sf::st_centroid() %>% sf::st_transform(CONFIG$crs_geo)
# nearest-port geodesic distance (s2 great-circle), metres -> km
D <- sf::st_distance(cent_geo, ports_sel)             # [n_hex x n_ports], metres
nearest_idx <- max.col(-as.matrix(D))
port_dist_km <- as.numeric(D[cbind(seq_len(nrow(D)), nearest_idx)]) / 1000

# centroids in EPSG:3338 too (analysis CRS), for the final table
cent_alb <- sf::st_coordinates(hex_grid %>% sf::st_centroid())
port_tab <- data.frame(
  hexID        = as.character(hex_grid$hexID),
  centroid_x   = cent_alb[, 1],
  centroid_y   = cent_alb[, 2],
  port_dist    = port_dist_km,
  port_dist_log = log1p(port_dist_km),               # documented transform: right-skewed
  nearest_port = ports_sel$name[nearest_idx],
  stringsAsFactors = FALSE
)
say("\nport_dist (km) summary:")
print(summary(port_tab$port_dist))
say("ENDOGENEITY CAVEAT: ports tend to sit where traffic already is, so")
say("port_dist partially encodes the response. Flagged; not modelled here.")

# join distance + centroids into the analytic table
df <- df %>% left_join(port_tab, by = "hexID")
df <- ledger(df, "after port-distance join", "ADD: port_dist, port_dist_log, centroid_x/y")

# Figure 11: port_dist distribution
p11 <- ggplot(port_tab, aes(port_dist)) +
  geom_histogram(bins = 50, fill = "#1F5FBF", colour = "white") +
  labs(title = "Distribution of distance-to-nearest-port (per hex)",
       x = "port_dist (km, geodesic)", y = "n hexes") + theme_minimal()
ggsave(fpath("fig11_port_dist_distribution.png"), p11, width = 8, height = 5, dpi = 110)

# Figure 12: map of hex centroids coloured by port_dist, ports overlaid (EPSG:3338)
tryCatch({
  ports_alb <- sf::st_transform(ports_sel, CONFIG$crs_ship)
  pa <- sf::st_coordinates(ports_alb)
  p12 <- ggplot(port_tab, aes(centroid_x, centroid_y, colour = port_dist)) +
    geom_point(size = 0.5) +
    scale_colour_viridis_c(option = "C", name = "port_dist (km)") +
    geom_point(data = data.frame(x = pa[,1], y = pa[,2]),
               aes(x, y), inherit.aes = FALSE, shape = 17, size = 2, colour = "red") +
    coord_fixed() + theme_minimal() +
    labs(title = "Hex centroids coloured by distance to nearest port",
         subtitle = "red triangles = Natural Earth ports used (EPSG:3338)",
         x = "centroid_x (m)", y = "centroid_y (m)")
  ggsave(fpath("fig12_map_port_dist_with_ports.png"), p12, width = 8, height = 7, dpi = 110)
  say("saved fig11/fig12 (port distance EDA)")
}, error = function(e) say("(fig12 skipped: %s)", conditionMessage(e)))

# =============================================================================
## 8. CONCURVITY (data-level COLLINEARITY only)
#
#  We quantify pairwise Pearson correlations among ice_conc, air_temp,
#  wind_speed, port_dist and show 2D-binned joint densities for the key pairs.
#  NB: this is data-level collinearity. True 'concurvity' is a model-stage
#  diagnostic and is OUT OF SCOPE here.
# =============================================================================
banner("8. CONCURVITY — data-level collinearity")

cc_vars <- c("ice_conc", "air_temp", "wind_speed", "port_dist")
cor_mat <- cor(df[, cc_vars], use = "complete.obs", method = "pearson")
say("Pearson correlation matrix (complete obs):")
print(round(cor_mat, 3))
say("\nNote: ice_conc vs air_temp ~ %.3f (strong negative, as expected).",
    cor_mat["ice_conc", "air_temp"])

# Figure 13: correlation heatmap
cor_long <- as.data.frame(as.table(cor_mat))
names(cor_long) <- c("v1", "v2", "r")
p13 <- ggplot(cor_long, aes(v1, v2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(limits = c(-1, 1), low = "#B2182B", mid = "white",
                       high = "#2166AC", midpoint = 0) +
  labs(title = "Pairwise Pearson correlation (data-level collinearity)",
       x = NULL, y = NULL) + theme_minimal()
ggsave(fpath("fig13_corr_matrix.png"), p13, width = 6.5, height = 5.5, dpi = 110)

# Figure 14: 2D density ice vs air_temp
p14 <- ggplot(df, aes(ice_conc, air_temp)) +
  geom_bin2d(bins = 70) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "count") +
  labs(title = "Joint density: ice_conc vs air_temp",
       subtitle = sprintf("Pearson r = %.3f", cor_mat["ice_conc","air_temp"]),
       x = "ice_conc (0-1)", y = "air_temp (C)") + theme_minimal()
ggsave(fpath("fig14_joint_density_ice_temp.png"), p14, width = 8, height = 5.5, dpi = 110)
say("saved fig13 (corr matrix) and fig14 (ice-temp joint density)")

# =============================================================================
## 9. ORTHOGONALIZED-TEMPERATURE TRANSFORMATION (data construction, NOT a model)
#
#  We construct air_temp_resid = residuals of a univariate penalised-spline
#  smoother of air_temp on ice_conc. RATIONALE: ice_conc and air_temp are
#  strongly collinear (Section 8); the residual is the part of temperature NOT
#  explained by ice, so downstream modelling can later choose raw temperature OR
#  its ice-orthogonalised residual. This is a DERIVED VARIABLE, not an analysis
#  of anything. A penalised spline (mgcv) is used (vs a polynomial) because the
#  ice->temperature relationship is clearly non-linear; smoothing parameter
#  chosen by REML.
# =============================================================================
banner("9. ORTHOGONALIZED TEMPERATURE — air_temp_resid (derived variable)")

if (HAS_MGCV) {
  suppressPackageStartupMessages(library(mgcv))
  # univariate smoother used ONLY to build a residual covariate (bam for speed)
  temp_helper <- bam(air_temp ~ s(ice_conc, bs = "ps", k = 20, m = c(3, 2)),
                     data = df, method = "fREML", discrete = TRUE, nthreads = 4)
  df$air_temp_resid <- as.numeric(residuals(temp_helper))
  var_expl <- summary(temp_helper)$dev.expl
  say("Fraction of air_temp variance explained by ice_conc (smoother): %.4f", var_expl)
  say("cor(air_temp_resid, ice_conc) after orthogonalisation: %.5f",
      cor(df$air_temp_resid, df$ice_conc, use = "complete.obs"))

  # Figure 15: residual vs ice (should be a flat cloud)
  set.seed(1)
  samp <- df[sample(nrow(df), min(40000, nrow(df))), ]
  p15 <- ggplot(samp, aes(ice_conc, air_temp_resid)) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 10),
                colour = "firebrick", se = FALSE) +
    labs(title = "air_temp_resid vs ice_conc (orthogonalisation check)",
         subtitle = sprintf("ice explains %.1f%% of air_temp; residual correlation ~ 0",
                            100 * var_expl),
         x = "ice_conc (0-1)", y = "air_temp_resid (C)") + theme_minimal()
  ggsave(fpath("fig15_air_temp_resid_vs_ice.png"), p15, width = 8, height = 5.5, dpi = 110)
  say("saved fig15 (orthogonalisation check)")
} else {
  say("mgcv not available — falling back to a 3rd-order polynomial for orthogonalisation.")
  temp_helper <- lm(air_temp ~ poly(ice_conc, 3), data = df)
  df$air_temp_resid <- as.numeric(residuals(temp_helper))
  var_expl <- summary(temp_helper)$r.squared
  say("R^2 of air_temp ~ poly(ice_conc,3): %.4f", var_expl)
}

# =============================================================================
## 10. DERIVED VARIABLES + FINAL TIDY TABLE
# =============================================================================
banner("10. DERIVED VARIABLES + final analytic table")

df <- df %>%
  mutate(
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    # ice_free_flag: conventional open-water indicator (ice_conc < 0.15)
    ice_free_flag = as.integer(ice_conc < CONFIG$ice_free_threshold),
    # post-2017 indicator — labelled honestly as a time break, NOT 'the Polar Code effect'
    post2017 = as.integer(date >= CONFIG$polar_code_date)
  )
say("Derived columns added:")
say("  ice_free_flag = 1{ ice_conc < %.2f }  (open-water indicator)", CONFIG$ice_free_threshold)
say("  post2017      = 1{ date >= %s }  (time break; NOT causal 'Polar Code effect')",
    CONFIG$polar_code_date)
say("  air_temp_resid = ice-orthogonalised temperature (Section 9)")
say("  centroid_x/y   = hex centroid in EPSG:3338")

keep <- c("hexID", "date", "year", "month",
          "nOprD_A", count_cols,
          "ice_conc", "air_temp", "wind_speed", "air_temp_resid",
          "port_dist", "port_dist_log", "nearest_port",
          "centroid_x", "centroid_y", "ice_free_flag", "post2017")
keep <- unique(keep[keep %in% names(df)])
analytic <- df[, keep]
analytic <- ledger(analytic, "FINAL analytic table", "ADD: derived flags; SELECT tidy columns")

# The final tidy table is the core DATA product of this stage, so it lives in
# Data/ (next to the inputs), while logs/audit/figures live in Data/output/.
out_rds <- file.path(CONFIG$data_dir, "analytic_dataset_built.rds")
saveRDS(analytic, out_rds)
say("Saved final analytic table: %s", out_rds)

# =============================================================================
## 11. EXPLORATION (figures for the LaTeX reader)
# =============================================================================
banner("11. EDA FIGURES")

# --- Fig 1: KEY CURVE — zero-proportion of nOprD_A by ice-concentration bin ---
zb <- analytic %>%
  mutate(ice_bin = cut(ice_conc, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(ice_bin) %>%
  summarise(zero_prop = mean(nOprD_A == 0), n = n(), .groups = "drop")
say("Zero-proportion of nOprD_A by ice bin:")
print(as.data.frame(zb))
p1 <- ggplot(zb[!is.na(zb$ice_bin), ], aes(ice_bin, zero_prop)) +
  geom_col(fill = "#0B2E6D") +
  labs(title = "Zero-proportion of nOprD_A by sea-ice bin",
       x = "ice concentration bin", y = "P(nOprD_A = 0)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(fpath("fig01_zero_prop_by_ice_bin.png"), p1, width = 8, height = 5, dpi = 110)

# --- Fig 2: monthly total traffic time series ---
mt <- analytic %>% group_by(date) %>%
  summarise(total = sum(nOprD_A), .groups = "drop")
p2 <- ggplot(mt, aes(date, total)) + geom_line(colour = "#1F5FBF") +
  labs(title = "Monthly total vessel-days (nOprD_A)", x = NULL, y = "sum nOprD_A") +
  theme_minimal()
ggsave(fpath("fig02_monthly_traffic_timeseries.png"), p2, width = 9, height = 4.5, dpi = 110)

# --- Fig 3 & 4: seasonal cycle + composition by vessel type -------------------
# class columns present (cargo/tanker/fishing/other)
cls <- c(C = "nOprD_C", T_ = "nOprD_T", O = "nOprD_O", F_ = "nOprD_F")
cls <- cls[cls %in% names(analytic)]
if (length(cls) >= 2) {
  comp <- analytic %>%
    summarise(across(all_of(unname(cls)), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "col", values_to = "total") %>%
    mutate(vessel = recode(col, nOprD_C = "Cargo", nOprD_T = "Tanker",
                           nOprD_O = "Other", nOprD_F = "Fishing"),
           pct = total / sum(total))
  say("Vessel-type composition (sum vessel-days):")
  print(as.data.frame(comp[, c("vessel", "total", "pct")]), row.names = FALSE)

  seas <- analytic %>%
    group_by(date) %>%
    summarise(across(all_of(unname(cls)), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-date, names_to = "col", values_to = "total") %>%
    mutate(vessel = recode(col, nOprD_C = "Cargo", nOprD_T = "Tanker",
                           nOprD_O = "Other", nOprD_F = "Fishing"))
  p3 <- ggplot(seas, aes(date, total, colour = vessel)) +
    geom_line() + facet_wrap(~ vessel, scales = "free_y") +
    scale_colour_viridis_d(option = "D", guide = "none") +
    labs(title = "Seasonal cycle by vessel type", x = NULL, y = "monthly vessel-days") +
    theme_minimal()
  ggsave(fpath("fig03_seasonal_cycle_by_vessel_type.png"), p3, width = 9, height = 6, dpi = 110)

  p4 <- ggplot(comp, aes(x = 2, y = total, fill = vessel)) +
    geom_col(width = 1, colour = "white") + coord_polar("y") + xlim(0.5, 2.5) +
    geom_text(aes(label = percent(pct, accuracy = 0.1)),
              position = position_stack(vjust = 0.5), colour = "white", size = 4) +
    scale_fill_viridis_d(option = "D") + theme_void() +
    labs(title = "Vessel-type composition (2015-2020)", fill = NULL)
  ggsave(fpath("fig04_vessel_type_composition.png"), p4, width = 6.5, height = 5.5, dpi = 110)
}

# --- Fig 5: raw-raster vs hex-aggregated ice (one month) ----------------------
# We only have the extracted ice table (no raw ice NetCDF cached here), so we
# render the hex-aggregated panel; the raw-raster panel needs ice_data/*.nc.
tryCatch({
  jan <- analytic %>% filter(date == as.Date("2015-01-01")) %>% select(hexID, ice_conc)
  jan_hex <- hex_grid %>% left_join(jan, by = "hexID")
  p5 <- ggplot(jan_hex) +
    geom_sf(aes(fill = ice_conc), colour = NA) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey85", limits = c(0,1),
                         name = "ice (0-1)") +
    theme_void() + labs(title = "Hex-aggregated sea-ice — Jan 2015",
                        subtitle = "(raw 25-km raster panel requires ice_data/*.nc)")
  ggsave(fpath("fig05_ice_hex_jan2015.png"), p5, width = 7, height = 7, dpi = 110)
}, error = function(e) say("(fig05 skipped: %s)", conditionMessage(e)))

# --- Fig 6/7/8: faceted env maps, winter vs summer, three years ---------------
make_facet_map <- function(var, title, opt, fname) {
  tryCatch({
    sel <- expand.grid(year = c(2015, 2018, 2020), mm = c("01", "07"),
                       stringsAsFactors = FALSE) %>%
      mutate(date = as.Date(sprintf("%d-%s-01", year, mm)),
             season = ifelse(mm == "01", "Winter (Jan)", "Summer (Jul)"))
    md <- analytic %>% filter(date %in% sel$date) %>%
      select(hexID, date, all_of(var)) %>%
      left_join(sel, by = "date")
    gg <- hex_grid %>% left_join(md, by = "hexID") %>% filter(!is.na(date))
    p <- ggplot(gg) + geom_sf(aes(fill = .data[[var]]), colour = NA) +
      scale_fill_viridis_c(option = opt, na.value = "grey85", name = var) +
      facet_grid(rows = vars(season), cols = vars(year)) +
      theme_void() + labs(title = title)
    ggsave(fpath(fname), p, width = 11, height = 7, dpi = 110)
  }, error = function(e) say("(%s skipped: %s)", fname, conditionMessage(e)))
}
make_facet_map("ice_conc",   "Sea-ice concentration — winter vs summer", "plasma",
               "fig06_faceted_ice_winter_summer.png")
make_facet_map("air_temp",   "2m air temperature — winter vs summer",    "inferno",
               "fig07_faceted_temp_winter_summer.png")
make_facet_map("wind_speed", "10m wind speed — winter vs summer",        "viridis",
               "fig08_faceted_wind_winter_summer.png")

# --- Fig 9/10: histograms of nOprD_A (all, and conditional on > 0) ------------
p9 <- ggplot(analytic, aes(nOprD_A)) +
  geom_histogram(bins = 60, fill = "#0B2E6D") + scale_y_log10() +
  labs(title = "nOprD_A — all hex-months (y log10)", x = "vessel-days", y = "count") +
  theme_minimal()
ggsave(fpath("fig09_hist_nOprD_A_all.png"), p9, width = 8, height = 5, dpi = 110)

p10 <- ggplot(filter(analytic, nOprD_A > 0), aes(nOprD_A)) +
  geom_histogram(bins = 60, fill = "#1F5FBF") + scale_x_log10() +
  labs(title = "nOprD_A | nOprD_A > 0 (x log10)", x = "vessel-days", y = "count") +
  theme_minimal()
ggsave(fpath("fig10_hist_nOprD_A_positive.png"), p10, width = 8, height = 5, dpi = 110)

say("EDA figures written to %s", CONFIG$fig_dir)

# =============================================================================
## 12. SUMMARY + SESSION INFO
# =============================================================================
banner("12. SUMMARY")

say("Final analytic rows         : %s", format(nrow(analytic), big.mark = ","))
say("Final unique hexes          : %s", format(length(unique(analytic$hexID)), big.mark = ","))
say("Final unique months         : %d", length(unique(analytic$date)))
say("Rows deleted (env NA)       : %s  (driver: %s layer, %d hexes)",
    format(n_drop_rows, big.mark = ","), driver, n_drop_hex)
say("Zero proportion of nOprD_A  : %.4f", mean(analytic$nOprD_A == 0))
say("cor(ice_conc, air_temp)     : %.3f", cor_mat["ice_conc", "air_temp"])
say("air_temp var expl by ice    : %.4f", var_expl)
say("port_dist (km) median       : %.1f", median(analytic$port_dist))
say("Ice product                 : %s", CONFIG$ice_product)

banner("sessionInfo() (key packages)")
si <- sessionInfo()
say("R: %s", si$R.version$version.string)
for (p in c("sf","terra","exactextractr","dplyr","ggplot2","mgcv","rnaturalearth"))
  if (requireNamespace(p, quietly = TRUE))
    say("  %-14s %s", p, as.character(packageVersion(p)))

cat("\n[build_dataset_and_eda] DONE.\n")

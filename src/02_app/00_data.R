#### Collate ECAD climate data ####

# Update to group by other seasons (Spring/Summer/Autumn/Winter) (done)
# Calculate SPI over these as well (done)
# TODO Investigate whether we should remove 0s before calculating SPI or not
# TODO Investigate NAs in `frank` due to season year being 1959! Should we
# allow this?

# TODO Write why we prefer non-airport/urban areas
# - airports and urban areas are more susceptible to increased building
# congestion over itme.

# TODO Add Portugal as well maybe??

# Keep only continental sites (and remove islands from plots) (done)
# Where there are multiple stations at the same location, keep the one with the longest record (done)
# Keep sites with good data records (done)
# Maybe look at how many sites have complete cases over a given
# period? Rather than having over 80% of observations (can do this weekly..) (done)

#### Libs ####

library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(lubridate)
library(parallel)
library(data.table)
library(sf)
library(ggplot2)
library(patchwork)
library(SPEI)

#### Functions ####

# function to convert DMS to decimal degrees
dms_to_decimal <- function(x) {
  sign <- ifelse(substr(x, 1, 1) == "-", -1, 1)

  x <- sub("^[+-]", "", x)

  deg <- as.numeric(sub(":.*", "", x))
  min <- as.numeric(substr(x, regexpr(":", x) + 1, regexpr(":", x) + 2))
  sec <- as.numeric(substr(x, nchar(x) - 1, nchar(x)))

  sign * (deg + min / 60 + sec / 3600)
}

# function to load and cleanup data
load_data <- \(files, spec_col, min_date = NULL, max_date = NULL) {
  # longer names for the columns we want to keep
  spec_col_name <- switch(spec_col,
    tx = "temp",
    fg = "wind_speed",
    rr = "rain"
  )
  q_col <- paste0("q_", spec_col)

  # load files
  data_lst <- mclapply(
    files, read.csv,
    skip = 20, header = TRUE
  )

  # Join files and cleanup data
  ret <- bind_rows(data_lst) |>
    janitor::clean_names() |>
    # match station id to station name
    left_join(select(stations, -c(hght, cn))) |>
    # remove missing or "suspect" data
    filter(!!sym(q_col) == 0) |>
    select(-c(souid, !!sym(q_col))) |>
    rename(
      station_id = staid,
      station_name = staname,
      # temp = tx
      !!spec_col_name := !!sym(spec_col)
    ) |>
    mutate(
      # change from all caps to capitalised words
      station_name = stringr::str_to_title(station_name),
      # chop off whitespace at end of station name
      station_name = stringr::str_trim(station_name),
      # convert date to date format
      # date = as.Date(date, format = "%Y%m%d")
      date = as.Date(as.character(date), format = "%Y%m%d")
    )

  if (!is.null(min_date)) {
    ret <- filter(ret, date >= min_date)
  }
  if (!is.null(max_date)) {
    ret <- filter(ret, date <= max_date)
  }

  message("Loaded and cleaned data for ", spec_col_name)
  return(ret)
}


#### Metadata & Shapefile ####

# load map of Spain to use as background in plots
# areas <- mapSpain::esp_get_munic_siane(moveCAN = TRUE) |>
areas <- sf::read_sf("data/02_app/spain_shapefile.geojson") |>
  filter(!ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla"))

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# dates over which to pull data
start_date <- as.Date("1960-01-01")
# start_date <- as.Date("1950-01-01")
end_date <- as.Date("2024-12-31")

# station file
stations_file <- "data/02_app/ECA_blend_tx/stations.txt"

# files with temperature, rain and wind speed data
temp_files <- list.files("data/02_app/ECA_blend_tx", full.names = TRUE)
rain_files <- list.files("data/02_app/ECA_blend_rr", full.names = TRUE)
ws_files <- list.files("data/02_app/ECA_blend_fg", full.names = TRUE)

# only interested in files matching these patterns
temp_files <- temp_files[grepl("TX", temp_files)]
rain_files <- rain_files[grepl("RR", rain_files)]
ws_files <- ws_files[grepl("FG", ws_files)]

# combine into a list for easier handling
files <- list(
  "temp" = temp_files,
  "rain" = rain_files,
  "ws"   = ws_files
)


#### Stations Data ####

# read stations data
stations <- read.csv(stations_file, header = TRUE, skip = 16) |>
  janitor::clean_names() |>
  mutate(cn = stringr::str_trim(cn))

# only interested in stations in Spain
stations <- filter(stations, cn == "ES")

# pull Spanish sites
st_id <- stations$staid
# prepend with 0s until of length 6
st_id <- sprintf("%06d", st_id)

# keep files from Spain
files_es <- lapply(files, \(x) {
  x[substr(x, nchar(x) - 9, nchar(x) - 4) %in% st_id]
})


#### Load Data and filter ####

# read.csv(temp_files[1], skip = 19, header = TRUE) |>
#   head()

# load clean data
# temp_data <- load_data(temp_files_es, "tx")

# load clean data for all variables
data_lst <- mapply(
  load_data,
  files_es,
  c("tx", "rr", "fg"),
  # min_date = start_date,
  min_date = start_date - years(1), # need year before for Winter data
  max_date = end_date,
  SIMPLIFY = FALSE
)

# join data
ecad_data <- reduce(
  data_lst,
  full_join,
  by = c("station_id", "station_name", "date", "lat", "lon")
) |>
  relocate(temp, .after = lon)

# dplyr too slow for this, so switch to data.table for faster processing
setDT(ecad_data)

# pull number of dates from start to finish
n_dates <- length(seq.Date(start_date, end_date, by = "day"))
n_weeks <- length(seq.Date(start_date, end_date, by = "week"))

# filter so that only stations with >80% of records are kept
station_coverage <- ecad_data[
  ,
  .(
    # n = .N,
    n_valid = sum(!is.na(temp)) # ,
    # prop_valid = mean(!is.na(temp))
    # prop_valid = n_valid / n_dates
  ),
  by = station_id
][
  ,
  prop_valid := n_valid / n_dates
][order(-prop_valid)]

ecad_filtered <- ecad_data[
  station_coverage[prop_valid > 0.8],
  on = "station_id"
]

# convert lon and lat to numeric
ecad_filtered[, `:=`(
  lon = dms_to_decimal(lon),
  lat = dms_to_decimal(lat)
)]

# Only keep continental sites
ecad_filtered <- ecad_filtered[
  !(
    # Canary Islands
    lon >= -18.5 & lon <= -13 &
      lat >= 27 & lat <= 30
  ) &
    !(
      # Balearic Islands
      lon >= 1 & lon <= 5 &
        lat >= 38 & lat <= 41
    )
  # remove the station Melilla (by name)
][
  !station_name %in% c("Melilla")
]

stopifnot(length(unique(ecad_filtered$station_id)) > 80) # ~86


#### Convert to weekly scale ####

# take weekly summaries of environmental variables
ecad_weekly <- ecad_filtered[
  ,
  .(
    # temp       = max(temp, na.rm = TRUE),
    temp_max   = as.numeric(max(temp, na.rm = TRUE)),
    # use high quantile rather than max to better reflect how hot a week is
    temp       = quantile(temp, probs = 0.9, na.rm = TRUE),
    rain       = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE)
  ),
  by = .(
    date = floor_date(date, "week"),
    station_id,
    station_name,
    lon,
    lat
  )
]

# replace NaN with NA in weekly data
ecad_weekly[is.nan(wind_speed), wind_speed := NA_real_]
# replace infs with NA
ecad_weekly[is.infinite(temp), temp := NA_real_]

# how many sites have full records (i.e all of n_weeks)?
# start_date_temp <- as.Date("1950-01-01")
# start_date_temp <- as.Date("1960-01-01")
start_date_temp <- start_date
all_weeks <- seq.Date(start_date_temp, end_date, by = "week")
n_weeks <- length(all_weeks)
complete_stations <- ecad_weekly[date >= start_date_temp][
  ,
  .(n_weeks_present = uniqueN(date)),
  by = .(station_id, station_name)
][n_weeks_present == n_weeks]
nrow(complete_stations)
# for start_date = "1950-01-01", n_sites = 54
# for start_date = "1950-01-01", n_sites = 71, maps looks more filled out!

# only keep these!
ecad_weekly_complete <- ecad_weekly[
  complete_stations,
  on = .(station_id, station_name)
]

# remove days with no rain, as they will mess up SPI estimates below
# TODO Only do when interested in modelling rain!
# ecad_weekly_complete <- ecad_weekly_complete |>
#   filter(rain > 0)


#### Calculate SPI ####

# # ecad_weekly_complete[
# #   ,
# #   season := fifelse(
# #     month(date) %in% c(10:12, 1:3),
# #     "Winter",
# #     "Summer"
# #   )
# # ]
#
# # add season
# ecad_weekly_complete[
#   ,
#   month := month(date)
# ]
#
# # TODO Rewrite in data.table
# ecad_weekly_complete <- ecad_weekly_complete |>
#   mutate(season = case_when(
#     month %in% c(12, 1:2) ~ "Winter",
#     month %in% 3:5 ~ "Spring",
#     month %in% 6:8 ~ "Summer",
#     TRUE ~ "Autumn"
#   ))
#
# # assign season-year:
# # Jan-Mar belong to the winter that started in the previous calendar year
# # ecad_weekly_complete[
# #   ,
# #   season_year := fifelse(
# #     season == "Winter" & month(date) %in% c(1, 2, 3),
# #     year(date) - 1L,
# #     year(date)
# #   )
# # ]
# # Jan-Feb belong to the winter that started in the previous calendar year
# ecad_weekly_complete[
#   ,
#   season_year := fifelse(
#     season == "Winter" & month(date) %in% c(1, 2),
#     year(date) - 1L,
#     year(date)
#   )
# ]
#
# # seasonal precipitation totals
# seasonal_rain <- ecad_weekly_complete[
#   ,
#   .(
#     rain_season = sum(rain, na.rm = TRUE)
#   ),
#   by = .(
#     station_id,
#     station_name,
#     lon,
#     lat,
#     season,
#     season_year
#   )
# ]
#
# # remove edge season at start
# seasonal_rain <- seasonal_rain[
#   !(season == "Winter" & season_year == year(start_date) - 1L)
# ]
#
# # # empirical drought index by station and season
# # seasonal_rain[
# #   ,
# #   rain_u := frank(rain_season, ties.method = "average") / (.N + 1),
# #   # by = .(station_id, season)
# #   by = .(season)
# # ]
#
# # empirical drought indexes:
# # # Calculate empirical/non-parametric seasonal SPI
# # Local: relative to each station's own seasonal climatology
# seasonal_rain[
#   ,
#   drought_local := 1 - frank(rain_season, ties.method = "average") / (.N + 1),
#   by = .(station_id, season)
# ]
#
# # Global: relative to all stations in the same season/year
# seasonal_rain[
#   ,
#   drought_global := 1 - frank(rain_season, ties.method = "average") / (.N + 1),
#   by = .(season, season_year)
# ]
#
# # standardise to normal distribution (allows for fitting normal GAM later)
# seasonal_rain <- seasonal_rain |>
#   mutate(
#     drought_local_norm = qnorm(drought_local),
#     drought_global_norm = qnorm(drought_global)
#   )
#
# # seasonal_rain[
# #   ,
# #   spi_emp := qnorm(rain_u)
# # ]
# #
# # # large positive is equivalent to dryness!
# # seasonal_rain[
# #   ,
# #   drought_spi := -spi_emp
# # ]
#
# # remove unneeded columns
# seasonal_rain[, c("rain_u", "spi_emp") := NULL]
#
# # join seasonal drought back to weekly data
# ecad_weekly_complete <- seasonal_rain[
#   ecad_weekly_complete,
#   # on = .(station_id, station_name, lon, lat, season, season_year)
#   on = .(station_id, station_name, lon, lat, season, season_year)
# ]

# make sure every station has every week
station_grid <- CJ(
  station_id = unique(ecad_weekly_complete$station_id),
  date = all_weeks
)

ecad_weekly_complete <- station_grid[
  ecad_weekly_complete,
  on = .(station_id, date)
]

# wide matrix: rows = weeks, columns = stations
rain_wide <- dcast(
  ecad_weekly_complete,
  date ~ station_id,
  value.var = "rain"
)

rain_dates <- rain_wide$date
rain_mat <- as.matrix(rain_wide[, -"date"])

# 13-week rolling SPI = approximately seasonal / 3-month SPI
spi_13 <- SPEI::spi(
  rain_mat,
  # 13 weeks, week 13 uses weeks 1-13, week 14 uses 2-14
  scale = 13,
  # assume rainfall follows Gamma
  # (automatically handles 0s)
  distribution = "Gamma",
  # Unbiased Probability Weighted Moments
  fit = "ub-pwm",
  # equal weight to every week
  kernel = list(type = "rectangular", shift = 0),
  na.rm = TRUE
)

# TODO Look into -Infs (23)
spi_13_fitted <- spi_13$fitted
as.data.frame(spi_13_fitted) |>
  filter(if_any(everything(), \(x) is.infinite(x))) |>
  nrow()
# 243 out of 3393 rows
# # give them a value of -4/4 (equivalent to 1 in 10,000 chance of rainfall)
# spi_13_fitted[spi_13_fitted == -Inf] <- -4
# spi_13_fitted[spi_13_fitted == Inf] <- 4
# change to NA, may cause problems in tails with Gaussian GAM!!
spi_13_fitted[spi_13_fitted == -Inf] <- NA
spi_13_fitted[spi_13_fitted == Inf] <- NA


# extract fitted SPI values
spi_13_wide <- as.data.table(spi_13_fitted)
spi_13_wide[, date := rain_dates]

# convert back to long format
spi_13_long <- melt(
  spi_13_wide,
  id.vars = "date",
  variable.name = "station_id",
  value.name = "spi_13w"
) |>
  # eww, mixing data.table and dplyr syntax here...
  rename(drought_local = spi_13w) |>
  mutate(station_id = as.integer(as.character(station_id)))

spi_13_long |>
  filter(is.na(drought_local)) |>
  nrow()
# 852 of 240903 rows

# join back onto weekly data
ecad_weekly_complete <- merge(
  ecad_weekly_complete,
  spi_13_long,
  by = c("station_id", "date"),
  all.x = TRUE
)


#### Remove neighbours ####

stations <- unique(
  ecad_weekly_complete[, .(station_id, station_name, lon, lat)]
)

stations_sf <- st_as_sf(
  stations,
  coords = c("lon", "lat"),
  # crs = 4326,
  crs    = st_crs(areas),
  remove = FALSE
)

stations_proj <- st_transform(stations_sf, 25830)

# find points within 25km of each other
nb <- st_is_within_distance(
  stations_proj,
  stations_proj,
  # dist = 25000
  dist = 25000
)

# create cluster IDs
edges <- rbindlist(lapply(seq_along(nb), function(i) {
  data.table(from = i, to = as.integer(nb[[i]]))
}))

g <- igraph::graph_from_data_frame(edges, directed = FALSE)

stations[, cluster_id := igraph::components(g)$membership]

stations[, is_airport := grepl(
  "Aeropuerto|Airport|Base Aerea|Aerodromo|Foronda|Cuatro Vientos",
  station_name,
  ignore.case = TRUE
)]

stations_ranked <- stations[
  ,
  .SD[order(is_airport, station_name)],
  by = cluster_id
]

stations_keep <- stations_ranked[
  ,
  .SD[1],
  by = cluster_id
]

ecad_selected <- ecad_weekly_complete[
  stations_keep,
  on = .(station_id, station_name)
]
# now left with 40 stations, with good coverage across Spain!
# Is this too many to cluster/use for changepoint detection?

# calculate season
ecad_selected <- ecad_selected |>
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12, 1:2) ~ "Winter",
      month %in% 3:5 ~ "Spring",
      month %in% 6:8 ~ "Summer",
      TRUE ~ "Autumn"
    ),
    season_year = case_when(
      season == "Winter" & month %in% c(1, 2) ~ year(date) - 1L,
      TRUE ~ year(date)
    )
  )

#### Save ####

# final changes
ecad_selected <- ecad_selected |>
  mutate(
    temp_max = temp_max / 10,
    temp = temp / 10,
    rain = rain / 10
  ) |>
  select(
    date,
    month,
    season,
    season_year,
    station_id,
    station_name,
    lon,
    lat,
    temp,
    temp_max,
    rain,
    wind_speed,
    # rain_season,
    drought_local
  )

# save
readr::write_csv(
  ecad_selected,
  "data/02_app/ecad_clean.csv.gz"
)

# Also save daily data (only keeping stations of interest)
ecad_daily <- semi_join(
  ecad_filtered,
  distinct(ecad_selected, station_id)
) |>
  mutate(
    temp = temp / 10,
    rain = rain / 10
  )
readr::write_csv(
  ecad_daily,
  "data/02_app/ecad_clean_daily.csv.gz"
)


#### Visualise ####

# convert weekly to sf
ecad_sf <- st_as_sf(
  # ecad_filtered,
  # ecad_weekly_complete,
  ecad_selected,
  coords = c("lon", "lat"),
  crs = st_crs(areas)
)

# plot a week of data (for week starting ???)
# ecad_sf |>
#   filter(date == "1998-09-13") |>
ecad_sf |>
  group_by(station_id) |>
  filter(if_all(everything(), ~ !is.na(.))) |>
  ungroup() |>
  filter(date == .data$date[1]) |>
  # pivot_longer(cols = c(temp, rain, wind_speed), names_to = "variable", values_to = "value") |>
  ggplot() +
  geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
  # geom_sf(aes(fill = temp), shape = 21, size = 8, colour = "black") +
  geom_sf(aes(fill = drought_local), shape = 21, size = 8, colour = "black") +
  # facet_wrap(~variable) +
  scale_colour_viridis_c() +
  # temperature in tenths of celcius
  # labs(fill = "°C") +
  labs(fill = "SPI") +
  scale_fill_viridis_c() +
  theme_bw()

# Plot each variable with patchwork
# vars <- c("temp", "rain", "wind_speed", "drought_global", "drought_local")
vars <- c("rain", "drought_global", "drought_local", "temp", "wind_speed")
# var_labs <- c("°C", "mm", "m/s", "Global drought index", "Local drought index")
var_labs <- c("mm", "Global drought index", "Local drought index", "°C", "m/s")


p_lst <- lapply(seq_along(vars), \(i) {
  ecad_sf |>
    filter(date == "1998-09-06") |>
    # pivot_longer(cols = c(temp, rain, wind_speed), names_to = "variable", values_to = "value") |>
    ggplot() +
    geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
    # geom_sf(aes(fill = temp), shape = 21, size = 5, colour = "black") +
    geom_sf(aes(fill = .data[[vars[[i]]]]), shape = 21, size = 5, colour = "black") +
    # facet_wrap(~variable) +
    scale_colour_viridis_c() +
    # temperature in tenths of celcius
    labs(fill = var_labs[[i]]) +
    scale_fill_viridis_c() +
    CeCl::cecl_theme(nejm_pal = FALSE) +
    theme(legend.text = element_text(angle = 45, hjust = 1))
})

wrap_plots(p_lst, nrow = 2)
# Interesting patterns here, particularly South/West for temp and rain

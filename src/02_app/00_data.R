#### Collate ECAD climate data ####

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
areas <- mapSpain::esp_get_munic_siane(moveCAN = TRUE) |>
  filter(!ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla"))

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# dates over which to pull data
start_date <- as.Date("1960-01-01")
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
  min_date = start_date,
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


#### Convert to weekly scale ####

# take weekly summaries of environmental variables
ecad_weekly <- ecad_filtered[
  ,
  .(
    temp       = max(temp, na.rm = TRUE),
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


#### Calculate SPI ####

# add season
ecad_weekly_complete[
  ,
  season := fifelse(
    month(date) %in% c(10:12, 1:3),
    "Winter",
    "Summer"
  )
]

# assign season-year:
# Jan-Mar belong to the winter that started in the previous calendar year
ecad_weekly_complete[
  ,
  season_year := fifelse(
    season == "Winter" & month(date) %in% c(1, 2, 3),
    year(date) - 1L,
    year(date)
  )
]

# seasonal precipitation totals
seasonal_rain <- ecad_weekly_complete[
  ,
  .(
    rain_season = sum(rain, na.rm = TRUE)
  ),
  by = .(
    station_id,
    station_name,
    lon,
    lat,
    season,
    season_year
  )
]

# remove edge season at start
seasonal_rain <- seasonal_rain[
  !(season == "Winter" & season_year == year(start_date) - 1L)
]

# # empirical drought index by station and season
# seasonal_rain[
#   ,
#   rain_u := frank(rain_season, ties.method = "average") / (.N + 1),
#   # by = .(station_id, season)
#   by = .(season)
# ]

# empirical drought indexes:
# Local: relative to each station's own seasonal climatology
seasonal_rain[
  ,
  drought_local := 1 - frank(rain_season, ties.method = "average") / (.N + 1),
  by = .(station_id, season)
]

# Global: relative to all stations in the same season/year
seasonal_rain[
  ,
  drought_global := 1 - frank(rain_season, ties.method = "average") / (.N + 1),
  by = .(season, season_year)
]

# # Calculate empirical/non-parametric seasonal SPI
# seasonal_rain[
#   ,
#   spi_emp := qnorm(rain_u)
# ]
#
# # large positive is equivalent to dryness!
# seasonal_rain[
#   ,
#   drought_spi := -spi_emp
# ]

# remove unneeded columns
seasonal_rain[, c("rain_u", "spi_emp") := NULL]

# join seasonal drought back to weekly data
ecad_weekly_complete <- seasonal_rain[
  ecad_weekly_complete,
  # on = .(station_id, station_name, lon, lat, season, season_year)
  on = .(station_id, station_name, lon, lat, season, season_year)
]


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

# final changes
ecad_selected <- ecad_selected |>
  mutate(
    temp = temp / 10,
    rain = rain / 10
  ) |>
  select(
    date,
    station_id,
    station_name,
    lon,
    lat,
    temp,
    rain,
    wind_speed,
    season,
    season_year,
    rain_season,
    drought_local,
    drought_global
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
# TODO Plot for both complete and non-complete cases
ecad_sf |>
  filter(date == "1998-09-06") |>
  # pivot_longer(cols = c(temp, rain, wind_speed), names_to = "variable", values_to = "value") |>
  ggplot() +
  geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
  geom_sf(aes(fill = temp), shape = 21, size = 8, colour = "black") +
  # facet_wrap(~variable) +
  scale_colour_viridis_c() +
  # temperature in tenths of celcius
  labs(fill = "°C") +
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

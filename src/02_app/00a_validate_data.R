#### Validate Data ####

# Check very large precipitation readings from ECAD observational data
# to ensure they match with AEMET OpenData set from Wikipedia

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
library(climaemet)
devtools::load_all("../CeCl/")

#### Metadata ####

Sys.getenv("AEMET_API_KEY")
# aemet_api_key("AEMET_API_KEY")
aemet_api_key(Sys.getenv("AEMET_API_KEY"))
Sys.getenv("AEMET_API_KEY")


#### Data ####


# TODO load daily data
# TODO Plot time series of days with > 100mm of rain
data <- readr::read_csv("data/02_app/ecad_clean_daily.csv.gz")

#### Calculate monthly data ####

data_month <- data |>
  mutate(month = floor_date(date, "month")) |>
  group_by(station_name, month) |>
  summarise(
    rain_sum = sum(rain, na.rm = TRUE),
    rain_mean = mean(rain, na.rm = TRUE),
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )


#### Plot time series ####

# plot days with lots of rain to identify sites of interest
data |>
  filter(
    rain > 150,
    # date >= "1980-01-01",
    # date < "2011-01-01"
  ) |>
  arrange(desc(rain)) |>
  mutate(station_name = factor(station_name)) |>
  ggplot(aes(x = date, y = rain, colour = station_name)) +
  geom_point(size = 5) +
  # guides(colour = "none") +
  cecl_theme(nejm_pal = FALSE) +
  NULL

# particular sites of interest (with days of > 200mm rain):
# - Badajoz Aeropuerto (flooded in 1989),
# - Allacant/Alicante, (flooded 1997)
# - Bilbao Aeropuerto, (flooded 1983)
# - Santioago De Compostela Aeropuerto
# San Javier Aeropuerto
# San Javioer Aeropuerto
# So days with large floods are accounted for!

# do the same for temperature
data |>
  filter(
    temp > 45
  ) |>
  arrange(desc(temp)) |>
  mutate(station_name = factor(station_name)) |>
  ggplot(aes(x = date, y = temp, colour = station_name)) +
  geom_point(size = 5) +
  # guides(colour = "none") +
  cecl_theme(nejm_pal = FALSE) +
  NULL
# hotter days becoming more frequent! But no days crazy above average or anything

# plot monthly rain data to identify potential strange ones
data_month |>
  filter(
    rain_mean > 15,
    month >= "1980-01-01",
    month < "2011-01-01"
  ) |>
  arrange(desc(rain_mean)) |>
  mutate(station_name = factor(station_name)) |>
  ggplot(aes(x = month, y = rain_mean, colour = station_name)) +
  geom_point(size = 5) +
  # guides(colour = "none") +
  cecl_theme(nejm_pal = FALSE) +
  NULL

data_month |>
  filter(
    rain_sum > 500,
    month >= "1980-01-01",
    month < "2011-01-01"
  ) |>
  arrange(desc(rain_sum)) |>
  mutate(station_name = factor(station_name)) |>
  ggplot(aes(x = month, y = rain_sum, colour = station_name)) +
  geom_point(size = 5) +
  # guides(colour = "none") +
  cecl_theme(nejm_pal = FALSE) +
  NULL



#### Load API Data ####

# Find station codes, e.g. Valencia
# TODO Match these stations with closest ones in daily data
stations <- aemet_stations()
stations |> filter(grepl("VALENCIA", nombre, ignore.case = TRUE))

# Normal climatology values; standard normal is 1981–2010
val_norm <- aemet_normal_clim("8414A") # replace with the Valencia station code

#### Fit marginal models to Spain data ####

# Idea: rolling empirical transformation:
# Take each month, fit linear model and remove trend.
#

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
devtools::load_all("../CeCl")

#### Metadata ####

decades <- seq(1960, 2010, by = 10)

#### Load Data ####

data <- readr::read_csv(
  "data/02_app/ecad_clean.csv.gz"
) |>
  # look at full years only
  filter(
    date > as_date("1960-01-01"),
    date < as_date("2024-01-01")
  ) |>
  # add useful columns for modelling
  mutate(
    year = as.numeric(substr(date, 1, 4)),
    month = lubridate::month(date),
    season_month = case_when(
      season == "Summer" ~ match(month, c(4, 5, 6, 7, 8, 9)),
      season == "Winter" ~ match(month, c(10, 11, 12, 1, 2, 3))
    )
  ) |>
  relocate(c(year, month), .after = date)

areas <- read_sf("data/02_app/spain_shapefile.geojson") |>
  filter(
    !ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla")
  )

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")


#### Initial Calculations ####

data <- data |>
  select(
    name = station_name, lon, lat,
    date, year, season_month, season,
    temp, rain, wind_speed, contains("drought")
  ) |>
  filter(
    if_all(
      c(temp, rain, drought_local, drought_global),
      \(x) !is.na(x)
    )
  ) |>
  mutate(month = lubridate::month(date))

# pull station names for looping over later
station_names <- unique(data$name)

#### Fit Gaussian distribution ####

#### Optional: Also fit GPD to exceedances ####
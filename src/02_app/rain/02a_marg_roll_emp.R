#### Fit marginal models to Spain data ####

# Idea: rolling empirical transformation:
# Take each month, fit linear model and remove trend.

# TODO Remove weeks with no rain and try again!

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
    date >= as_date("1960-01-01"),
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
  relocate(c(year, month), .after = date) |>
  filter(rain != 0) # ignore weeks with no rain??

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

#### Detrend ####

detrend_monthly <- \(data, var) {
  data |>
    group_by(name, month) |>
    group_modify(~ {
      fit <- lm(
        reformulate("year", response = var),
        data = .x
      )

      beta <- coef(fit)[["year"]]
      year_bar <- mean(.x$year, na.rm = TRUE)

      .x |>
        mutate(
          # fitted = predict(fit, .x),
          !!var := .data[[var]] -
            beta * (.data$year - year_bar)
        )
    }) |>
    ungroup()
}

temp_dt <- detrend_monthly(data, "temp")
rain_dt <- detrend_monthly(data, "rain")


#### Transform by season ####

# TODO Functionalise/join
temp_laplace_season <- temp_dt |>
  # group_by(name, month) |>
  group_by(name, season) |>
  mutate(
    F_hat = rank(temp) / (n() + 1),
    # temp_laplace = dlaplace(F_hat)[, 1]
    temp_laplace = qlaplace(F_hat)
  ) |>
  ungroup()

rain_laplace_season <- rain_dt |>
  # group_by(name, month) |>
  group_by(name, season) |>
  mutate(
    F_hat = rank(rain) / (n() + 1),
    # rain_laplace = dlaplace(F_hat)[, 1]
    rain_laplace = qlaplace(F_hat)
  ) |>
  ungroup()

# laplace_season <- cbind(
#   select(temp_laplace_season, name, date, year, month, season, temp_laplace),
#   select(rain_laplace_season, rain_laplace)
# )
laplace_season <- temp_laplace_season |>
  select(name, date, year, month, season, temp_laplace) |>
  inner_join(
    rain_laplace_season |>
      select(name, date, rain_laplace),
    by = c("name", "date")
  )

# laplace_monthly |>
#   filter(name == "Santander", season == "Summer") |>
#   # filter(name == "Santander", month == 12) |>
#   ggplot(aes(x = rain_laplace, y = temp_laplace)) +
#   geom_point()

# TODO Functionalise
# TODO Add other transformations to plot to compare
pdf("plots/02_app/roll_emp/01_laplace_trans_winter.pdf", width = 12, height = 8)
laplace_season |>
  filter(season == "Winter") |>
  mutate(name = as.factor(name)) |>
  group_split(name) |>
  lapply(\(x) {
    x |>
      ggplot(
        aes(x = temp_laplace, y = rain_laplace)
      ) +
      geom_point(alpha = 0.8) +
      cecl_theme() +
      ggtitle(paste0(x$name[[1]], " - Winter"))
  })
dev.off()

pdf("plots/02_app/roll_emp/01_laplace_trans_summer.pdf", width = 12, height = 8)
laplace_season |>
  filter(season == "Summer") |>
  mutate(name = as.factor(name)) |>
  group_split(name) |>
  lapply(\(x) {
    x |>
      ggplot(
        aes(x = temp_laplace, y = rain_laplace)
      ) +
      geom_point(alpha = 0.8) +
      cecl_theme() +
      ggtitle(paste0(x$name[[1]], " - Summer"))
  })
dev.off()

# convert to `cecl_marg` object so we can fit CE model easily using CeCl
marg_season_roll_emp <- laplace_season |>
  select(name, season, contains("laplace")) |>
  mutate(name = as.factor(name), season = as.factor(season)) |>
  group_split(season) |>
  lapply(\(x) {
    ret <- x |>
      group_split(name) |>
      lapply(\(y) {
        ret <- y |>
          select(contains("laplace")) |>
          as.matrix()
        colnames(ret) <- c("temp", "rain")
        ret
      })
    names(ret) <- unique(laplace_season$name)
    ret <- as_cecl_marg(ret)
    ret$season <- as.character(x$season[[1]])
    ret
  })
names(marg_season_roll_emp) <- unique(laplace_season$season)

saveRDS(
  marg_season_roll_emp,
  file = "data/02_app/marg_season_roll_emp.rds"
)

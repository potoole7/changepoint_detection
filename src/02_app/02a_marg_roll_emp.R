#### Fit marginal models to Spain data ####

# Idea: rolling empirical transformation:
# Take each month, fit linear model and remove trend.

# TODO Compare resulting Laplace transformation to "naive" empirical transform

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
# seasons <- c("Winter", "Summer")
seasons <- c("Winter", "Spring", "Summer", "Autumn")
# var_dep <- "rain"
# var_dep <- "drought_local"
var_dep <- "drought_local_rev"

# temp_var <- "temp"
temp_var <- "temp_max"


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
    ),
    drought_local_rev = -drought_local
  ) |>
  relocate(c(year, month), .after = date) |>
  # filter(rain != 0) # ignore weeks with no rain??
  # filter(!is.na(var_dep))
  filter(!is.na(.data[[var_dep]]))

# if specified, use maximum temperature rather than 90th quantile
if (temp_var == "temp_max") {
  data <- data |>
    mutate(
      temp_max = ifelse(is.infinite(temp_max), NA, temp_max),
      temp     = temp_max
    ) |>
    filter(!is.na(temp))
}

# reverse drought_local variable to give positive alpha values, if desired
if (var_dep == "drought_local_rev") {
  data <- data |>
    mutate(drought_local = -drought_local)
  var_dep <- "drought_local"
}

# check percentage of dates available across all locations for each season
station_count <- n_distinct(data$station_name)

date_coverage <- bind_rows(lapply(seasons, \(s) {
  season_data <- data |>
    filter(.data$season == s)

  all_season_dates <- season_data |>
    distinct(.data$date) |>
    pull(.data$date)

  valid_dates <- season_data |>
    # Add filters for missing measurements here if required
    distinct(.data$date, .data$station_name) |>
    count(.data$date, name = "n_stations") |>
    filter(.data$n_stations == station_count) |>
    pull(.data$date)

  tibble(
    season        = s,
    n_valid_dates = length(valid_dates),
    n_dates       = length(all_season_dates),
    perc          = n_valid_dates / n_dates
  )
}))

date_coverage # 94% - 80%: Fine!

# only keep dates available at every station (by season)
valid_dates <- data %>%
  distinct(season, date, station_name) %>%
  group_by(season) %>%
  mutate(n_stations_in_season = n_distinct(station_name)) %>%
  group_by(season, date) %>%
  filter(n_distinct(station_name) == first(n_stations_in_season)) %>%
  distinct(season, date)

data <- data %>%
  semi_join(valid_dates, by = c("season", "date"))

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
      # c(temp, rain, drought_local, drought_global),
      c(temp, !!var_dep),
      \(x) !is.na(x)
    )
  ) |>
  mutate(month = lubridate::month(date))

# pull station names for looping over later
station_names <- unique(data$name)

#### Detrend ####

# TODO Look into to see if this removes much!
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
          !!var := .data[[var]] -
            beta * (.data$year - year_bar)
        )
    }) |>
    ungroup()
}

temp_dt <- detrend_monthly(data, "temp")
# rain_dt <- detrend_monthly(data, "rain")
dep_dt <- detrend_monthly(data, var_dep)

# plot detrended data for a single station to check
wrap_plots(lapply(list(data, temp_dt), \(x) {
  x |>
    filter(name == "Valencia") |>
    ggplot(aes(x = date, y = temp)) +
    geom_line() +
    facet_wrap(~season) +
    cecl_theme()
}))

wrap_plots(lapply(list(data, dep_dt), \(x) {
  x |>
    filter(name == "Valencia") |>
    ggplot(aes(x = date, y = get(var_dep))) +
    geom_line() +
    facet_wrap(~season) +
    labs(y = var_dep) +
    cecl_theme()
}))

# data looks very similar for both original and "de-trended" data??

x <- data |>
  filter(name == data$name[1]) |>
  pull(!!var_dep)
y <- temp_dt |>
  filter(name == data$name[1]) |>
  pull(!!var_dep)

plot(x - y) # okay, apparently not the same!!

# check fitted slopes
wrap_plots(lapply(c("temp", var_dep), \(x) {
  data |>
    group_by(name, month) |>
    summarise(
      beta = coef(lm(.data[[x]] ~ year))[2],
      .groups = "drop"
    ) |>
    ggplot(aes(beta)) +
    geom_histogram() +
    cecl_theme() +
    ggtitle(x)
}))

# also overlay trend over a given year
ggplot(
  data |> filter(name == "Valencia"),
  aes(year, temp)
) +
  geom_point(alpha = 0.2) +
  facet_wrap(~month) +
  geom_smooth(method = "lm")
# clear positive trend in every month

# and look at the same for "detrended" data
temp_dt |>
  group_by(name, month) |>
  summarise(
    beta = coef(lm(temp ~ year))[2],
    .groups = "drop"
  )
# looks to have very flat lines!

# verify slopes are tiny
temp_dt |>
  group_by(name, month) |>
  summarise(
    beta = coef(lm(temp ~ year))[2],
    .groups = "drop"
  ) |>
  summarise(
    min = min(beta),
    max = max(beta),
    mean = mean(beta),
    sd = sd(beta)
  )


#### Transform by month ####

# Previously season, now month!

# temp_laplace_season <- temp_dt |>
#   # group_by(name, month) |>
#   group_by(name, season) |>
#   mutate(
#     F_hat = rank(temp) / (n() + 1),
#     laplace = qlaplace(F_hat)
#   ) |>
#   ungroup()
# dep_laplace_season <- dep_dt |>
#   group_by(name, season) |>
#   mutate(
#     F_hat = rank(
#       .data[[var_dep]],
#       ties.method = "average",
#       na.last = "keep"
#     ) / (sum(!is.na(.data[[var_dep]])) + 1),
#     laplace = qlaplace(F_hat)
#   ) |>
#   ungroup()


# function to transform by season to Laplace margins using rank transform
trans_fun <- \(df, var) {
  dep_laplace_season <- df |>
    # group_by(name, season) |>
    group_by(name, month) |>
    mutate(
      n_nonmissing = sum(!is.na(.data[[var]])),
      F_hat = rank(
        .data[[var]],
        ties.method = "average",
        na.last = "keep"
      ) / (n_nonmissing + 1),
      laplace = qlaplace(F_hat)
    ) |>
    ungroup() |>
    select(-n_nonmissing)
}

# temp_laplace_season <- trans_fun(df, "temp")
# dep_laplace_season <- trans_fun(df, var_dep)
temp_laplace_month <- trans_fun(temp_dt, "temp")
dep_laplace_month <- trans_fun(dep_dt, var_dep)


# join data
# laplace_season <- temp_laplace_season |>
#   select(name, date, year, month, season, temp = laplace) |>
#   inner_join(
#     dep_laplace_season |>
#       select(name, date, dep = laplace),
#     # select(name, date, year, month, season, dep = laplace),
#     by = c("name", "date")
#   ) |>
#   # rename dep to var_dep
#   rename(!!var_dep := dep) |>
#   identity()
laplace_month <- temp_laplace_month |>
  select(name, date, year, month, season, temp = laplace) |>
  inner_join(
    dep_laplace_month |>
      select(name, date, dep = laplace),
    by = c("name", "date"),
    relationship = "one-to-one"
  ) |>
  rename("{var_dep}" := dep)

laplace_season <- laplace_month |>
  filter(season %in% seasons)

# Check that each month at each station spans approximately a standard Laplace distribution (medians should be ~0)
laplace_month |>
  group_by(name, month) |>
  summarise(
    temp_median = median(temp, na.rm = TRUE),
    dep_median = median(.data[[var_dep]], na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
# successful!

# Function to plot bivariate Laplace data
plot_trans_lst <- \(df, season) {
  p_lst <- df |>
    filter(season == !!season) |>
    mutate(name = as.factor(name)) |>
    group_split(name) |>
    lapply(\(x) {
      x |>
        ggplot(
          # aes(x = temp_laplace, y = rain_laplace)
          aes(x = temp, y = get(var_dep))
        ) +
        geom_point(alpha = 0.8) +
        cecl_theme() +
        labs(x = "temp", y = var_dep) +
        ggtitle(paste0(x$name[[1]], " - ", season))
    })
  return(p_lst)
}


# TODO Add other transformations to plot to compare
pdf("plots/02_app/roll_emp/01_laplace_trans_winter.pdf", width = 12, height = 8)
plot_trans_lst(laplace_season, "Winter")
dev.off()

pdf("plots/02_app/roll_emp/01_laplace_trans_spring.pdf", width = 12, height = 8)
plot_trans_lst(laplace_season, "Spring")
dev.off()

pdf("plots/02_app/roll_emp/01_laplace_trans_summer.pdf", width = 12, height = 8)
plot_trans_lst(laplace_season, "Summer")
dev.off()

pdf("plots/02_app/roll_emp/01_laplace_trans_autumn.pdf", width = 12, height = 8)
plot_trans_lst(laplace_season, "Autumn")
dev.off()

# convert to `cecl_marg` object so we can fit CE model easily using CeCl
# marg_season_roll_emp <- laplace_season |>
#   # select(name, season, contains("laplace")) |>
#   select(name, season, temp, all_of(var_dep)) |>
#   mutate(name = as.factor(name), season = factor(season, levels = seasons)) |>
#   group_split(season) |>
#   lapply(\(x) {
#     ret_names <- unique(x$name)
#     ret <- x |>
#       group_split(name) |>
#       # lapply(\(y) {
#       mclapply(\(y) {
#         # browser()
#         # print(y$season[1])
#         # print(y$name[1])
#         ret1 <- y |>
#           # select(contains("laplace")) |>
#           select(temp, all_of(var_dep)) |>
#           as.matrix()
#         # colnames(ret) <- c("temp", "rain")
#         colnames(ret1) <- c("temp", var_dep)
#         ret1
#       })
#     names(ret) <- ret_names
#     ret <- as_cecl_marg(ret)
#     ret$season <- as.character(x$season[[1]])
#     ret
#   })
# names(marg_season_roll_emp) <- unique(laplace_season$season)
# safer version
marg_season_roll_emp <- laplace_season |>
  select(name, date, season, temp, all_of(var_dep)) |>
  mutate(
    name = factor(name),
    season = factor(season, levels = seasons)
  ) |>
  group_split(season) |>
  lapply(\(x) {
    station_groups <- x |>
      group_by(name)

    ret_names <- station_groups |>
      group_keys() |>
      pull(name) |>
      as.character()

    station_data <- station_groups |>
      group_split()
    
    original_dates <- lapply(
      station_data,
      \(y) {
        y |>
          arrange(date) |>
          pull(date)
      }
    )

    ret <- mclapply(
      station_data,
      \(y) {
        ret1 <- y |>
          arrange(date) |>
          select(temp, all_of(var_dep)) |>
          as.matrix()

        colnames(ret1) <- c("temp", var_dep)
        ret1
      }
    )

    names(ret) <- ret_names
    names(original_dates) <- ret_names

    ret <- as_cecl_marg(ret)
    # keep dates and season to ensure we can match later
    ret$dates <- original_dates
    ret$season <- as.character(x$season[[1]])

    ret
  })
season_names <- laplace_season |>
  mutate(season = factor(season, levels = seasons)) |>
  distinct(season) |>
  arrange(season) |>
  pull(season) |>
  as.character()

names(marg_season_roll_emp) <- season_names

saveRDS(
  marg_season_roll_emp,
  file = "data/02_app/marg_season_roll_emp.rds"
)

# if (var_dep == "drought_local") {
#   saveRDS(
#     marg_season_roll_emp,
#     file = "data/02_app/marg_season_roll_emp.rds"
#   )
# } else if (var_dep == "drought_local_rev") {
#   saveRDS(
#     marg_season_roll_emp,
#     file = "data/02_app/marg_season_roll_emp_rev.rds"
#   )
# }

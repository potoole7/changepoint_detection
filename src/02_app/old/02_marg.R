#### Fit marginal models to Spain data ####

# Perform empirical transformation to Laplace distribution (done)

# Marginal TODOs:
# Add (season) month as well as year as predictor (done)
# TODO Use quantile regression with qgam for choosing dependence threshold (see notes)
#      over time (year + season month)
# TODO Compare models with varying vs fixed shape parameters (and fixed vs varying thresh)
# TODO Look into outliers for Daroca and Zamora (Winter) for evgam transformation!
# TODO Implement QQ plots
# TODO Mess around with threshold if they're not great ..


# Fit GPD with evgam to temperature and rain (done)
# TODO Select marginal thresholds
# TODO - Test for q = 0.9
# TODO - Look at mean residual life plots and stability plots for stationary GPD
# TODO Add year as predictor
# TODO Plot:
# TODO - Parameter values (map for given year, also maybe time series by year??)
# TODO - Return levels on map
# TODO - Uncertainty (somehow??)
# TODO - Diagnostic plots (QQ)

# TODO Perform empirical transformation as per Jordan's idea

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
library(qgam)
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
  )

#### Empirical transformation ####

# Split marginal models by season
marg_season <- data |>
  bind_rows(
    data |>
      mutate(season = "Year")
  ) |>
  mutate(season = factor(season, levels = c("Winter", "Summer", "Year"))) |>
  group_split(season) |>
  lapply(\(x) {
    # browser()
    print(x$name[1])
    print(x$season[1])
    CeCl::cecl_marg(
      # select(x, -rain),
      select(x, name, temp, rain),
      mult_col = "name",
      # vars = c("temp", "rain", "wind_speed", "drought_local", "drought_global"),
      # vars = c("temp", "rain", "drought_local", "drought_global"),
      vars = c("temp", "rain"),
      thresh_method = "none",
      marg_method = "ecdf"
    )
  })
names(marg_season) <- c("Winter", "Summer", "Year")

saveRDS(
  marg_season,
  file = "data/02_app/marg_season_emp.rds"
)

# also do for each decade
marg_season_decade <- data |>
  bind_rows(
    data |>
      mutate(season = "Year")
  ) |>
  mutate(
    season = factor(season, levels = c("Winter", "Summer", "Year")),
    decade = factor(floor(year(date) / 10) * 10, levels = decades)
  ) |>
  filter(!is.na(decade)) |>
  group_split(decade, season) |>
  lapply(\(x) {
    ret <- cecl_marg(
      x,
      mult_col = "name",
      vars = c("temp", "rain", "drought_local", "drought_global"),
      thresh_method = "none",
      marg_method = "ecdf"
    )
    # label object with season and decade to ensure it's identifiable
    ret$season <- x$season[[1]]
    ret$decade <- x$decade[[1]]
    ret
  })
# names(marg_season_decade) <- c("Winter", "Summer", "Year")
names(marg_season_decade) <- crossing(
  "decade" = decades,
  "season" = c("Winter", "Summer", "Year")
) |>
  mutate(label = paste(decade, season, sep = "_")) |>
  pull(label)

saveRDS(
  marg_season_decade,
  file = "data/02_app/marg_season_decade_emp.rds"
)


#### Marginal evgam model: Choose threshold ####

# TODO Do more involved analysis, but for now just use reasonable quantile
q <- 0.9 # for 60 years of weekly data, gives 156 exceedances

# function to threshold data
# TODO Allow for location-specific quantiles in CeCl
thresh_data <- \(data, var, season, q = 0.9) {
  data_spec <- data |>
    filter(season == !!season)

  ret <- data_spec |>
    group_by(name) |>
    mutate(thresh = quantile(.data[[var]], q)) |>
    ungroup() |>
    mutate(excess = .data[[var]] - thresh) |>
    filter(excess > 0) |>
    select(name, lon, lat, date, year, season_month, season, thresh, excess)

  # names(ret)[names(ret) == "excess"] <- var
  return(ret)
}

# test
thresh_data(data, "rain", "Winter", q = 0.9)

season_var_df <- tidyr::crossing(
  "var" = c("rain", "temp"),
  "season" = c("Winter", "Summer")
) |>
  mutate(lst_name = paste0(var, "_", season))

thresh_data_season_var <- season_var_df |>
  mutate(var = as.factor(var), season = as.factor(season)) |>
  group_split(var, season) |>
  lapply(\(x) {
    with(x, thresh_data(data, as.character(var), as.character(season)))
  })

names(thresh_data_season_var) <- season_var_df$lst_name


#### evgam: Fit ####

# TODO: Fix!!
# debugonce(cecl_marg)
# marg_evgam <- cecl_marg(
#   x_spec,
#   thresh_method = "none", # already thresholded
#   marg_method = "evgam",
#   # model formulae for GPD scale and shape parameters
#   marg_args = list(f = list("excess ~ year", "~ year")),
#   ncores = 1,
#   ret_obj = TRUE
# )

# x <- thresh_data_season_var[[1]]
# function to fit GPD with evgam and make predictions for each year
evgam_fit <- \(
  x,
  # TODO Test different models
  f = list(excess ~ s(year) + s(season_month, k = 5), ~1),
  name,
  ret_mod = FALSE
) {
  # pull useful names
  x_spec <- x |>
    select(name, year, season_month, thresh, excess) |>
    filter(name == !!name)

  # formula
  # f <- list(excess ~ year, ~year)

  # fit evgam model
  m <- evgam::evgam(f, data = x_spec, family = "gpd")

  # create predictions for unique rows in pred_data (ensures one pred per loc)
  predictors <- m$predictor.names
  pred_dat_distinct <- x_spec |>
    distinct(across(all_of(predictors)))
  predictions <- cbind(
    pred_dat_distinct,
    predict(m, pred_dat_distinct, type = "response")
  )

  # TODO Also return uncertainty!!
  predictions$thresh <- x_spec$thresh[[1]]
  predictions <- predictions |>
    mutate(name = x_spec$name[1], thresh = x_spec$thresh[[1]]) |>
    relocate(name)
  rownames(predictions) <- NULL
  if (ret_mod) {
    return(list(
      "evgam_fit" = m,
      "marginal"  = predictions
    ))
  } else {
    return(predictions)
  }
}

name <- "Valencia"
# name = "Avila"
# name = "Lleida"

# # # test for Valencia (and other locations)
# obj1 <- evgam_fit(
#   thresh_data_season_var$temp_Winter,
#   f = list(excess ~ s(year) + s(season_month, k = 3), ~1),
#   name = name,
#   ret_mod = TRUE
# )
# fit1 <- obj1[[1]]
# summary(fit1)
# plot(fit1, pages = 1)
#
# obj2 <- evgam_fit(
#   thresh_data_season_var$temp_Winter,
#   f = list(excess ~ s(year) + s(season_month, k = 3), ~ s(year)),
#   name = name,
#   ret_mod = TRUE
# )
# fit2 <- obj2[[1]]
# summary(fit2)
# plot(fit2, pages = 1)
#
# AIC(fit1, fit2) #
#
# # plot
# bind_rows(
#   mutate(obj1$marginal, type = "fixed_shape"),
#   mutate(obj2$marginal, type = "varying_shape")
# ) |>
#   pivot_longer(c(scale, shape), names_to = "parameter") |>
#   ggplot(aes(x = year, y = value, colour = type)) +
#   geom_line() +
#   facet_wrap(~parameter, scale = "free") +
#   cecl_theme()


# fit location-specific GPDs for each variable-season setup
station_names <- unique(data$name)
gpd_season_var <- lapply(thresh_data_season_var, \(x) {
  bind_rows(mclapply(station_names, \(y) {
    dat <- filter(x, name == y)
    n_month <- length(unique(dat$season_month))
    k_month <- min(5, n_month - 1) # must be greater than number of months avail

    f <- list(
      as.formula(paste0(
        "excess ~ s(year) + s(season_month, k = ", k_month, ")"
      )),
      ~1
    )

    evgam_fit(x, f = f, name = y)
  }))
})
names(gpd_season_var) <- names(thresh_data_season_var)

# TODO Plot time series
# TODO Plot maps for some years


#### Laplace transform ####

# function to perform marginal transformation
trans_marg_ns_one <- \(data_orig,
  gpd_fit,
  var,
  season,
  station) {
  # take data for given station and season
  dat <- data_orig |>
    filter(name == station, season == !!season) |>
    arrange(date) |>
    select(name, lon, lat, date, season, all_of(var)) |>
    filter(!is.na(.data[[var]]))

  # take GPD fit to given station
  gpd <- gpd_fit |>
    filter(name == station)

  # Perform piecewise transformation, calculating estimated dist function
  F_hat <- p_gpd_ecdf_ns(
    dat = dat,
    gpd = gpd,
    var = var
  )

  # add to data
  dat |>
    mutate(
      F_hat = F_hat,
      laplace = dlaplace(F_hat)[, 1] # calculate laplace data
    )
}

# Semi-parametric CDF with year-varying GPD parameters
p_gpd_ecdf_ns <- \(dat, gpd, var) {
  x <- dat[[var]]
  u <- gpd$thresh[1]

  # Join year-varying parameters
  dat_par <- dat |>
    mutate(year = as.numeric(substr(date, 1, 4))) |>
    left_join(
      gpd |>
        select(year, scale, shape),
      by = "year"
    )

  # Empirical CDF for all observations
  ord <- order(x)
  x_sort <- x[ord]
  m <- length(x)

  ecdf_vals <- seq_len(m) / (m + 1)
  ecdf_orig <- numeric(m)
  ecdf_orig[ord] <- ecdf_vals

  # Initialise
  F_hat <- ecdf_orig

  # Parametric GPD tail above threshold
  above <- x > u

  if (any(above, na.rm = TRUE)) {
    excess <- x[above] - u
    sigma <- dat_par$scale[above]
    xi <- dat_par$shape[above]

    surv_gpd <- ifelse(
      abs(xi) < 1e-8,
      exp(-excess / sigma),
      pmax(1 + xi * excess / sigma, 0)^(-1 / xi)
    )

    p_exc <- mean(x > u, na.rm = TRUE)

    F_hat[above] <- 1 - p_exc * surv_gpd
  }

  F_hat
}

# perform Laplace transformation for each variable and season
laplace_season_var <- lapply(seq_along(gpd_season_var), \(i) {
  # bind_rows(mclapply(station_names, \(x) { # loop through stations
  bind_rows(lapply(station_names, \(x) { # loop through stations
    trans_marg_ns_one(
      data_orig = data,
      gpd_fit = gpd_season_var[[i]],
      var = season_var_df$var[[i]],
      season = season_var_df$season[[i]],
      station = x
    )
  })) |>
    select(name, laplace)
})
names(laplace_season_var) <- names(gpd_season_var)

# function to join together rain and temperature transformed data for a season
trans_fun <- \(x, season) {
  # join for each variable
  laplace_df <- cbind(
    # laplace_season_var$rain_ |>
    x[[paste0("temp_", season)]] |>
      rename(temp = laplace),
    # laplace_season_var$temp_ |>
    x[[paste0("rain_", season)]] |>
      select(-name) |>
      rename(rain = laplace)
  ) |>
    mutate(name = factor(name))


  laplace_lst <- group_split(laplace_df, name, .keep = FALSE) |>
    lapply(\(x) {
      ret <- as.matrix(x)
      colnames(ret) <- names(x)
      ret
    })
  names(laplace_lst) <- unique(laplace_df$name)

  # convert to cecl_marg format so it works with CE
  return(as_cecl_marg(laplace_lst))
}

marg_season_evgam <- list(
  "Winter" = trans_fun(laplace_season_var, "Winter"),
  "Summer" = trans_fun(laplace_season_var, "Summer")
)

saveRDS(
  marg_season_evgam,
  file = "data/02_app/marg_season_evgam.rds"
)

# appears to be better signal/stronger alpha (note scale of temperature)
plot_diff <- \(season) {
  lapply(station_names, \(x) {
    df_evgam <- data.frame(marg_season_evgam[[season]]$transformed[[x]])
    df_ecdf <- data.frame(marg_season[[season]]$transformed[[x]])

    df_plot <- bind_rows(
      mutate(df_evgam, type = "evgam"),
      mutate(df_ecdf, type = "ecdf")
    )

    df_plot |>
      ggplot(aes(x = temp, y = rain, colour = type)) +
      geom_point(alpha = 0.8) +
      cecl_theme() +
      ggtitle(paste0(x, " - ", season))
  })
}


pdf("plots/02_app/03_laplace_trans_compare_winter.pdf", width = 12, height = 8)
plot_diff("Winter")
dev.off()

pdf("plots/02_app/03_laplace_trans_compare_summer.pdf", width = 12, height = 8)
plot_diff("Summer")
dev.off()


#### Plot strange countries ####

data |>
  filter(name %in% c("Daroca", "Zamora")) |>
  ggplot(aes(x = temp, y = rain)) +
  geom_point() +
  facet_wrap(name ~ season, scale = "free") +
  cecl_theme()

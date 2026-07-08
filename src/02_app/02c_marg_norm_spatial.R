#### Fit marginal models to Spain data ####

# TODO Could add a station specific random effect?
# s(name, bs = "re")

# Idea: Fit marginal models to the Spain data for temperature and rainfall,
# using Gaussian and Gamma distributions respectively.
# The models will account for seasonal variations and temporal trends.

# This time, also include a spatial component in the marginal models, using
# the station coordinates (lon, lat) as covariates.


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
    date < as_date("2024-01-01"),
    rain > 0 # only keep days with rain for rainfall modelling
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


#### Fit Gaussian and Gamma distributions ####

# function to find the minimum valid k pair for year and season_month
find_min_k_pair <- \(
  # k_spatial_vals = c(10, 20, 30, 40),
  k_spatial_vals = c(40),
  k_year_vals = 3:9,
  k_season_month_vals = 3:5,
  response_name,
  x_spec,
  fam,
  ...
) {
  # for (ky in k_year_vals) {
  #   for (km in k_season_month_vals) {
  #     f_try <- list(
  #       reformulate(
  #         termlabels = c(
  #           sprintf("s(year, k = %s)", ky),
  #           sprintf("s(season_month, k = %s)", km)
  #         ),
  #         response = response_name
  #       ),
  #       ~1
  #     )
  #
  #     fit <- try(
  #       mgcv::gam(f_try, data = x_spec, family = fam, method = "REML", ...),
  #       silent = TRUE
  #     )
  #
  #     if (!inherits(fit, "try-error")) {
  #       return(list(
  #         k_year = ky,
  #         k_season_month = km,
  #         formula = f_try,
  #         fit = fit
  #       ))
  #     }
  #   }
  # }

  for (ks in k_spatial_vals) {
    for (ky in k_year_vals) {
      for (km in k_season_month_vals) {
        f_try <- list(
          reformulate(
            termlabels = c(
              sprintf("s(lon, lat, k = %s)", ks),
              sprintf("s(year, k = %s)", ky),
              sprintf("s(season_month, k = %s)", km)
            ),
            response = response_name
          ),
          ~1
        )

        fit <- try(
          mgcv::gam(f_try, data = x_spec, family = fam, method = "REML", ...),
          silent = TRUE
        )

        if (!inherits(fit, "try-error")) {
          return(list(
            k_spatial = ks,
            k_year = ky,
            k_season_month = km,
            formula = f_try,
            fit = fit
          ))
        }
      }
    }
  }
  stop("No valid k pair found.")
}

# Function to compute quantile-quantile data for mgcv models
# qq_data_mgcv <- \(m, x_spec, response_name, family) {
#   y <- x_spec[[response_name]]
#   p <- predict(m, x_spec, type = "response")
#
#   if (family == "gaulss") {
#     mu <- p[, 1]
#     sigma <- p[, 2]
#
#     qres <- qnorm(pnorm(y, mean = mu, sd = sigma))
#   }
#
#   if (family == "gammals") {
#     mu <- p[, 1]
#     scale <- exp(p[, 2]) # gammals second parameter is log scale
#
#     shape <- 1 / scale
#     rate <- shape / mu
#
#     qres <- qnorm(pgamma(y, shape = shape, rate = rate))
#   }
#
#   tibble(
#     observed = sort(qres),
#     theoretical = qnorm(ppoints(length(qres)))
#   )
# }

# qq_data_mgcv <- function(m, x_spec, response_name, family) {
#   y <- x_spec[[response_name]]
#   p <- predict(m, x_spec, type = "response")
#
#   if (family == "gaulss") {
#     mu <- p[, 1]
#     sigma <- 1 / p[, 2]
#
#     qres <- (y - mu) / sigma
#   }
#
#   if (family == "gammals") {
#     mu <- p[, 1]
#     scale <- exp(p[, 2])
#
#     shape <- 1 / scale
#     rate <- shape / mu
#
#     qres <- qnorm(pgamma(y, shape = shape, rate = rate))
#   }
#
#   tibble(
#     observed = sort(qres),
#     theoretical = qnorm(ppoints(length(qres)))
#   )
# }

qq_data_mgcv <- function(m, x_spec, response_name, family) {
  p <- predict(m, x_spec, type = "response")
  y <- x_spec[[response_name]]

  if (family == "gaulss") {
    mu <- p[, 1]
    sigma <- 1 / p[, 2]
    qres <- (y - mu) / sigma
  }

  if (family == "gammals") {
    mu <- p[, 1]
    scale <- exp(p[, 2])

    shape <- 1 / scale
    rate <- shape / mu

    qres <- qnorm(pgamma(y, shape = shape, rate = rate))
  }

  x_spec |>
    mutate(qres = qres) |>
    group_by(name) |>
    arrange(qres, .by_group = TRUE) |>
    mutate(
      i = row_number(),
      n = n(),
      observed = qres,
      theoretical = qnorm(ppoints(n))
    ) |>
    ungroup() |>
    select(name, theoretical, observed)
}

# # for testing
# x <- data |> filter(season == "Winter")
# # # TODO Vary sigma as well as mean?
# f <- list(temp ~ s(year, k = 3) + s(season_month, k = 3), ~1)
# # name <- "Valencia"
# name <- "Badajoz Aeropuerto"
# # ret_mod <- TRUE
# family <- "gaulss"
# response_name <- "temp"
# k_year_vals <- 3:9
# k_season_month_vals <- 3:5

# Function to fit either gaulss or gammals model for a given station name
mgcv_fit <- \(
  x,
  # name,
  response_name = c("temp", "rain"),
  ret_mod = FALSE,
  family = c("gaulss", "gammals"),
  # k_spatial_vals = 3:9,
  # k_spatial_vals = c(10, 20, 30, 40),
  k_spatial_vals = c(40),
  k_year_vals = 3:9,
  k_season_month_vals = 3:5,
  ...
) {
  family <- match.arg(family)
  response_name <- match.arg(response_name)

  stopifnot(
    (response_name == "temp" && family == "gaulss") ||
      (response_name == "rain" && family == "gammals")
  )

  x_spec <- x |>
    select(name, lon, lat, year, season_month, all_of(response_name)) |>
    # filter(name == !!name) |>
    identity()

  if (nrow(x_spec) == 0) {
    return(NULL)
  }

  if (family == "gammals") {
    x_spec <- x_spec |>
      filter(.data[[response_name]] > 0)

    if (nrow(x_spec) == 0) {
      return(NULL)
    }
  }

  fam <- switch(family,
    gaulss = mgcv::gaulss(),
    gammals = mgcv::gammals()
  )

  min_fit <- find_min_k_pair(
    k_spatial_vals = k_spatial_vals,
    k_year_vals = k_year_vals,
    k_season_month_vals = k_season_month_vals,
    response_name = response_name,
    x_spec = x_spec,
    fam = fam,
    ...
  )

  m <- min_fit$fit
  f <- min_fit$formula

  predictors <- attr(terms(m), "term.labels")

  pred_dat_distinct <- x_spec |>
    # distinct(across(all_of(predictors)))
    distinct(name, across(all_of(predictors)))

  pred <- predict(m, pred_dat_distinct, type = "response")

  predictions <- bind_cols(
    pred_dat_distinct,
    as.data.frame(pred)
  )

  # names(predictions) <- c(predictors, "mu", "sigma")
  if (family == "gaulss") {
    names(predictions) <- c("name", predictors, "mu", "inv_sigma")
  }

  if (family == "gammals") {
    names(predictions) <- c("name", predictors, "mu", "log_scale")
  }

  predictions <- predictions |>
    mutate(
      # name = x_spec$name[1],
      k_spatial = min_fit$k_spatial,
      k_year = min_fit$k_year,
      k_season_month = min_fit$k_season_month
    ) |>
    relocate(name, k_spatial, k_year, k_season_month)

  rownames(predictions) <- NULL

  # compute quantile-quantile data for the fitted model
  # qq_dat <- qq_data_mgcv(
  #   m = m,
  #   x_spec = x_spec,
  #   response_name = response_name,
  #   family = family
  # )
  qq_dat <- qq_data_mgcv(
    m = m,
    x_spec = x_spec,
    response_name = response_name,
    family = family
  )

  if (ret_mod) {
    list(
      mgcv_fit = m,
      marginal = predictions,
      qq = qq_dat,
      formula = f,
      k_year = min_fit$k_year,
      k_season_month = min_fit$k_season_month
    )
  } else {
    predictions
  }
}

# # test
# debugonce(mgcv_fit)
# mgcv_fit(
#   x = data |> filter(season == "Winter"),
#   # name = "Badajoz Aeropuerto",
#   response_name = "temp",
#   family = "gaulss",
#   ret_mod = TRUE
# )

# fit Gaussian marginal models to temperature for each station
# debugonce(mgcv_fit)
gaulss_season_var_spat <- data |>
  mutate(season = factor(season, levels = c("Winter", "Summer"))) |>
  group_split(season) |>
  # mclapply(\(x) {
  lapply(\(x) {
    # ret <- lapply(station_names, \(y) {
    #   mgcv_fit(
    #     x,
    #     name = y,
    #     response_name = "temp",
    #     family = "gaulss",
    #     ret_mod = TRUE
    #   )
    # })
    # names(ret) <- station_names
    # ret
    # })
    mgcv_fit(
      x,
      response_name = "temp",
      family = "gaulss",
      ret_mod = TRUE
    )
  })
names(gaulss_season_var_spat) <- c("Winter", "Summer")

# fit Gamma marginal models to rain for each station
# TODO Note that for Gamma we may have log scale/sigma parameter?
gammals_season_var_spat <- data |>
  mutate(season = factor(season, levels = c("Winter", "Summer"))) |>
  group_split(season) |>
  lapply(\(x) {
    # ret <- lapply(station_names, \(y) {
    #   mgcv_fit(
    #     x,
    #     name = y,
    #     response_name = "rain",
    #     family = "gammals",
    #     ret_mod = TRUE
    #   )
    # })
    # names(ret) <- station_names
    # ret
    mgcv_fit(
      x,
      response_name = "rain",
      family = "gammals",
      ret_mod = TRUE
    )
  })
names(gammals_season_var_spat) <- c("Winter", "Summer")


#### Validate ####

# function to compute QQ envelope for a given sample size n,
# number of simulations nsim, and confidence level conf
qq_envelope <- function(n, nsim = 1000, conf = 0.95) {
  sims <- replicate(nsim, sort(rnorm(n)))

  tibble(
    theoretical = qnorm(ppoints(n)),
    lower = apply(sims, 1, quantile, probs = (1 - conf) / 2),
    upper = apply(sims, 1, quantile, probs = 1 - (1 - conf) / 2)
  )
}

# pull QQ data for each sesaon and station
# qq_data <- bind_rows(lapply(seq_along(gaulss_season_var_spat), \(i) {
#   season_name <- names(gaulss_season_var_spat)[i]
#   bind_rows(lapply(seq_along(station_names), \(j) {
#     bind_rows(
#       mutate(gaulss_season_var_spat[[i]][[j]]$qq, family = "gaulss"),
#       mutate(gammals_season_var_spat[[i]][[j]]$qq, family = "gammals")
#     ) |>
#       mutate(
#         season = season_name,
#         name = station_names[j]
#       )
#   }))
# })) |>
#   # calculate uncertainty envelopes for each station and season
#   group_by(name, family, season) |>
#   group_modify(\(.x, .y) {
#     bind_cols(
#       .x |> arrange(theoretical),
#       qq_envelope(nrow(.x)) |> select(lower, upper)
#     )
#   }) |>
#   ungroup()
qq_data_spat <- bind_rows(lapply(names(gaulss_season_var_spat), function(season_name) {
  bind_rows(
    gaulss_season_var_spat[[season_name]]$qq |>
      mutate(family = "gaulss", season = season_name),
    gammals_season_var_spat[[season_name]]$qq |>
      mutate(family = "gammals", season = season_name)
  )
})) |>
  group_by(name, family, season) |>
  group_modify(\(.x, .y) {
    bind_cols(
      .x |> arrange(theoretical),
      qq_envelope(nrow(.x)) |> select(lower, upper)
    )
  }) |>
  ungroup()

qq_data_spat |>
  filter(name == "Santander") |>
  mutate(
    ind = paste0(family, ", ", season),
    ind = factor(ind, levels = c(
      "gaulss, Winter", "gaulss, Summer",
      "gammals, Winter", "gammals, Summer"
    ))
  ) |>
  ggplot(aes(theoretical, observed)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_point(alpha = 0.4) +
  facet_wrap(~ind) +
  cecl_theme()


# run for all stations and save results
pdf("plots/02_app/mgcv_spat/marginal_qq_plots.pdf", width = 10, height = 8)
lapply(station_names, \(x) {
  qq_data_spat |>
    filter(name == x) |>
    mutate(
      ind = paste0(family, ", ", season),
      ind = factor(ind, levels = c(
        "gaulss, Winter", "gaulss, Summer",
        "gammals, Winter", "gammals, Summer"
      ))
    ) |>
    ggplot(aes(x = theoretical, y = observed)) +
    geom_point() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    facet_wrap(~ind, scales = "free") +
    cecl_theme() +
    labs(title = x)
})

dev.off()

# conclusions:
# Way worse than individual models for each location! Some particular 
# stinkers as well


#### Laplace transformation ####

p_mgcv_spat_marg <- \(dat, fit, var, family) {
  p <- predict(fit, newdata = dat, type = "response")
  x <- dat[[var]]

  if (family == "gaulss") {
    mu <- p[, 1]
    sigma <- 1 / p[, 2]

    F_hat <- pnorm(x, mean = mu, sd = sigma)
  }

  if (family == "gammals") {
    mu <- p[, 1]
    scale <- exp(p[, 2])

    shape <- 1 / scale
    rate <- shape / mu

    F_hat <- pgamma(x, shape = shape, rate = rate)
  }

  pmin(pmax(F_hat, 1e-6), 1 - 1e-6)
}

trans_mgcv_spat_one <- \(data_orig,
  mgcv_fit,
  var,
  family,
  station,
  season
) {
  dat <- data_orig |>
    select(name, lon, lat, date, year, season_month, season, all_of(var)) |>
    filter(name == !!station, season == !!season) |>
    arrange(date) |>
    filter(!is.na(.data[[var]]))

  if (family == "gammals") {
    dat <- dat |> filter(.data[[var]] > 0)
  }

  stopifnot(nrow(dat) > 0)

  # fit_obj <- mgcv_fit[[station]]
  # fit_obj <- mgcv_fit[[season]]
  fit_obj <- mgcv_fit

  F_hat <- p_mgcv_spat_marg(
    dat = dat,
    fit = fit_obj$mgcv_fit,
    var = var,
    family = family
  )

  dat |>
    mutate(
      F_hat = F_hat,
      laplace = dlaplace(F_hat)[, 1]
    )
}

# for (s in names(gaulss_season_var_spat)) {
#   names(gaulss_season_var_spat[[s]]) <- station_names
#   names(gammals_season_var_spat[[s]]) <- station_names
# }

laplace_mgcv_season_var_spat <- list(
  temp_Winter = bind_rows(mclapply(station_names, \(st) {
    trans_mgcv_spat_one(
      data_orig = data,
      mgcv_fit = gaulss_season_var_spat$Winter,
      var = "temp",
      family = "gaulss",
      season = "Winter",
      station = st
    )
  })) |> select(name, laplace),
  temp_Summer = bind_rows(mclapply(station_names, \(st) {
    trans_mgcv_spat_one(
      data_orig = data,
      mgcv_fit = gaulss_season_var_spat$Summer,
      var = "temp",
      family = "gaulss",
      season = "Summer",
      station = st
    )
  })) |> select(name, laplace),
  rain_Winter = bind_rows(mclapply(station_names, \(st) {
    trans_mgcv_spat_one(
      data_orig = data,
      mgcv_fit = gammals_season_var_spat$Winter,
      var = "rain",
      family = "gammals",
      season = "Winter",
      station = st
    )
  })) |> select(name, laplace),
  rain_Summer = bind_rows(mclapply(station_names, \(st) {
    trans_mgcv_spat_one(
      data_orig = data,
      mgcv_fit = gammals_season_var_spat$Summer,
      var = "rain",
      family = "gammals",
      season = "Summer",
      station = st
    )
  })) |> select(name, laplace)
)

trans_fun <- \(x, season) {
  laplace_df <- cbind(
    x[[paste0("temp_", season)]] |>
      rename(temp = laplace),
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

  as_cecl_marg(laplace_lst)
}

marg_season_mgcv_spat <- list(
  Winter = trans_fun(laplace_mgcv_season_var_spat, "Winter"),
  Summer = trans_fun(laplace_mgcv_season_var_spat, "Summer")
)

saveRDS(
  marg_season_mgcv_spat,
  file = "data/02_app/marg_season_mgcv_spat.rds"
)

# plot one transformed dataset
# data.frame(marg_season_mgcv$Winter$transformed$`Badajoz Aeropuerto`) |>
data.frame(marg_season_mgcv_spat$Winter$transformed$Santander) |>
  ggplot(aes(x = temp, y = rain)) +
  geom_point(alpha = 0.8) +
  cecl_theme() +
  ggtitle(paste0("Winter"))


#### GPD: threshold ####

q <- 0.9
# q <- 0.92
# q <- 0.95

# Idea: For exceedances over 90th quantile (from fitted Gaussian/Gamma marginal models),
# fit GPD to the exceedances

# pull 90th quantile from fitted marginal models for each station and season
q_mgcv_marg <- function(fit, newdata, family, prob = 0.9) {
  p <- predict(fit, newdata = newdata, type = "response")

  if (family == "gaulss") {
    mu <- p[, 1]
    sigma <- 1 / p[, 2]

    return(qnorm(prob, mean = mu, sd = sigma))
  }

  if (family == "gammals") {
    mu <- p[, 1]
    scale <- exp(p[, 2])

    shape <- 1 / scale
    rate <- shape / mu

    return(qgamma(prob, shape = shape, rate = rate))
  }

  stop("Unknown family")
}

# function to compute threshold and excesses for each station and season
thresh_mgcv <- function(
    data,
    mgcv_season_var,
    var,
    family,
    season,
    q = 0.9) {
  bind_rows(lapply(station_names, function(station) {
    dat <- data |>
      filter(name == station, season == !!season) |>
      arrange(date) |>
      filter(!is.na(.data[[var]]))

    if (family == "gammals") {
      dat <- dat |> filter(.data[[var]] > 0)
    }

    if (nrow(dat) == 0) {
      return(NULL)
    }

    # fit <- mgcv_season_var[[season]][[station]]$mgcv_fit
    fit <- mgcv_season_var[[season]]$mgcv_fit

    dat |>
      mutate(
        thresh = q_mgcv_marg(
          fit = fit,
          newdata = dat,
          family = family,
          prob = q
        ),
        excess = .data[[var]] - thresh
      ) |>
      filter(excess > 0) |>
      select(name, lon, lat, date, year, season_month, season, thresh, excess)
  }))
}

# threshold
thresh_data_season_var_mgcv_spat <- list(
  rain_Winter = thresh_mgcv(
    data = data,
    mgcv_season_var = gammals_season_var_spat,
    var = "rain",
    family = "gammals",
    season = "Winter",
    q = q
  ),
  rain_Summer = thresh_mgcv(
    data = data,
    mgcv_season_var = gammals_season_var_spat,
    var = "rain",
    family = "gammals",
    season = "Summer",
    q = q
  ),
  temp_Winter = thresh_mgcv(
    data = data,
    mgcv_season_var = gaulss_season_var_spat,
    var = "temp",
    family = "gaulss",
    season = "Winter",
    q = q
  ),
  temp_Summer = thresh_mgcv(
    data = data,
    mgcv_season_var = gaulss_season_var_spat,
    var = "temp",
    family = "gaulss",
    season = "Summer",
    q = q
  )
)


#### Fit GPD ####

find_min_k_pair_evgam <- function(
    x_spec,
    # k_spatial_vals = 3:9,
    k_spatial_vals = c(40),
    k_year_vals = 3:9,
    k_season_month_vals = 3:5,
    ...) {
  for (ks in k_spatial_vals) {
    for (ky in k_year_vals) {
      for (km in k_season_month_vals) {
        f_try <- list(
          excess ~ s(lon, lat, k = ks) + 
            s(year, k = ky) + 
            s(season_month, k = km),
          ~1
        )
  
        fit <- try(
          evgam::evgam(
            f_try,
            data = x_spec,
            family = "gpd",
            ...
          ),
          silent = TRUE
        )
  
        if (!inherits(fit, "try-error")) {
          return(list(
            k_spat = ks,
            k_year = ky,
            k_season_month = km,
            formula = f_try,
            fit = fit
          ))
        }
      }
    }
  }

  stop("No valid k pair found for evgam GPD.")
}

evgam_fit <- function(
    x,
    # name,
    ret_mod = FALSE,
    # k_spatial_vals = 3:9,
    k_spatial_vals = c(40),
    k_year_vals = 3:9,
    k_season_month_vals = 3:5,
    ...) {
  x_spec <- x |>
    select(name, lon, lat, year, season_month, thresh, excess) |>
    # filter(.data$name == .env$name) |>
    filter(
      !is.na(excess),
      excess > 0,
      !is.na(year),
      !is.na(season_month),
      !is.na(thresh)
    )

  if (nrow(x_spec) == 0) {
    return(NULL)
  }

  min_fit <- find_min_k_pair_evgam(
    x_spec = x_spec,
    k_spatial_vals = k_spatial_vals,
    k_year_vals = k_year_vals,
    k_season_month_vals = k_season_month_vals,
    ...
  )

  m <- min_fit$fit

  predictors <- m$predictor.names

  pred_dat_distinct <- x_spec |>
    # distinct(across(all_of(predictors)), thresh)
    distinct(name, across(all_of(predictors)), thresh)

  pred <- predict(m, pred_dat_distinct, type = "response")

  predictions <- bind_cols(
    pred_dat_distinct,
    as.data.frame(pred)
  ) |>
    mutate(
      # name = x_spec$name[1],
      k_year = min_fit$k_year,
      k_season_month = min_fit$k_season_month
    ) |>
    relocate(name, k_year, k_season_month)

  rownames(predictions) <- NULL

  if (ret_mod) {
    list(
      evgam_fit = m,
      marginal = predictions,
      formula = min_fit$formula,
      k_spatial = min_fit$k_spat,
      k_year = min_fit$k_year,
      k_season_month = min_fit$k_season_month
    )
  } else {
    predictions
  }
}

# fit GPD to the excesses for each station and season
gpd_mgcv_season_var_spat <- lapply(names(thresh_data_season_var_mgcv_spat), \(nm) {
  x <- thresh_data_season_var_mgcv_spat[[nm]]

  # ret <- mclapply(station_names, \(st) {
  #   evgam_fit(
  #     x = x,
  #     name = st,
  #     ret_mod = TRUE
  #   )
  # })
  # ret <- lapply(station_names, \(st) {
  ret <- evgam_fit(
      x = x,
      ret_mod = TRUE
    )
  # })
  # names(ret) <- station_names
  ret
})
names(gpd_mgcv_season_var_spat) <- names(thresh_data_season_var_mgcv_spat)


#### QQ plots ####

# fitted gpd cdf for excesses
pgpd_excess <- function(excess, scale, shape) {
  ifelse(
    abs(shape) < 1e-8,
    1 - exp(-excess / scale),
    1 - pmax(1 + shape * excess / scale, 0)^(-1 / shape)
  )
}

# function to compute QQ data for a given station and season
gpd_qq_data_one <- \(data_thresh, gpd_fit, station) {
  d <- data_thresh |>
    filter(name == station) |>
    select(name, lon, lat, date, year, season_month, thresh, excess)

  # gpd <- gpd_fit[[station]]$marginal |>
  #   select(year, season_month, thresh, scale, shape)
  gpd <- gpd_fit$marginal |>
    filter(name == station) |> 
    select(lon, lat, year, season_month, thresh, scale, shape)

  qq_dat <- d |>
    left_join(
      gpd # ,
      # by = c("year", "season_month", "thresh")
    ) |>
    filter(
      !is.na(excess),
      !is.na(scale),
      !is.na(shape),
      scale > 0
    ) |>
    mutate(
      g_hat = pgpd_excess(excess, scale, shape),
      g_hat = pmin(pmax(g_hat, 1e-10), 1 - 1e-10),
      exp_resid = -log(1 - g_hat)
    ) |>
    arrange(exp_resid) |>
    mutate(
      i = row_number(),
      n = n(),
      theoretical = qexp((i - 0.5) / n),
      observed = exp_resid,
      lower = qexp(qbeta(0.025, i, n - i + 1)),
      upper = qexp(qbeta(0.975, i, n - i + 1))
    )

  qq_dat
}

season_var_df <- tidyr::crossing(
  "var" = c("rain", "temp"),
  "season" = c("Winter", "Summer")
) |>
  mutate(lst_name = paste0(var, "_", season))


qq_df <- bind_rows(lapply(seq_len(nrow(season_var_df)), \(i) {
  lst_nm <- season_var_df$lst_name[[i]]

  # bind_rows(mclapply(station_names, \(st) {
  bind_rows(lapply(station_names, \(st) {
    # fit_obj <- gpd_mgcv_season_var_spat[[lst_nm]][[st]]
    fit_obj <- gpd_mgcv_season_var_spat[[lst_nm]]
    
    if (is.null(fit_obj)) {
      return(NULL)
    }

    gpd_qq_data_one(
      data_thresh = thresh_data_season_var_mgcv_spat[[lst_nm]],
      gpd_fit = gpd_mgcv_season_var_spat[[lst_nm]],
      station = st
    ) |>
      mutate(
        season = !!season_var_df$season[[i]],
        var = !!season_var_df$var[[i]]
      ) |>
      relocate(season, var, .after = name)
  }))
}))

qq_df |>
  filter(name == qq_df$name[1]) |>
  mutate(ind = paste0(var, ", ", season)) |>
  ggplot(aes(theoretical, observed)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_point() +
  facet_wrap(~ind) +
  cecl_theme() +
  labs(title = qq_df$name[1], x = "Theoretical", y = "Observed")

qq_plots <- mclapply(station_names, \(x) {
  qq_df |>
    filter(name == x) |>
    mutate(ind = paste0(var, ", ", season)) |>
    ggplot(aes(theoretical, observed)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point() +
    facet_wrap(~ind) +
    cecl_theme() +
    labs(title = x, x = "Theoretical", y = "Observed")
})

pdf("plots/02_app/mgcv_spat/05_qq_plots.pdf", width = 12, height = 8)
qq_plots
dev.off()


#### Piecewise Laplace transformation for GPD fitting ####

p_gpd_mgcv_ns <- function(dat, bulk_fit, gpd_fit, var, family) {
  x <- dat[[var]]

  # fitted bulk CDF at observed value
  F_bulk_x <- p_mgcv_spat_marg(
    dat = dat,
    fit = bulk_fit,
    var = var,
    family = family
  )

  # join GPD threshold/params
  dat_par <- dat |>
    left_join(
      gpd_fit |>
        select(name, lon, lat, year, season_month, thresh, scale, shape) |>
        distinct(name, lon, lat, year, season_month, .keep_all = TRUE) # ,
      # by = c("year", "season_month")
    )

  u <- dat_par$thresh

  # fitted bulk CDF at threshold
  dat_u <- dat_par
  dat_u[[var]] <- u

  F_bulk_u <- p_mgcv_spat_marg(
    dat = dat_u,
    fit = bulk_fit,
    var = var,
    family = family
  )

  F_hat <- F_bulk_x

  above <- !is.na(x) & !is.na(u) & x > u

  if (any(above)) {
    excess <- x[above] - u[above]
    sigma <- dat_par$scale[above]
    xi <- dat_par$shape[above]

    G <- ifelse(
      abs(xi) < 1e-8,
      1 - exp(-excess / sigma),
      1 - pmax(1 + xi * excess / sigma, 0)^(-1 / xi)
    )

    # splice GPD tail onto fitted bulk CDF
    F_hat[above] <- F_bulk_u[above] +
      (1 - F_bulk_u[above]) * G
  }

  pmin(pmax(F_hat, 1e-6), 1 - 1e-6)
}

trans_marg_mgcv_gpd_one <- function(
    data_orig,
    bulk_fit,
    gpd_fit,
    var,
    family,
    season,
    station) {
  dat <- data_orig |>
    select(name, lon, lat, date, year, season_month, season, all_of(var)) |>
    filter(name == station, season == !!season) |>
    arrange(date) |>
    filter(!is.na(.data[[var]]))

  if (family == "gammals") {
    dat <- dat |> filter(.data[[var]] > 0)
  }

  stopifnot(nrow(dat) > 0)

  # bulk_obj <- bulk_fit[[season]][[station]]$mgcv_fit
  # gpd_obj <- gpd_fit[[paste0(var, "_", season)]][[station]]$marginal
  bulk_obj <- bulk_fit[[season]]$mgcv_fit
  gpd_obj <- gpd_fit[[paste0(var, "_", season)]]$marginal

  F_hat <- p_gpd_mgcv_ns(
    dat = dat,
    bulk_fit = bulk_obj,
    gpd_fit = gpd_obj,
    var = var,
    family = family
  )

  dat |>
    mutate(
      F_hat = F_hat,
      laplace = dlaplace(F_hat)[, 1]
    )
}

laplace_mgcv_gpd_season_var_spat <- list(
  temp_Winter = bind_rows(lapply(station_names, \(st) {
    trans_marg_mgcv_gpd_one(data, gaulss_season_var_spat, gpd_mgcv_season_var_spat,
      var = "temp", family = "gaulss", season = "Winter", station = st
    )
  })) |> select(name, date, laplace),
  temp_Summer = bind_rows(lapply(station_names, \(st) {
    trans_marg_mgcv_gpd_one(data, gaulss_season_var_spat, gpd_mgcv_season_var_spat,
      var = "temp", family = "gaulss", season = "Summer", station = st
    )
  })) |> select(name, date, laplace),
  rain_Winter = bind_rows(lapply(station_names, \(st) {
    trans_marg_mgcv_gpd_one(data, gammals_season_var_spat, gpd_mgcv_season_var_spat,
      var = "rain", family = "gammals", season = "Winter", station = st
    )
  })) |> select(name, date, laplace),
  rain_Summer = bind_rows(lapply(station_names, \(st) {
    trans_marg_mgcv_gpd_one(data, gammals_season_var_spat, gpd_mgcv_season_var_spat,
      var = "rain", family = "gammals", season = "Summer", station = st
    )
  })) |> select(name, date, laplace)
)

trans_fun2 <- function(x, season) {
  laplace_df <- x[[paste0("temp_", season)]] |>
    rename(temp = laplace) |>
    left_join(
      x[[paste0("rain_", season)]] |>
        rename(rain = laplace),
      by = c("name", "date")
    ) |>
    select(name, temp, rain) |>
    mutate(name = factor(name))

  laplace_lst <- group_split(laplace_df, name, .keep = FALSE) |>
    lapply(\(x) {
      ret <- as.matrix(x)
      colnames(ret) <- names(x)
      ret
    })

  names(laplace_lst) <- unique(laplace_df$name)

  as_cecl_marg(laplace_lst)
}

# plot transformed data for a given location
laplace_mgcv_gpd_season_var_spat$rain_Summer |>
  filter(name == "Valencia") |>
  rename(rain = laplace) |>
  left_join(
    laplace_mgcv_gpd_season_var_spat$temp_Summer |>
      filter(name == "Valencia") |>
      rename(temp = laplace)
  ) |>
  ggplot(aes(x = temp, y = rain)) +
  geom_point()

marg_season_mgcv_gpd_spat <- list(
  Winter = trans_fun2(laplace_mgcv_gpd_season_var_spat, "Winter"),
  Summer = trans_fun2(laplace_mgcv_gpd_season_var_spat, "Summer")
)

saveRDS(
  marg_season_mgcv_gpd_spat,
  # file = "data/02_app/marg_season_mgcv_gpd.rds"
  file = "data/02_app/marg_season_mgcv_gpd_spat.rds"
)

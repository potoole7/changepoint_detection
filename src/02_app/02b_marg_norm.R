#### Fit marginal models to Spain data ####

# Change to fit to drought_local (will take Gaussian model) (done)
# Change to use cyclical spline for month (done)
# Transform to Laplace (done)

# TODO Change to some form of truncated normal distribution or something for fits the bulk well but does poorly at the tails where we have

# TODO Also fit GPD to upper tail of temperature data, but only join with bulk
# model for drought_local

# TODO Plot differences for laplace transformed data fitted using evgam and
# mgcv packages

# Idea: Fit marginal models to the Spain data for temperature and rainfall,
# using Gaussian and Gamma distributions respectively.
# The models will account for seasonal variations and temporal trends.


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

# different seasons and variables to fit margianl model for
season_var_df <- tidyr::crossing(
  "var" = c(var_dep, "temp"),
  "season" = seasons
) |>
  mutate(lst_name = paste0(var, "_", season))


#### Load Data ####

data <- readr::read_csv(
  "data/02_app/ecad_clean.csv.gz"
) |>
  # look at full years only
  filter(
    date >= as_date("1960-01-01"),
    date < as_date("2024-01-01"),
    # rain > 0 # only keep days with rain for rainfall modelling
    !is.na(!!var_dep)
  ) |>
  # add useful columns for modelling
  mutate(
    year = as.numeric(substr(date, 1, 4)),
    month = lubridate::month(date) #   ,
    #   s e a s o n _month = case_when(
    #   season == "Summer" ~ match(month, c(4, 5, 6, 7, 8, 9)),
    #   season == "Winter" ~ match(month, c(10, 11, 12, 1, 2, 3))
    # )
  ) |>
  relocate(c(year, month), .after = date)

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
    # date, year, season_month, season,
    date, year, season_year, month, season,
    temp, rain, wind_speed, contains("drought")
  ) |>
  filter(
    if_all(
      c(temp, rain, drought_local),
      \(x) !is.na(x)
    )
  ) |>
  mutate(month = lubridate::month(date))

# pull station names for looping over later
station_names <- unique(data$name)

# # add weights for drought_local to account for different number of
# # observations per station and season
# data <- data |>
#   group_by(name, season, season_year) |>
#   mutate(
#     drought_weight = 1 / n()
#   ) |>
#   ungroup()


#### Fit Gaussian and mma distributions ####

# function to find the minimum valid k pair for year and season_month
find_min_k_pair <- \(
  k_year_vals = 3:9,
  # k_season_month_vals = 3:5,
  response_name,
  x_spec,
  fam,
  ...
) {
  # w <- if (
  #   response_name == "drought_local" && "drought_weight" %in% names(x_spec)
  # ) {
  #   x_spec$drought_weight
  # } else {
  #   rep(1, nrow(x_spec))
  # }

  for (ky in k_year_vals) {
    # for (km in k_season_month_vals) {
    # if (response_name == "drought_local") {
    #   f_try <- list(
    #     reformulate(
    #       termlabels = c(
    #         sprintf("s(season_year, k = %s)", ky)
    #       ),
    #       response = response_name
    #     ),
    #     ~1
    #   )
    # } else {
    f_try <- list(
      reformulate(
        termlabels = c(
          sprintf("s(year, k = %s)", ky),
          sprintf("s(month, bs = 'cc', k = %s)", 11)
        ),
        response = response_name
      ),
      ~1
    )
    # }
    fit <- try(
      mgcv::gam(
        f_try,
        data = x_spec,
        family = fam,
        method = "REML",
        # weights = w,
        ...
      ),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      return(list(
        k_year = ky,
        # k_season_month = km,
        formula = f_try,
        fit = fit
      ))
    }
    # }
  }
  stop("No valid k pair found.")
}

# function to compute quantile-quantile data for a given fitted mgcv model and
# data
qq_data_mgcv <- function(m, x_spec, response_name, family) {
  y <- x_spec[[response_name]]
  p <- predict(m, x_spec, type = "response")

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

  tibble(
    observed = sort(qres),
    theoretical = qnorm(ppoints(length(qres)))
  )
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
  name,
  response_name = c("temp", "rain", "drought_local"),
  ret_mod = FALSE,
  family = c("gaulss", "gammals"),
  k_year_vals = 3:9,
  k_season_month_vals = 3:5,
  ...
) {
  family <- match.arg(family)
  response_name <- match.arg(response_name)

  # stopifnot(
  #   (response_name == "temp" && family == "gaulss") ||
  #     (response_name == "rain" && family == "gammals")
  # )

  x_spec <- x |>
    # select(name, year, season_month, all_of(response_name)) |>
    select(name, year, season_year, contains("month"), contains("weight"), all_of(response_name)) |>
    filter(name == !!name)

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
    k_year_vals = k_year_vals,
    # k_season_month_vals = k_season_month_vals,
    response_name = response_name,
    x_spec = x_spec,
    fam = fam,
    ...
  )

  m <- min_fit$fit
  f <- min_fit$formula

  predictors <- attr(terms(m), "term.labels")

  pred_dat_distinct <- x_spec |>
    distinct(across(all_of(predictors)))

  pred <- predict(m, pred_dat_distinct, type = "response")

  predictions <- bind_cols(
    pred_dat_distinct,
    as.data.frame(pred)
  )

  # names(predictions) <- c(predictors, "mu", "sigma")
  if (family == "gaulss") {
    names(predictions) <- c(predictors, "mu", "inv_sigma")
  }

  if (family == "gammals") {
    names(predictions) <- c(predictors, "mu", "log_scale")
  }

  predictions <- predictions |>
    mutate(
      name = x_spec$name[1],
      # k_season_month = min_fit$k_season_month
      k_year = min_fit$k_year
    ) |>
    # relocate(name, k_year, k_season_month)
    relocate(name, k_year)

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
      # k_season_month = min_fit$k_season_month,
      k_year = min_fit$k_year
    )
  } else {
    predictions
  }
}

# name <- "Badajoz Aeropuerto" # not great for drought_local
name <- "Valencia" # pretty good for both, might need a GPD upper tail

# test
plot(mgcv_fit(
  # x = data |> filter(season == "Winter"),
  x = data,
  name = name,
  response_name = "temp",
  family = "gaulss",
  ret_mod = TRUE
)$qq, type = "o")
abline(a = 0, b = 1, col = "red", lwd = 2)
# looks great!! Maybe some slight deviation at the tai that a GPD fited above
# a threshold could help with

plot(mgcv_fit(
  # x = data |> filter(season == "Winter"),
  x = data,
  name = name,
  response_name = var_dep,
  family = "gaulss",
  ret_mod = TRUE
)$qq, type = "o")
abline(a = 0, b = 1, col = "red", lwd = 2)

# fit gaulss distribution to temperature and drought_local for each station
gaulss_var <- lapply(c("temp", var_dep), \(var) {
  ret <- mclapply(station_names, \(station) {
    # if (var == "temp") {
    #   data_var <- data
    # } else if (var == var_dep) {
    #   data_var <- data_drought
    # } else {
    #   stop("Unknown variable")
    # }
    mgcv_fit(
      data,
      # data_var,
      name = station,
      response_name = var,
      family = "gaulss",
      ret_mod = TRUE
    )
  })
  names(ret) <- station_names
  ret
})
names(gaulss_var) <- c("temp", var_dep)


#### Validate ####

spec <- gaulss_var[[1]][[1]]
pred <- predict(spec$mgcv_fit, type = "response")

str(pred)
head(pred)
summary(pred[, 1])
summary(pred[, 2])

# function to compute QQ envelope for a given sample size n,
# number of simulations nsim, and confidence level conf
# TODO NOTE: QQ envelope assumes stationarity, which is not the case here,
# so this is just a rough guide
qq_envelope <- function(n, nsim = 1000, conf = 0.95) {
  sims <- replicate(nsim, sort(rnorm(n)))

  tibble(
    theoretical = qnorm(ppoints(n)),
    lower = apply(sims, 1, quantile, probs = (1 - conf) / 2),
    upper = apply(sims, 1, quantile, probs = 1 - (1 - conf) / 2)
  )
}

# pull QQ data for each sesaon and station
# qq_data <- bind_rows(lapply(seq_along(gaulss_var), \(i) {
#   season_name <- names(gaulss_var)[i]
#   bind_rows(lapply(seq_along(station_names), \(j) {
#     bind_rows(
#       mutate(gaulss_var[[i]][[j]]$qq, family = "gaulss"),
#       mutate(gammals_season_var[[i]][[j]]$qq, family = "gammals")
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
qq_data <- bind_rows(lapply(names(gaulss_var), \(var) {
  bind_rows(lapply(station_names, \(station) {
    mutate(
      gaulss_var[[var]][[station]]$qq,
      family = "gaulss",
      var = !!var,
      name = !!station
    ) |>
      relocate(name, var, family)
  }))
})) |>
  # calculate uncertainty envelopes for each station and season
  group_by(name, var, family) |>
  group_modify(\(.x, .y) {
    bind_cols(
      .x |> arrange(theoretical),
      qq_envelope(nrow(.x)) |> select(lower, upper)
    )
  }) |>
  ungroup()

# plot for each station and save results
pdf("plots/02_app/mgcv/marginal_qq_plots.pdf", width = 10, height = 8)
lapply(station_names, \(station) {
  qq_data |>
    filter(name == !!station) |>
    # mutate(
    #   ind = paste0(family, ", ", season),
    #   ind = factor(ind, levels = c(
    #     "gaulss, Winter", "gaulss, Summer",
    #     "gammals, Winter", "gammals, Summer"
    #   ))
    # ) |>
    ggplot(aes(x = theoretical, y = observed)) +
    geom_point() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    # facet_wrap(~ind, scales = "free") +
    facet_wrap(~var, scales = "free") +
    cecl_theme() +
    labs(title = station)
})
dev.off()

# conclusions:
# Across stations, QQ residual plots from the gaulss model for temperature are
# generally close to the theoretical normal distribution, with only minor
# departures in the extreme upper tail (particularly A Coruna Badajoz, Almeria)
# Fit really good to SPI/drought local! Few deviations but hopefully fixed by
# fitting GPD to upper tail

#### Laplace transformation ####

p_mgcv_marg <- \(dat, fit, var, family) {
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

trans_mgcv_one <- \(
  data_orig,
  mgcv_fit,
  var,
  family,
  season,
  station) {
  dat <- data_orig |>
    # select(name, lon, lat, date, year, season_month, season, all_of(var)) |>
    select(name, lon, lat, date, season_year, year, month, season, all_of(var)) |>
    filter(name == !!station, season == !!season) |>
    arrange(date) |>
    filter(!is.na(.data[[var]]))

  if (family == "gammals") {
    dat <- dat |> filter(.data[[var]] > 0)
  }

  stopifnot(nrow(dat) > 0)

  fit_obj <- mgcv_fit[[station]]

  F_hat <- p_mgcv_marg(
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

# for (s in names(gaulss_var)) {
#   names(gaulss_var[[s]]) <- station_names
#   names(gammals_season_var[[s]]) <- station_names
# }

# laplace_mgcv_season_var <- list(
#   temp_Winter = bind_rows(mclapply(station_names, \(st) {
#     trans_mgcv_one(
#       data_orig = data,
#       mgcv_fit = gaulss_var$Winter,
#       var = "temp",
#       family = "gaulss",
#       season = "Winter",
#       station = st
#     )
#   })) |> select(name, laplace),
#   temp_Summer = bind_rows(mclapply(station_names, \(st) {
#     trans_mgcv_one(
#       data_orig = data,
#       mgcv_fit = gaulss_var$Summer,
#       var = "temp",
#       family = "gaulss",
#       season = "Summer",
#       station = st
#     )
#   })) |> select(name, laplace),
#   rain_Winter = bind_rows(mclapply(station_names, \(st) {
#     trans_mgcv_one(
#       data_orig = data,
#       mgcv_fit = gammals_season_var$Winter,
#       var = "rain",
#       family = "gammals",
#       season = "Winter",
#       station = st
#     )
#   })) |> select(name, laplace),
#   rain_Summer = bind_rows(mclapply(station_names, \(st) {
#     trans_mgcv_one(
#       data_orig = data,
#       mgcv_fit = gammals_season_var$Summer,
#       var = "rain",
#       family = "gammals",
#       season = "Summer",
#       station = st
#     )
#   })) |> select(name, laplace)
# )
laplace_mgcv_season_var <- lapply(seq_len(nrow(season_var_df)), \(i) {
  bind_rows(mclapply(station_names, \(station) {
    with(
      season_var_df,
      trans_mgcv_one(
        data_orig = data,
        # mgcv_fit = gammals_season_var$Summer,
        mgcv_fit = gaulss_var[[var[[i]]]],
        var = var[[i]],
        family = "gaulss",
        # season = "Summer",
        season = season[[i]],
        station = station
      )
    )
  })) %>%
    select(name, laplace)
})
names(laplace_mgcv_season_var) <- season_var_df$lst_name

# convert to cecl_marg structure/object
trans_fun <- \(x, season) {
  laplace_df <- cbind(
    x[[paste0("temp_", season)]] |>
      rename(temp = laplace),
    # x[[paste0("rain_", season)]] |>
    x[[paste0(var_dep, "_", season)]] |>
      select(-name) |>
      # TODO Make more programatic
      # rename(rain = laplace)
      rename(drought_local = laplace)
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

# marg_season_mgcv <- list(
#   Winter = trans_fun(laplace_mgcv_season_var, "Winter"),
#   Summer = trans_fun(laplace_mgcv_season_var, "Summer")
# )
marg_season_mgcv <- lapply(seasons, \(season) {
  trans_fun(laplace_mgcv_season_var, season)
})
names(marg_season_mgcv) <- seasons

saveRDS(
  marg_season_mgcv,
  file = "data/02_app/marg_season_mgcv.rds"
)

# plot one transformed dataset
data.frame(marg_season_mgcv$Winter$transformed$`Badajoz Aeropuerto`) |>
  # data.frame(marg_season_mgcv$Winter$transformed$Santander) |>
  # ggplot(aes(x = temp, y = drought_local)) +
  ggplot(aes(x = temp, y = get(var_dep))) +
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
    mgcv_var,
    var,
    family,
    # season,
    q = 0.9) {
  bind_rows(lapply(station_names, function(station) {
    dat <- data |>
      # filter(name == station, season == !!season) |>
      filter(name == station) |>
      arrange(date) |>
      filter(!is.na(.data[[var]]))

    if (family == "gammals") {
      dat <- dat |> filter(.data[[var]] > 0)
    }

    if (nrow(dat) == 0) {
      return(NULL)
    }

    # fit <- mgcv_season_var[[season]][[station]]$mgcv_fit
    fit <- mgcv_var[[var]][[station]]$mgcv_fit

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
      # select(name, lon, lat, date, year, season_month, season, thresh, excess)
      select(name, lon, lat, date, year, month, season, thresh, excess)
  }))
}

# threshold
# thresh_data_season_var_mgcv <- list(
#   rain_Winter = thresh_mgcv(
#     data = data,
#     mgcv_season_var = gammals_season_var,
#     var = "rain",
#     family = "gammals",
#     season = "Winter",
#     q = q
#   ),
#   rain_Summer = thresh_mgcv(
#     data = data,
#     mgcv_season_var = gammals_season_var,
#     var = "rain",
#     family = "gammals",
#     season = "Summer",
#     q = q
#   ),
#   temp_Winter = thresh_mgcv(
#     data = data,
#     mgcv_season_var = gaulss_var,
#     var = "temp",
#     family = "gaulss",
#     season = "Winter",
#     q = q
#   ),
#   temp_Summer = thresh_mgcv(
#     data = data,
#     mgcv_season_var = gaulss_var,
#     var = "temp",
#     family = "gaulss",
#     season = "Summer",
#     q = q
#   )
# )

thresh_data_var <- lapply(c("temp", var_dep), \(var) {
  thresh_mgcv(
    data = data,
    mgcv_var = gaulss_var,
    var = var,
    family = "gaulss",
    q = q
    # q = 0.85
  )
})
names(thresh_data_var) <- c("temp", var_dep)


#### Fit GPD ####

find_min_k_pair_evgam <- \(
  x_spec,
  # k_year_vals = 3:9,
  k_year_vals = 1:9,
  # k_season_month_vals = 3:5,
  ...) {
  for (ky in k_year_vals) {
    # for (km in k_season_month_vals) {
    f_try <- list(
      # excess ~ s(year, k = ky) + s(season_month, k = km),
      excess ~ s(year, k = ky) + s(month, bs = "cc", k = 11),
      ~1
    )

    fit <- try(
      evgam::evgam(
        formula = f_try,
        data = x_spec,
        family = "gpd",
        ...
      ),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {
      return(list(
        k_year = ky,
        # k_season_month = km,
        formula = f_try,
        fit = fit
      ))
    }
    # }
  }

  stop("No valid k pair found for evgam GPD.")
}

evgam_fit <- \(
  x,
  name,
  ret_mod = FALSE,
  k_year_vals = 2:9,
  # k_year_vals = 1:10,
  # k_season_month_vals = 3:5,
  ...) {
  x_spec <- x |>
    # select(name, year, season_month, thresh, excess) |>
    select(name, year, month, season, thresh, excess) |>
    filter(.data$name == .env$name) |>
    filter(
      !is.na(excess),
      excess > 0,
      !is.na(year),
      # !is.na(season_month),
      !is.na(month),
      !is.na(thresh)
    )

  if (nrow(x_spec) == 0) {
    return(NULL)
  }

  min_fit <- find_min_k_pair_evgam(
    x_spec = x_spec,
    k_year_vals = k_year_vals,
    # k_season_month_vals = k_season_month_vals,
    ...
  )

  m <- min_fit$fit

  predictors <- m$predictor.names

  pred_dat_distinct <- x_spec |>
    distinct(across(all_of(predictors)), season, thresh)

  pred <- predict(m, pred_dat_distinct, type = "response")

  predictions <- bind_cols(
    pred_dat_distinct,
    as.data.frame(pred)
  ) |>
    mutate(
      name = x_spec$name[1],
      k_year = min_fit$k_year,
      # k_season_month = min_fit$k_season_month
    ) |>
    # relocate(name, k_year, k_season_month)
    relocate(name, k_year)

  rownames(predictions) <- NULL

  if (ret_mod) {
    list(
      evgam_fit = m,
      marginal = predictions,
      formula = min_fit$formula,
      # k_season_month = min_fit$k_season_month,
      k_year = min_fit$k_year
    )
  } else {
    predictions
  }
}

# fit GPD to the excesses for each station and season
# gpd_mgcv_season_var <- lapply(names(thresh_data_season_var_mgcv), \(nm) {
gpd_mgcv_var <- lapply(names(thresh_data_var), \(var) {
  # x <- thresh_data_season_var_mgcv[[nm]]
  print(var)
  x <- thresh_data_var[[var]]

  ret <- lapply(station_names, \(station) {
    print(station)
    evgam_fit(
      x = x,
      name = station,
      ret_mod = TRUE
    )
  })

  names(ret) <- station_names
  ret
})
# names(gpd_mgcv_season_var) <- names(thresh_data_season_var_mgcv)
names(gpd_mgcv_var) <- names(thresh_data_var)


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
    # select(name, date, year, season_month, thresh, excess)
    select(name, date, year, month, season, thresh, excess)

  gpd <- gpd_fit[[station]]$marginal |>
    # select(year, season_month, thresh, scale, shape)
    select(year, month, season, thresh, scale, shape)

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
  # "var" = c("rain", "temp"),
  "var" = c(var_dep, "temp"),
  # "season" = c("Winter", "Summer")
  season = seasons
) |>
  mutate(lst_name = paste0(var, "_", season))

# qq_df <- bind_rows(lapply(seq_len(nrow(season_var_df)), \(i) {
#   lst_nm <- season_var_df$lst_name[[i]]
#
#   # bind_rows(mclapply(station_names, \(st) {
#   bind_rows(lapply(station_names, \(st) {
#     # fit_obj <- gpd_mgcv_season_var[[lst_nm]][[st]]
#     fit_obj <- gpd_mgcv_season_var[[lst_nm]][[st]]
#
#     if (is.null(fit_obj)) {
#       return(NULL)
#     }
#
#     gpd_qq_data_one(
#       data_thresh = thresh_data_season_var_mgcv[[lst_nm]],
#       gpd_fit = gpd_mgcv_season_var[[lst_nm]],
#       station = st
#     ) |>
#       mutate(
#         season = !!season_var_df$season[[i]],
#         var = !!season_var_df$var[[i]]
#       ) |>
#       relocate(season, var, .after = name)
#   }))
# }))

qq_df <- bind_rows(lapply(names(thresh_data_var), \(var) {
  bind_rows(mclapply(station_names, \(station) {
    # bind_rows(lapply(station_names, \(station) {
    fit_obj <- gpd_mgcv_var[[var]][[station]]

    if (is.null(fit_obj)) {
      return(NULL)
    }

    gpd_qq_data_one(
      data_thresh = thresh_data_var[[var]],
      gpd_fit     = gpd_mgcv_var[[var]],
      station     = station
    ) |>
      mutate(var = !!var) |>
      relocate(season, var, .after = name)
  }))
}))

# test
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
  facet_wrap(~ind, nrow = 2, scales = "free") +
  cecl_theme() +
  labs(title = qq_df$name[1], x = "Theoretical", y = "Observed")
# looks great!!

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
    facet_wrap(~ind, nrow = 2, scales = "free") +
    cecl_theme() +
    labs(title = x, x = "Theoretical", y = "Observed")
})

pdf("plots/02_app/mgcv/05_qq_plots.pdf", width = 12, height = 8)
qq_plots
dev.off()


#### Piecewise Laplace transformation for GPD fitting ####

p_gpd_mgcv_ns <- \(dat, bulk_fit, gpd_fit, var, family) {
  x <- dat[[var]]

  # fitted bulk CDF at observed value
  F_bulk_x <- p_mgcv_marg(
    dat = dat,
    fit = bulk_fit,
    var = var,
    family = family
  )

  # join GPD threshold/params
  dat_par <- dat |>
    left_join(
      gpd_fit |>
        # select(year, season_month, thresh, scale, shape) |>
        # distinct(year, season_month, .keep_all = TRUE) # ,
        select(year, month, season, thresh, scale, shape) |>
        distinct(year, month, season, .keep_all = TRUE)
      # by = c("year", "season_month")
    )

  u <- dat_par$thresh

  # fitted bulk CDF at threshold
  dat_u <- dat
  dat_u[[var]] <- u

  F_bulk_u <- p_mgcv_marg(
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

trans_marg_mgcv_gpd_one <- \(
  data_orig,
  bulk_fit,
  gpd_fit,
  var,
  family,
  season,
  station) {
  dat <- data_orig |>
    # select(name, lon, lat, date, year, season_month, season, all_of(var)) |>
    select(name, lon, lat, date, year, month, season, all_of(var)) |>
    filter(name == station, season == !!season) |>
    arrange(date) |>
    filter(!is.na(.data[[var]]))

  if (family == "gammals") {
    dat <- dat |> filter(.data[[var]] > 0)
  }

  stopifnot(nrow(dat) > 0)

  # bulk_obj <- bulk_fit[[season]][[station]]$mgcv_fit
  # gpd_obj <- gpd_fit[[paste0(var, "_", season)]][[station]]$marginal
  bulk_obj <- bulk_fit[[var]][[station]]$mgcv_fit
  gpd_obj <- gpd_fit[[var]][[station]]$marginal

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

# # transform data for each station, season and variable
# laplace_mgcv_gpd_season_var <- list(
#   temp_Winter = bind_rows(lapply(station_names, \(st) {
#     trans_marg_mgcv_gpd_one(data, gaulss_season_var, gpd_mgcv_season_var,
#       var = "temp", family = "gaulss", season = "Winter", station = st
#     )
#   })) |> select(name, date, laplace),
#   temp_Summer = bind_rows(lapply(station_names, \(st) {
#     trans_marg_mgcv_gpd_one(data, gaulss_season_var, gpd_mgcv_season_var,
#       var = "temp", family = "gaulss", season = "Summer", station = st
#     )
#   })) |> select(name, date, laplace),
#   rain_Winter = bind_rows(lapply(station_names, \(st) {
#     trans_marg_mgcv_gpd_one(data, gammals_season_var, gpd_mgcv_season_var,
#       # var = "rain", family = "gammals", season = "Winter", station = st
#       var = var_dep, family = "gammals", season = "Winter", station = st
#     )
#   })) |> select(name, date, laplace),
#   rain_Summer = bind_rows(lapply(station_names, \(st) {
#     trans_marg_mgcv_gpd_one(data, gammals_season_var, gpd_mgcv_season_var,
#       # var = "rain", family = "gammals", season = "Summer", station = st
#       var = var_dep, family = "gammals", season = "Summer", station = st
#     )
#   })) |> select(name, date, laplace)
# )

# transform data for each station, season and variable
laplace_mgcv_gpd_season_var <- lapply(seq_len(nrow(season_var_df)), \(i) {
  bind_rows(mclapply(station_names, \(station) {
    with(
      season_var_df,
      trans_marg_mgcv_gpd_one(
        data, gaulss_var, gpd_mgcv_var,
        var = var[i],
        family = "gaulss",
        # season = "Winter",
        season = season[i],
        station = station
      )
    )
  })) |>
    select(name, date, laplace)
})
names(laplace_mgcv_gpd_season_var) <- season_var_df$lst_name

trans_fun2 <- \(x, season, dep_var) {
  laplace_df <- x[[paste0("temp_", season)]] |>
    rename(temp = laplace) |>
    left_join(
      # x[[paste0("rain_", season)]] |>
      x[[paste0(dep_var, "_", season)]] |>
        # rename(rain = laplace),
        rename(!!dep_var := laplace),
      by = c("name", "date")
    ) |>
    # select(name, temp, rain) |>
    select(name, temp, !!dep_var) |>
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
# laplace_mgcv_gpd_season_var$rain_Summer |>
laplace_mgcv_gpd_season_var[[paste0(var_dep, "_Summer")]] |>
  filter(name == "Valencia") |>
  rename(!!var_dep := laplace) |>
  left_join(
    laplace_mgcv_gpd_season_var$temp_Summer |>
      filter(name == "Valencia") |>
      rename(temp = laplace)
  ) |>
  # ggplot(aes(x = temp, y = rain)) +
  ggplot(aes(x = temp, y = get(var_dep))) +
  geom_point()

# marg_season_mgcv_gpd <- list(
#   Winter = trans_fun2(laplace_mgcv_gpd_season_var, "Winter"),
#   Summer = trans_fun2(laplace_mgcv_gpd_season_var, "Summer")
# )
marg_season_mgcv_gpd <- lapply(seasons, \(season) {
  print(season)
  trans_fun2(laplace_mgcv_gpd_season_var, season, var_dep)
})

saveRDS(
  marg_season_mgcv_gpd,
  file = "data/02_app/marg_season_mgcv_gpd.rds"
)

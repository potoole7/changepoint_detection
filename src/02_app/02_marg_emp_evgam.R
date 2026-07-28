#### Fit marginal models to Spain data ####

# TODO Allow for location-specific quantiles in CeCl

# TODO Fit for drought_local vs temp

# Threshold with qgam model fit to all months together with# cyclical spline (don't split) (done)
# Fit evgam model in the same way (done)
# TODO QQ plots
# TODO Use method to include lowest number of k possible from other scripts
# TODO What to do about empirical transformation in this case???


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

seed <- 123

decades <- seq(1960, 2010, by = 10)

seasons <- c("Winter", "Spring", "Summer", "Autumn")

# temp_var <- "temp"
temp_var <- "temp_max"

# var_dep <- "rain"
# var_dep <- "drought_local"
var_dep <- "drought_local_rev"

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
    month = lubridate::month(date) # ,
    # season_month = case_when(
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

# remove 0s in rain if interested in rain, since they are not extreme events
if (var_dep == "rain") {
  data <- data |>
    filter(rain > 0)
}

# variables of interest
vars <- c("temp", var_dep)

# different seasons and variables to fit margianl model for
season_var_df <- tidyr::crossing(
  # "var" = c("rain", "temp"),
  "var" = c(var_dep, "temp"),
  # "season" = c("Winter", "Summer")
  "season" = c("Winter", "Spring", "Summer", "Autumn")
) |>
  mutate(lst_name = paste0(var, "_", season))


areas <- read_sf("data/02_app/spain_shapefile.geojson") |>
  filter(
    !ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla")
  )

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")


#### Initial Calculations ####

data2 <- data |>
  select(
    name = station_name, lon, lat,
    # date, year, season_month, season,
    date, year, season_year, month, season,
    temp, rain, wind_speed, contains("drought")
  ) |>
  filter(
    if_all(
      # c(temp, rain, drought_local, drought_global),
      c(temp, rain, drought_local),
      \(x) !is.na(x)
    )
  )

# pull station names for looping over later
station_names <- unique(data2$name)

data2 |>
  filter(name == data2$name[1]) |>
  ggplot(aes(date, drought_local)) +
  geom_point(alpha = 0.7)

# data_drought |>
#   filter(name == data2$name[1]) |>
#   ggplot(aes(season_year, drought_local)) +
#   geom_point(alpha = 0.7)

#### Convert temperature to seasonal max ####

# data2 <- data2 |>
#   group_by(name, lon, lat, season, season_year) |>
#   summarise(
#     temp = max(temp, na.rm = TRUE),
#     drought_local = max(drought_local, na.rm = TRUE),
#     .groups = "drop"
#   )


#### Empirical transformation ####

# Split marginal models by season
data2 <- data2 |>
  # bind_rows(
  #   data2 |>
  #     mutate(season = "Year")
  # ) |>
  # mutate(season = factor(season, levels = c("Winter", "Summer", "Year"))) |>
  mutate(
    season = factor(
      # season, levels = c("Winter", "Spring", "Summer", "Autumn", "Year")
      season,
      # levels = c("Winter", "Spring", "Summer", "Autumn")
      levels = seasons
    ),
    decade = factor(floor(year(date) / 10) * 10, levels = decades)
  )

# store for later
seasons_all <- levels(data2$season)

marg_season <- data2 |>
  group_split(season) |>
  mclapply(\(x)   {
    # lapply(\(x) {
    # browser()
    # print(x$name[1])
    # print(x$season[1])
    CeCl::cecl_marg(
      # select(x, -rain),
      # select(x, name, temp, rain),
      select(x, name, temp, !!sym(var_dep)),
      mult_col = "name",
      # vars = c("temp", "rain", "wind_speed", "drought_local", "drought_global"),
      # vars = c("temp", "rain", "drought_local", "drought_global"),
      # vars = c("temp", "rain"),
      vars = c("temp", var_dep),
      thresh_method = "none",
      marg_method = "ecdf"
    )
  })
# names(marg_season) <- c("Winter", "Summer", "Year")
names(marg_season) <- seasons_all

saveRDS(
  marg_season,
  file = "data/02_app/marg_season_emp.rds"
)

# also do for each decade
marg_season_decade <- data2 |>
  # bind_rows(
  #   data2 |>
  #     mutate(season = "Year")
  # ) |>
  filter(!is.na(decade)) |>
  group_split(decade, season) |>
  mclapply(\(x) {
    ret <- cecl_marg(
      x,
      mult_col = "name",
      # vars = c("temp", "rain", "drought_local", "drought_global"),
      vars = c("temp", var_dep),
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
  # "season" = c("Winter", "Summer", "Year")
  "season" = seasons_all
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
    # select(name, lon, lat, date, year, season_month, season, thresh, excess)
    select(name, lon, lat, date, year, month, season, thresh, excess)

  # names(ret)[names(ret) == "excess"] <- var
  return(ret)
}

# test
# thresh_data(data2, "rain", "Winter", q = 0.9)
thresh_data(data2, var_dep, "Winter", q = 0.9)


# apply thresholding function to each season and variable combination
thresh_data_season_var <- season_var_df |>
  mutate(var = as.factor(var), season = as.factor(season)) |>
  group_split(var, season) |>
  lapply(\(x) {
    with(x, thresh_data(data2, as.character(var), as.character(season)))
  })

names(thresh_data_season_var) <- season_var_df$lst_name


#### Vary threshold with qgam ####

q <- 0.9

# function to threshold data using qgam, for given target threshold
thresh_qgam <- \(
  data,
  # TODO Fit with same method as before to use as small a k as possible
  # f = target ~ s(year) + s(season_month, k = 5),
  f = target ~ s(year, k = 3) + s(month, k = 11, bs = "cc"),
  var,
  q = 0.9,
  ret_mod = FALSE
) {
  data_spec <- data

  names(data_spec)[names(data_spec) == var] <- "target"

  fit <- qgam::qgam(f, data = data_spec, qu = 0.9)

  # summary(fit)
  predictions <- predict(fit, data = data_spec, type = "response")
  # summary(predictions)
  # plot(fit, pages = 1)

  ret <- data_spec |>
    mutate(
      thresh = predictions,
      excess = target - thresh
    ) |>
    filter(excess > 0) |>
    # select(name, lon, lat, date, year, season_month, season, thresh, excess)
    select(name, lon, lat, date, year, month, season, thresh, excess)

  # ret |>
  #   select(thresh, date, year, season_month) |>
  #   distinct() |>
  #   ggplot(aes(
  #     x = date, y = predictions, colour = season_month
  #     # x = date, y = predictions, colour = year
  #   )) +
  #   # geom_point() +
  #   geom_line() +
  #   cecl_theme(nejm_pal = FALSE) +
  #   # theme(legend.title = element_blank(), legend.text = element_text(angle = 45, hjust = 1)) +
  #   NULL

  # data_spec |>
  #   mutate(thresh = predictions) |>
  #   ggplot(aes(date)) +
  #   geom_point(aes(y = target), alpha = 0.2) +
  #   geom_line(aes(y = thresh), linewidth = 1) +
  #   cecl_theme()

  if (ret_mod) {
    return(list(
      "qgam_fit" = fit,
      "thresh"   = ret
    ))
  } else {
    return(ret)
  }
}

# TODO Change so that it no longer loops through seasons
# fit location-specific GPDs for each variable-season setup
# i <- 1
# y <- station_names[1]
# thresh_data_season_var_qgam <- lapply(seq_len(nrow(season_var_df)), \(i) {
#   print(season_var_df$lst_name[[i]])
#   x <- data2 |>
#     filter(season == season_var_df$season[[i]])
#   bind_rows(lapply(station_names, \(y) {
#     print(y)
#     data_spec <- filter(x, name == y)
#     n_month <- length(unique(data_spec$season_month))
#     k_month <- min(5, n_month - 1) # must be greater than number of months aail
#
#     f <- as.formula(
#       paste0("target ~ s(year) + s(season_month, k = ", k_month, ")")
#    )
#
#     thresh_qgam(
#       data_spec,
#       f,
#       var = season_var_df$var[[i]],
#       q = q
#     )
#   }))
# })
# names(thresh_data_season_var_qgam) <- season_var_df$season
# names(thresh_data_season_var_qgam) <- season_var_df$lst_name
thresh_data_var_qgam <- lapply(vars, \(var) {
  print(var)
  bind_rows(mclapply(station_names, \(station) {
    print(station)
    data_spec <- filter(data2, name == station)

    f <- target ~ s(year, k = 3) + s(month, bs = "cc", k = 11)

    thresh_qgam(
      data_spec,
      f,
      var = var,
      q = q
    )
  }))
})
names(thresh_data_var_qgam) <- vars

# check
thresh_data_var_qgam[[1]] |>
  filter(name == "Badajoz Aeropuerto", year == 1961)


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
  # f = list(excess ~ s(year) + s(season_month, k = 5), ~1),
  f = list(excess ~ s(year) + s(month, k = 11, bs = "cc"), ~1),
  name,
  ret_mod = FALSE
) {
  # pull useful names
  x_spec <- x |>
    # select(name, year, season_month, thresh, excess) |>
    select(name, year, month, season, thresh, excess) |>
    filter(name == !!name)

  # formula
  # f <- list(excess ~ year, ~year)

  # fit evgam model
  m <- evgam::evgam(f, data = x_spec, family = "gpd")

  # create predictions for unique rows in pred_data (ensures one pred per loc)
  predictors <- m$predictor.names
  pred_dat_distinct <- x_spec |>
    # keep season to label properly
    # distinct(across(all_of(predictors)), thresh)
    distinct(across(all_of(predictors)), season, thresh)
  # distinct(across(all_of(predictors)), thresh)
  predictions <- cbind(
    pred_dat_distinct,
    predict(m, pred_dat_distinct, type = "response")
  )

  # TODO Also return uncertainty!!
  predictions <- predictions |>
    mutate(name = x_spec$name[1]) |>
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

# # test for Valencia (and other locations)
# TODO Accidentally used
obj1 <- evgam_fit(
  # thresh_data_season_var$temp_Winter,
  thresh_data_var_qgam$temp,
  # f = list(excess ~ s(year) + s(season_month, k = 3), ~1),
  f = list(excess ~ s(year, k = 3) + s(month, k = 11, bs = "cc"), ~1),
  name = name,
  ret_mod = TRUE
)
fit1 <- obj1[[1]]
summary(fit1)
plot(fit1, pages = 1)

obj2 <- evgam_fit(
  # thresh_data_season_var$temp_Winter,
  # thresh_data_season_var$temp,
  thresh_data_var_qgam$temp,
  f = list(excess ~ s(year, k = 3) + s(month, k = 11, bs = "cc"), ~ s(year)),
  name = name,
  ret_mod = TRUE
)
fit2 <- obj2[[1]]
summary(fit2)
plot(fit2, pages = 1)

AIC(fit1, fit2) #

# plot
bind_rows(
  mutate(obj1$marginal, type = "fixed_shape"),
  mutate(obj2$marginal, type = "varying_shape")
) |>
  pivot_longer(c(scale, shape), names_to = "parameter") |>
  ggplot(aes(x = year, y = value, colour = type)) +
  geom_line() +
  facet_wrap(~parameter, scale = "free") +
  cecl_theme()

# var <- "drought_local"
# # station <- "Getafe"
# station <- "Barcelona Fabra"

# fit location-specific GPDs for each variable-season setup
set.seed(seed)
gpd_var <- lapply(vars, \(var) {
  x <- thresh_data_var_qgam[[var]]
  print(var)
  # bind_rows(lapply(station_names, \(station) {
  bind_rows(mclapply(station_names, \(station) {
    print(station)
    dat <- x |>
      filter(name == station) |>
      # ensure distinct works on thresh
      mutate(thresh = round(thresh, 4))
    # n_month <- length(unique(dat$season_month))
    # k_month <- min(5, n_month - 1) # must be greater than number of months avail

    # TODO Use minimum k allowed
    f <- list(
      as.formula(paste0(
        # "excess ~ s(year) + s(season_month, k = ", k_month, ")"
        "excess ~ s(year, k = 3) + s(month, k = 11, bs = 'cc')"
      )),
      ~1
    )

    ret <- evgam_fit(dat, f = f, name = station)
  }))
})
names(gpd_var) <- names(thresh_data_var_qgam)


# TODO Plot time series of scale parameter
# TODO Plot maps for some years (end of decade ones)


#### Laplace transform ####

# function to perform marginal transformation
trans_marg_ns_one <- \(data_orig,
  gpd_fit,
  var,
  season,
  station) {
  # take data for given station and season
  dat <- data_orig |>
    # select(name, lon, lat, date, year, season_month, season, all_of(var)) |>
    select(name, lon, lat, date, year, month, season, all_of(var)) |>
    filter(name == !!station, season == !!season) |>
    arrange(date) |>
    filter(!is.na(.data[[var]]))

  stopifnot(nrow(dat) > 0)

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

# # Semi-parametric CDF with year-varying GPD parameters
# p_gpd_ecdf_ns <- \(dat, gpd, var) {
#   x <- dat[[var]]
#   u <- gpd$thresh[1]
#
#   # Join year-varying parameters
#   dat_par <- dat |>
#     mutate(year = as.numeric(substr(date, 1, 4))) |>
#     left_join(
#       gpd |>
#         select(year, scale, shape),
#       by = "year"
#     )
#
#   # Empirical CDF for all observations
#   ord <- order(x)
#   x_sort <- x[ord]
#   m <- length(x)
#
#   ecdf_vals <- seq_len(m) / (m + 1)
#   ecdf_orig <- numeric(m)
#   ecdf_orig[ord] <- ecdf_vals
#
#   # Initialise
#   F_hat <- ecdf_orig
#
#   # Parametric GPD tail above threshold
#   above <- x > u
#
#   if (any(above, na.rm = TRUE)) {
#     excess <- x[above] - u
#     sigma <- dat_par$scale[above]
#     xi <- dat_par$shape[above]
#
#     surv_gpd <- ifelse(
#       abs(xi) < 1e-8,
#       exp(-excess / sigma),
#       pmax(1 + xi * excess / sigma, 0)^(-1 / xi)
#     )
#
#     p_exc <- mean(x > u, na.rm = TRUE)
#
#     F_hat[above] <- 1 - p_exc * surv_gpd
#   }
#
#   F_hat
# }

# Semi-parametric CDF with varying threshold and GPD parameters
p_gpd_ecdf_ns <- \(dat, gpd, var) {
  dat_par <- dat |>
    mutate(
      year = as.numeric(substr(date, 1, 4))
    ) |>
    left_join(
      gpd |>
        # select(year, season_month, thresh, scale, shape)
        select(year, month, thresh, scale, shape)
      # by = c("year", "season_month")
    )

  x <- dat_par[[var]]
  u <- dat_par$thresh

  # empirical CDF for observations below their own threshold
  ord <- order(x)
  m <- length(x)
  ecdf_vals <- seq_len(m) / (m + 1)

  ecdf_orig <- numeric(m)
  ecdf_orig[ord] <- ecdf_vals

  # initialise
  F_hat <- ecdf_orig

  # Parametric GPD tail above threshold
  # above <- x > u
  # NAs for u indicate values below threshold
  above <- !is.na(x) & !is.na(u) & x > u


  if (any(above, na.rm = TRUE)) {
    excess <- x[above] - u[above]
    sigma <- dat_par$scale[above]
    xi <- dat_par$shape[above]

    surv_gpd <- ifelse(
      abs(xi) < 1e-8,
      exp(-excess / sigma),
      pmax(1 + xi * excess / sigma, 0)^(-1 / xi)
    )

    # with a q-quantile threshold, exceedance prob is approx 1-q
    # estimate it empirically in case qgam gives slightly diff exceedance rate
    p_exc <- mean(above, na.rm = TRUE)

    F_hat[above] <- 1 - p_exc * surv_gpd
  }

  F_hat
}

# perform Laplace transformation for each variable and season
# debugonce(trans_marg_ns_one)
# debugonce(p_gpd_ecdf_ns)
# laplace_season_var <- lapply(seq_along(gpd_var), \(i) {
#   bind_rows(mclapply(station_names, \(x) { # loop through stations
#     # bind_rows(lapply(station_names, \(x) { # loop through stations
#     # print(i)
#     # print(x)
#     trans_marg_ns_one(
#       data_orig = data2,
#       gpd_fit = gpd_var[[i]],
#       var = season_var_df$var[[i]],
#       season = season_var_df$season[[i]],
#       station = x
#     )
#   })) |>
#     select(name, laplace)
# })
# names(laplace_season_var) <- names(gpd_var)

laplace_season_var <- lapply(seq_len(nrow(season_var_df)), \(i) {
  bind_rows(lapply(station_names, \(station) { # loop through stations
    # bind_rows(lapply(station_names, \(x) { # loop through stations
    # print(i)
    # print(x)
    trans_marg_ns_one(
      data_orig = data2,
      # gpd_fit = gpd_var[[i]],
      gpd_fit = gpd_var[[season_var_df$var[[i]]]],
      var = season_var_df$var[[i]],
      season = season_var_df$season[[i]],
      station = station
    )
  })) |>
    select(name, laplace)
})
# names(laplace_season_var) <- names(gpd_var)
names(laplace_season_var) <- season_var_df$lst_name

# function to join together rain and temperature transformed data for a season
trans_fun <- \(x, season) {
  # join for each variable
  laplace_df <- cbind(
    # laplace_season_var$rain_ |>
    x[[paste0("temp_", season)]] |>
      rename(temp = laplace),
    # laplace_season_var$temp_ |>
    # x[[paste0("rain_", season)]] |>
    x[[paste0(var_dep, "_", season)]] |>
      select(-name) |>
      # rename(rain = laplace)
      rename(!!var_dep := laplace)
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

# split by season to make it easier to use with CE later
marg_season_evgam <- lapply(seasons, \(season) {
  trans_fun(laplace_season_var, season)
})
names(marg_season) <- seasons

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
      # ggplot(aes(x = temp, y = rain, colour = type)) +
      ggplot(aes(x = temp, y = get(var_dep), colour = type)) +
      geom_point(alpha = 0.8) +
      cecl_theme() +
      ggtitle(paste0(x, " - ", season))
  })
}

pdf("plots/02_app/03_laplace_trans_compare_winter.pdf", width = 12, height = 8)
plot_diff("Winter")
dev.off()

pdf("plots/02_app/03_laplace_trans_compare_spring.pdf", width = 12, height = 8)
plot_diff("Spring")
dev.off()

pdf("plots/02_app/03_laplace_trans_compare_summer.pdf", width = 12, height = 8)
plot_diff("Summer")
dev.off()

pdf("plots/02_app/03_laplace_trans_compare_autumn.pdf", width = 12, height = 8)
plot_diff("Autumn")
dev.off()


#### Plot strange countries ####

# Inspect country with outlier
df <- trans_marg_ns_one(
  data2,
  gpd_var$temp,
  "temp",
  "Summer",
  "Vigo Aeropuerto"
)

df |>
  arrange(desc(abs(laplace))) |>
  select(
    date,
    temp,
    F_hat,
    laplace
  ) |>
  head(20)

gpd_var$temp |>
  filter(name == "Vigo Aeropuerto") |>
  summarise(
    min_xi = min(shape),
    max_xi = max(shape)
  )

gpd_var$temp |>
  filter(name == "Vigo Aeropuerto") |>
  summarise(
    min_sigma = min(scale),
    max_sigma = max(scale)
  )

# calculate empirical maximum
n <- sum(
  data2$name == "Vigo Aeropuerto" &
    data2$season == "Summer"
)

# rough scale of largest Laplace value for empirical transformation
log((n + 1) / 2)


#### QQ plots ####

# fitted gpd cdf for excesses
pgpd_excess <- function(excess, scale, shape) {
  ifelse(
    abs(shape) < 1e-8,
    1 - exp(-excess / scale),
    1 - pmax(1 + shape * excess / scale, 0)^(-1 / shape)
  )
}

# build gpd residuals for one station / variable / season
gpd_qq_data_one <- function(data_thresh,
                            gpd_fit,
                            station) {
  d <- data_thresh |>
    filter(name == station) |>
    # select(name, date, year, season_month, thresh, excess)
    select(name, date, year, month, thresh, excess)

  gpd <- gpd_fit |>
    filter(name == station) |>
    # select(year, season_month, scale, shape)
    select(year, month, scale, shape)

  qq_dat <- d |>
    # left_join(gpd, by = c("year", "season_month")) |>
    left_join(gpd, by = c("year", "month")) |>
    filter(
      !is.na(excess),
      !is.na(scale),
      !is.na(shape),
      scale > 0
    ) |>
    mutate(
      # evalulate fitted GPD CDF for each excess
      g_hat = pgpd_excess(excess, scale, shape),
      # clamp probabilities
      g_hat = pmin(pmax(g_hat, 1e-10), 1 - 1e-10),
      # transform to exponential(1) obs
      exp_resid = -log(1 - g_hat)
    ) |>
    arrange(exp_resid) |>
    mutate(
      i = row_number(),
      n = n(),
      # Compute theoretical Exp(1) quantiles at plotting positions
      theoretical = qexp((i - 0.5) / n),
      observed = exp_resid
    ) |>
    # Add a pointwise 95% simulation-free QQ envelope
    mutate(
      lower = qexp(qbeta(0.025, 1:n(), n():1)),
      upper = qexp(qbeta(0.975, 1:n(), n():1))
    )

  # ggplot(qq_dat, aes(theoretical, observed)) +
  #   geom_ribbon(
  #     aes(ymin = lower, ymax = upper),
  #     alpha = 0.2
  #   ) +
  #   geom_abline(intercept = 0, slope = 1, linetype = 2) +
  #   geom_point() +
  #   cecl_theme() +
  #   labs(x = "theoretical", y = "observed")

  qq_dat
}

# qq_df <- bind_rows(lapply(seq_len(nrow(season_var_df)), \(i) {
#   print(season_var_df[var[[i]]])
#   bind_rows(lapply(station_names, \(x) {
#   # bind_rows(mclapply(station_names, \(x) {
#     with(
#       season_var_df,
#       gpd_qq_data_one(
#         # data_thresh = thresh_data_var_qgam[[lst_name[[i]]]],
#         # gpd_fit = gpd_var[[lst_name[[i]]]],
#         data_thresh = thresh_data_var_qgam[[var[[i]]]],
#         gpd_fit = gpd_var[[var[[i]]]],
#         station = x
#       ) |>
#         mutate(season = !!season[[i]], var = !!var[[i]]) |>
#         relocate(season, var, .after = name)
#     )
#   }))
# }))
qq_df <- bind_rows(lapply(vars, \(var) {
  bind_rows(mclapply(station_names, \(station) {
    gpd_qq_data_one(
      data_thresh = thresh_data_var_qgam[[var]],
      gpd_fit = gpd_var[[var]],
      station = station
    ) |>
      mutate(var = !!var) |>
      relocate(var, .after = name)
  }))
}))

qq_df |>
  filter(name == unique(qq_df$name)[2]) |>
  ggplot(aes(theoretical, observed)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_point() +
  facet_wrap(~var) +
  cecl_theme() +
  labs(title = qq_df$name[1], x = "Theoretical", y = "Observed")

qq_plots <- mclapply(station_names, \(station) {
  qq_df |>
    filter(name == station) |>
    ggplot(aes(theoretical, observed)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point() +
    facet_wrap(~var) +
    cecl_theme() +
    labs(title = station, x = "Theoretical", y = "Observed")
})

pdf("plots/02_app/05_qq_plots.pdf", width = 12, height = 8)
qq_plots
dev.off()

# All look great, mercifully!!
# For drought local, looks terrible though !!!

#### Shape plots ####

shape_df <- bind_rows(lapply(seq_along(gpd_var), \(i) {
  gpd_var[[i]] |>
    distinct(name, shape) |>
    mutate(
      # var = season_var_df$var[[i]],
      # season = season_var_df$season[[i]]
      var = names(gpd_var)[[i]]
    )
}))

shape_sf <- shape_df |>
  left_join(
    distinct(data2, name, lon, lat)
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas_ccaa))


shape_sf |>
  # filter(var == "rain", season == "Summer") |>
  # mutate(ind = paste0(var, ", ", season)) |>
  ggplot() +
  geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
  geom_sf(
    aes(fill = shape),
    colour = "black",
    stroke = 1,
    size = 5,
    pch = 21
  ) +
  # scale_size_continuous(
  #   range = c(1, 6),
  #   guide = "legend"
  # ) +
  # facet_wrap(~ind) +
  facet_wrap(~var) +
  cecl_theme(nejm_pal = FALSE, legend.position = "right") +
  # scale_fill_gradient2(
  #   low = "blue3",
  #   mid = "white",
  #   high = "red3",
  #   midpoint = 0,
  #   # limits = c(-rng, rng),
  #   na.value = "grey80"
  # )
  scale_fill_viridis_b()


#### Scale plots ####

scale_df <- bind_rows(lapply(seq_along(gpd_var), \(i) {
  gpd_var[[i]] |>
    # distinct(name, year, season_month, scale) |>
    distinct(name, year, month, season, scale) |>
    mutate(
      # var = season_var_df$var[[i]],
      # season = season_var_df$season[[i]]
      var = names(gpd_var)[[i]]
    )
}))

scale_sf <- scale_df |>
  left_join(
    distinct(data2, name, lon, lat)
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas_ccaa))

# TODO Plot by season here?
scale_sf |>
  # filter(var == "rain", season == "Summer") |>
  filter(year == 2000, var %in% c(vars), season %in% c("Summer", "Winter")) |>
  # filter(year == 2020, var == "temp") |>
  # group_by(name, var, year, season) |>
  group_by(name, var, year) |>
  filter(scale == max(abs(scale))) |>
  ungroup() |>
  group_split(var) |>
  lapply(\(x) {
    x |>
      mutate(ind = paste0(var, ", ", season)) |>
      ggplot() +
      geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
      geom_sf(
        aes(fill = scale),
        colour = "black",
        stroke = 1,
        size = 5,
        pch = 21
      ) +
      facet_wrap(~ind, nrow = 2) +
      cecl_theme(nejm_pal = FALSE, legend.position = "bottom") +
      # scale_fill_gradient2(
      #   low = "blue3",
      #   mid = "white",
      #   high = "red3",
      #   midpoint = 0,
      #   na.value = "grey80"
      # )
      scale_fill_viridis_b() +
      theme(legend.text = element_text(size = 10, angle = 45, hjust = 1))
    # scale_fill_viridis_b()
  }) |>
  wrap_plots()


#### Return level plots ####

# TODO Requires model for P(X > u) (i.e. logistic regression for
# exceedance probability) to calculate return levels properly)

# Calculate GPD return levels from fitted threshold, scale, shape
calc_return_levels <- function(gpd_fit,
                               return_periods = c(10, 20, 50),
                               obs_per_year = 26,
                               q = 0.9) {
  lambda <- obs_per_year * (1 - q)

  tidyr::crossing(
    gpd_fit,
    return_period = return_periods
  ) |>
    mutate(
      return_level = ifelse(
        abs(shape) < 1e-8,
        thresh + scale * log(return_period * lambda),
        thresh + (scale / shape) * ((return_period * lambda)^shape - 1)
      )
    )
}

rl_df <- bind_rows(lapply(seq_along(gpd_var), \(i) {
  calc_return_levels(
    gpd_fit = gpd_var[[i]],
    return_periods = c(10, 20, 50),
    obs_per_year = 26,
    q = 0.9
  ) |>
    mutate(
      # var = season_var_df$var[[i]],
      # season = season_var_df$season[[i]]
      var = names(gpd_var)[[i]]
    )
}))

rl_sf <- rl_df |>
  left_join(
    distinct(data2, name, lon, lat)
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas_ccaa))

rl_sf |>
  filter(
    year == 2018,
    # season == "Summer",
    month == 3,
    return_period == 20
  ) |>
  # mutate(ind = paste0(var, ", ", season)) |>
  ggplot() +
  geom_sf(data = areas_ccaa, fill = NA, colour = "black") +
  geom_sf(
    aes(fill = return_level),
    colour = "black",
    stroke = 1,
    size = 5,
    pch = 21
  ) +
  # facet_wrap(~ind) +
  facet_wrap(~var) +
  cecl_theme(nejm_pal = FALSE, legend.position = "right") +
  scale_fill_viridis_b()

### Fit CE and Cluster ####

# TODO Why does plot_resid have values > 1 on the x-axis??

# TODO Pulling in this data properly???

# TODO Investigate in heatmaps if the "total" heatmap is correct; in some
# plots it just looks the same as the second one (i.e. for temp|rain)
# TODO Change scale on heatmaps to reflect actual dissimilarity values

# TODO Investigate failed CE fits (Getafe Summer is one)
# coef(dep_var_season[[1]]) |> filter(is.na(a))
# - 2 for Summer, 4 for Winter; not that big a deal I guess ?!? Try diff DQU?

# - Threshold stability plots for each location (using bootstrap) (done)
# - Plot CE diagnostic plots for each location (and variable) (done)
#   - Do for all marginal versions (Empirical and GPD) (done)
# TODO - Plot clustering solutions for different thresholds (i.e. sensitivity)
# TODO   - Plot for each conditioning variable also!!
# TODO   - Also do diagnostic plots?

# TODO - Add boostrapping function to package (with methods etc)
# TODO   - Have working with given marginal method/model

# Look at outlying precipitation and temperature, match up with Wikipedia (done)

# TODO Comment on fits + maps + plots!!!

# Find number of clusters (done)
# Plot map of estimated CE parameter values (done)

# Plot heatmap (after clustering) (done)

# Plot clustering maps & heatmaps for both models as well as agg matrix
# - See if perhaps there is a better signal for one model over the other? (done)

# Also cluster for different decades (done)
# - Plot clustering solution for each decade (done)
# - Plot heatmaps for each decade (done)
# - Make dataframe of difference in clustering for both (i) 1960 and
# the decade in question and (ii) the current decade and the previous one (done)
# - Could have map plot with size/shape determined by whether
# point changes! (done)
# - Plot difference in heatmap from 1960 (done)

# TODO Run again with `evgam` and other marginal models/estimates

# TODO Change cecl_clust so that you can specify which var and cond_var to use
# TODO Tidy up plot_scree function, should be OOP like rest of CeCl package

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
library(terra)
devtools::load_all("../CeCl")
library(glue)
library(fs)

source("src/00_functions.R")


#### Metadata ####

# dep_vars <- c("rain", "drought_local", "drought_global")
# dep_vars <- c("rain")
# dep_vars <- c("drought_local_norm")
# dep_vars <- c("drought_local")
dep_vars <- c("drought_local_rev")

# temp_var <- "temp"
temp_var <- "temp_max"

# k_vals <- 1:5
k_vals <- 1:6

decades <- seq(1960, 2010, by = 10)

# seasons <- c("Summer", "Winter", "Year")
# seasons <- c("Summer", "Winter")
seasons <- c("Winter", "Spring", "Summer", "Autumn")

cond_prob <- 0.8 # most stable based on looking at clustering plots
# cond_prob <- 0.85 # default dependence threshold quantile
# cond_prob <- 0.9 # try higher, CE diagnostics aren't good here!
# cond_prob <- 0.95

# TODO Make a file that just sources this one and changes these variables!
# save_dir <- "plots/02_app/"
# save_dir <- "plots/02_app/evgam/"
save_dir <- "plots/02_app/roll_emp/"
# save_dir <- "plots/02_app/mgcv/"
# save_dir <- "plots/02_app/mgcv/gpd/"
# save_dir <- "plots/02_app/mgcv_spat/"
# save_dir <- "plots/02_app/mgcv_spat/gpd/"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# file <- "data/02_app/marg_season_emp.rds"
# file <- "data/02_app/marg_season_evgam.rds"
file <- "data/02_app/marg_season_roll_emp.rds"
# file <- "data/02_app/marg_season_mgcv.rds"
# file <- "data/02_app/marg_season_mgcv_gpd.rds"
# file <- "data/02_app/marg_season_mgcv_spat.rds"
# file <- "data/02_app/marg_season_mgcv_gpd_spat.rds"

# marg_season_decade <- readRDS("data/02_app/marg_season_decade_emp.rds")

marg_season <- readRDS(file)
if (!is.null(names(marg_season))) {
  seasons <- names(marg_season)
}
# marg_season_decade <- readRDS("data/02_app/marg_season_decade_emp.rds")


#### Functions ####

# pull dependence paramaters into df
dep_to_df <- \(dep) {
  # pull a and b values for each location and variable
  ab_vals <- lapply(dep, \(x) {
    list(x[[1]][1:2], x[[2]][1:2])
  })

  vars <- names(dep[[1]])

  # convert to dataframe
  bind_rows(lapply(seq_along(ab_vals), \(i) {
    # browser()
    data.frame(
      "name" = names(ab_vals)[i],
      # "vars" = c("rain", "rain", "wind_speed", "wind_speed"),
      "vars" = c(vars[[1]], vars[[1]], vars[[2]], vars[[2]]),
      "parameter" = c("a", "b", "a", "b"),
      "value" = c(
        ab_vals[[i]][[1]][1], ab_vals[[i]][[1]][2],
        ab_vals[[i]][[2]][1], ab_vals[[i]][[2]][2]
      )
    )
  })) |>
    # TODO: Investigate NAs here
    filter(!is.na(name))
}

# function to convert to ECDF
calc_ecdf <- \(x) {
  apply(x, 2, \(dat_spec) {
    dat_spec_ord <- order(dat_spec)
    dat_spec_sort <- dat_spec[dat_spec_ord]

    # calculate ECDF
    m <- length(dat_spec)
    ecdf_vals <- (seq_len(m)) / (m + 1)
    # convert back to original order
    ecdf_dat_ord <- numeric(m)
    ecdf_dat_ord[dat_spec_ord] <- ecdf_vals
    ecdf_dat_ord
  })
}

# perform PIT on ECDF to get data on Laplace margins
trans_fun <- \(x) {
  # laplace_trans(calc_ecdf(x))
  dlaplace(calc_ecdf(x))
}

# pull dependence parameters and residuals out separately
pull_element <- \(x, element) {
  if (element %in% names(x)) {
    return(x[[element]])
  } else {
    lapply(x, pull_element, element)
  }
}

# function to perform bootstrapping for CE model for ECDF
boot_ce_ecdf <- \(
  dep,
  orig,
  transformed,
  cond_prob = 0.9,
  R = 100,
  trace = 10,
  ncores = 1,
  fixed_b = FALSE,
  ...
) {
  # pull data
  # TODO Should pull fixed_b from this!
  arg_vals <- list("cond_prob" = cond_prob) # TODO Add any others here
  dependence <- dep # dependence parameters

  # calculate marginal threshold
  # TODO Will have to change for real data, as can't use qgpd
  # marg_val <- do.call(evd::qgpd, c(list(marg_prob), marg_pars))

  # extract marginal and dependence thresholds
  # TODO Implement for multiple locations
  thresh_dep <- lapply(dependence, \(x) {
    res <- vapply(x, \(y) {
      y[rownames(y) == "dth", , drop = FALSE]
    }, numeric(length(x) - 1))
    if (!is.matrix(res)) {
      res <- matrix(res, nrow = 1)
      colnames(res) <- names(x)
    }
    res[1, , drop = FALSE]
  })

  # Parallel setup
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # Pull start values
  start <- lapply(dependence, \(x) {
    lapply(x, \(y) {
      # scale back towards zero in case point est on edge of original parameter
      # space and falls off edge of constrained space for bootstrap sample
      if (fixed_b) {
        y[c("a", "b"), , drop = FALSE] * c(0.75, 1)
      } else {
        y[c("a", "b"), , drop = FALSE] * 0.75
      }
    })
  })

  # Function to prepare bootstrapped data for each location and conditioned var
  # prep_boot_loc_dat <- \(
  #   dep_loc, trans_loc, orig_loc, thresh_dep, n_pass = 3
  # ) {
  #   # pull bootstrap sample
  #   # (must be different for each loc as nrows may differ)
  #   indices <- sample(seq_len(nrow(trans_loc)), replace = TRUE)
  #   trans_loc_boot <- trans_loc[indices, ]
  #
  #   # Reorder bootstrap sample to have the same order as original data
  #   test <- FALSE
  #   # which <- which(names(marg_loc) %in% cond_var)
  #   while (test == FALSE) {
  #     for (j in seq_along(dep_loc)) {
  #       # replace ordered Yi with ordered sample from standard Laplace CDF
  #       u <- matrix(runif(nrow(trans_loc_boot)))
  #       trans_loc_boot[
  #         order(trans_loc_boot[, j]), j
  #         # ] <- sort(laplace_trans(u)[, 1])
  #       ] <- sort(dlaplace(u)[, 1])
  #     }
  #     # need exceedance in cond. var, and also in other vars for these rows
  #     if (any(colSums(sweep(trans_loc_boot, 2, thresh_dep, FUN = ">")) > 0)) {
  #       test <- TRUE
  #     }
  #   }
  #
  #   # convert variables to original scale
  #   # orig_loc_boot <- inv_semi_par_cdf(
  #   #   # inverse Laplace transform to CDF
  #   #   inv_laplace_trans(trans_loc_boot),
  #   #   # original data for where semiparametric thresholding occurs
  #   #   orig_loc,
  #   #   # TODO Have to change if parameters are different
  #   #   lapply(seq_len(ncol(orig_loc)), \(i) {
  #   #     setNames(marg_pars[c("scale", "shape", "loc")], c("sigma", "xi", "thresh"))
  #   #   })
  #   # )
  #
  #   # convert from Laplace scale back to ECDF (no need to get original data)
  #   # ecdf_loc_boot <- inv_laplace_trans(trans_loc_boot)
  #   ecdf_loc_boot <- plaplace(trans_loc_boot)
  #
  #   # test for no marg exceedances over sampled points, if so resample w/ nPass
  #   # max_vals <- apply(orig_loc_boot, 2, max, na.rm = TRUE)
  #   # marg_thresh <- marg_val # TODO Will have to change if vars don't have same margins
  #   # if (!all(max_vals > marg_thresh)) {
  #   #   return(list(NA))
  #   # }
  #   return(ecdf_loc_boot)
  # }

  prep_boot_loc_dat <- \(trans_loc, thresh_dep, n_pass = 3) {
    vars <- colnames(trans_loc)

    for (pass in seq_len(n_pass)) {
      # nonparametric bootstrap sample
      idx <- sample(seq_len(nrow(trans_loc)), replace = TRUE)
      trans_boot_rank_ref <- trans_loc[idx, , drop = FALSE]

      # independent standard Laplace sample
      z_boot <- matrix(
        rlaplace(length(trans_boot_rank_ref)),
        nrow = nrow(trans_boot_rank_ref),
        ncol = ncol(trans_boot_rank_ref),
        dimnames = dimnames(trans_boot_rank_ref)
      )

      # reorder each simulated margin to match bootstrap ranks
      trans_boot <- z_boot
      for (v in vars) {
        trans_boot[order(trans_boot_rank_ref[, v]), v] <-
          sort(z_boot[, v])
      }

      # require exceedances for each conditioning variable
      has_exceed <- colSums(
        sweep(trans_boot, 2, thresh_dep[1, vars], FUN = ">")
      ) > 0

      if (all(has_exceed)) {
        return(trans_boot)
      }
    }

    NA
  }

  # perform bootstrapping
  boot_fits <- loop_fun(seq_len(R), \(i) {
    # browser()
    if (i %% trace == 0) {
      system(sprintf("echo %s", paste(i, "replicates done")))
    }
    # dat_ecdf_boot <- lapply(seq_along(dependence), \(j) { # loop through locations
    #   # prepare bootstrapped data for location j
    #   # TODO Change argument order?
    #   dat_spec <- prep_boot_loc_dat(
    #     dependence[[j]],
    #     transformed[[j]],
    #     orig[[j]],
    #     thresh_dep[[j]]
    #   )
    #   # check if NA, if so then rerun
    #   if (all(is.na(dat_spec)) && n_pass > 1) {
    #     for (i in seq_len(n_pass - 1)) {
    #       dat_spec <- prep_boot_loc_dat(
    #         marginal[[j]], transformed[[j]], dependence[[j]]
    #       )
    #       if (!all(is.na(dat_spec))) {
    #         break
    #       }
    #     }
    #   } else if (all(is.na(dat_spec))) {
    #     stop(paste0(
    #       "Failed to generate bootstrapped data after ", n_pass, " attempts"
    #     ))
    #   }
    #   return(dat_spec)
    # })
    #
    # # transform to Laplace margins
    # # dat_boot_trans <- lapply(dat_boot, \(x) {
    # #   laplace_trans(do.call(evd::pgpd, c(list(x), marg_pars)))
    # # })
    # # dat_boot_trans <- lapply(dat_ecdf_boot, laplace_trans)
    # dat_boot_trans <- lapply(dat_ecdf_boot, dlaplace)

    dat_boot_trans <- lapply(seq_along(dependence), \(j) {
      prep_boot_loc_dat(
        trans_loc = transformed[[j]],
        thresh_dep = thresh_dep[[j]]
      )
    })

    # refit CE model using same dependence quantile
    # TODO Could change to mapply?
    vars <- colnames(dat_boot_trans[[1]])
    fit_boot <- lapply(seq_along(dat_boot_trans), \(j) {
      # o <- CeCl:::ce_optim(
      #   Y = dat_boot_trans[[j]],
      #   dqu = cond_prob,
      #   start = start[[j]],
      #   control = list(maxit = 1e6),
      #   fixed_b = fixed_b,
      #   ...
      # )
      # loop over conditioning variables and use specific start values
      o_spec <- lapply(vars, \(cond_var) {
        # start_spec <- setNames(as.vector(start[[j]][[cond_var]]), c("a", "b"))
        # TODO Check if this is the right starting value
        start_spec <- start[[j]][[cond_var]][c("a", "b"), ]
        names(start_spec) <- c("a", "b")

        CeCl:::ce_optim(
          Y = dat_boot_trans[[j]],
          dqu = cond_prob,
          # start = start[[j]],
          start = start_spec,
          cond_vars = cond_var,
          control = list(maxit = 1e6),
          fixed_b = fixed_b,
          nruns = 2,
          ...
        )
      })

      # join together (as if ran in `ce_optim` across all variables)
      Reduce(`c`, o_spec)
    })


    # fit_boot <- mapply(\(dat, start_val) {
    #   o <- ce_optim(
    #     # Y = dat_boot_trans[[j]],
    #     Y = dat,
    #     dqu = cond_prob,
    #     # start = start[[j]],
    #     start = start_val,
    #     control = list(maxit = 1e6),
    #     ...
    #   )
    # }, dat_boot_trans, start)

    # extract just parameters
    fit_boot_pars <- lapply(fit_boot, \(x) lapply(x, `[[`, "params"))
    names(fit_boot_pars) <- names(dep)

    # Calculate conditional expectation at marg_val
    # return(lapply(fit_boot_pars, \(x) {
    #   lapply(x, \(y) {
    #     y["a", ] * marg_val + (marg_val^(y["b", ])) * y["m", ]
    #   })
    # }))
    # dep_out <- lapply(boot_fits, `[[`, "dependence") |>
    dep_out <- purrr::map_dfr(names(fit_boot_pars), function(loc_name) {
      # loc <- boot_sample[[1]]
      loc <- fit_boot_pars[[loc_name]]

      # loop over variables
      purrr::map_dfr(names(loc), function(cond_var) {
        # Get the matrix for the current conditioning variable
        var_mat <- loc[[cond_var]]

        # be careful if var_mat is not a matrix
        if (is.matrix(var_mat)) {
          var_names <- colnames(var_mat)
        } else {
          var_names <- names(loc)[!names(loc) == cond_var]
          var_mat <- as.matrix(var_mat)
        }

        # Create a tidy data frame for this location's matrix
        tibble(
          parameter = rep(c("a", "b"), each = ncol(var_mat)),
          # vars      = rep(colnames(var_mat), times = 2),
          vars      = rep(var_names, times = 2),
          value     = c(var_mat["a", ], var_mat["b", ]),
          cond_var  = rep(cond_var, ncol(var_mat) * 2),
          name      = rep(loc_name, ncol(var_mat) * 2)
        )
      })
    })

    # Combine all the lists of data frames into one final data frame
    dep_out <- bind_rows(dep_out)
  })

  return(boot_fits)
}

# function to plot bootstrapped parameter estimates for different quantiles
plot_boot_quant <- \(
  data_lst,
  quantiles = seq(0.5, 0.9, by = 0.1),
  constrain = TRUE,
  R = 10,
  ncores = 1,
  county_key_df = NULL,
  start = c("a" = 0.01, "b" = 0.01),
  ...
) {
  # transform data
  Y_lst <- lapply(data_lst, trans_fun)

  # fit model at each threshold
  ce_fit_quant <- mclapply(quantiles, \(q) {
    ret <- lapply(Y_lst, \(y) {
      o <- CeCl:::ce_optim(
        Y = y,
        dqu = q,
        constrain = constrain,
        start = start,
        # nruns = 3,
        nruns = 2,
        ...
      )
    })
    locs <- names(ret)
    ret <- lapply(c("resid", "params"), \(x) {
      setNames(pull_element(ret, x), locs)
    })
    names(ret) <- c("residual", "dependence")
    return(ret)
  }, mc.cores = ncores)

  # extract parameter estimates and label by quantile (for plotting)
  ab_df <- bind_rows(lapply(seq_along(ce_fit_quant), \(i) {
    dep_to_df(ce_fit_quant[[i]]$dependence) |>
      mutate(quantile = quantiles[i])
  }))

  # if (!is.null(county_key_df)) {
  #   ab_df <- ab_df |>
  #     left_join(county_key_df, by = c("name"))
  # }

  # now for each model fit, bootstrap
  # TODO Again, need to make these parameters
  boot_fit_quant <- mclapply(seq_along(ce_fit_quant), \(i) {
    boot_ce_ecdf(
      dep = ce_fit_quant[[i]]$dependence,
      orig = data_lst,
      transformed = Y_lst,
      cond_prob = quantiles[i],
      R = R,
      trace = R + 1,
      constrain = constrain,
      ...
    )
  }, mc.cores = ncores)

  # join all together
  boot_ab_df <- bind_rows(lapply(seq_along(boot_fit_quant), \(i) {
    ret <- bind_rows(boot_fit_quant[[i]], .id = "run") |>
      mutate(quantile = quantiles[i])
  }))

  loc_df <- ab_df |>
    distinct(across(matches("name") | matches("county")))
  plots <- lapply(seq_along(locs), \(i) {
    # browser()
    # title <- ifelse(
    #   is.null(county_key_df),
    #   loc_df$name[[i]],
    #   paste0(loc_df$name[i], " - ", loc_df$county[i])
    # )
    title <- loc_df$name[[i]]

    ab_df |>
      filter(name == loc_df$name[[i]]) |>
      mutate(vars = paste(parameter, vars, sep = " - ")) |>
      ggplot(aes(x = quantile, y = value)) +
      geom_point(
        data = boot_ab_df |>
          filter(name == loc_df$name[[i]]) |>
          mutate(vars = ifelse(vars == "windspeed", "wind_speed", vars)) |>
          mutate(vars = paste(parameter, vars, sep = " - ")),
        # mutate(
        #   vars = ifelse(
        #     vars == "windspeed",
        #     "Precipitation | Wind Speed",
        #     "Wind Speed | Precipitation"
        #   )
        # ),
        colour = "black",
        size = 2,
        alpha = 0.7
      ) +
      geom_line() +
      geom_point(size = 5, colour = "orange", alpha = 0.7) +
      # facet_wrap(~ parameter + vars, scales = "free") +
      facet_wrap(
        ~vars,
        ncol = 2,
        labeller = as_labeller(c(
          # "a - rain"       = "alpha ~ ' - ' ~ 'Precipitation | Temperature'",
          "a - rain" = "alpha ~ ' - ' ~ 'Precipitation | Temperature '",
          "a - temp" = "alpha ~ ' - ' ~ 'Temperature | Precipitation'",
          "b - rain" = "beta ~ ' - ' ~ 'Precipitation | Temperature'",
          "b - temp" = "beta ~ ' - ' ~ 'Temperature | Precipitation'"
        ), default = label_parsed)
      ) +
      # theme +
      cecl_theme() +
      scale_x_continuous(breaks = quantiles) +
      labs(
        x     = "q",
        y     = "Parameter estimate",
        title = title
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  names(plots) <- loc_df$name

  return(plots)
}


#### Load Data ####

# data
data <- readr::read_csv(
  "data/02_app/ecad_clean.csv.gz"
) |>
  mutate(decade = factor(floor(year(date) / 10) * 10, levels = decades))

if (dep_vars == "rain") {
  data <- data |>
    filter(rain > 0)
}

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
if (dep_vars == "drought_local_rev") {
  data <- data |>
    mutate(drought_local = -drought_local)
  dep_vars <- c("drought_local")
}

# marginal fits
marg_season <- readRDS(file)
if (!is.null(names(marg_season))) {
  seasons <- names(marg_season)
}
# marg_season_decade <- readRDS("data/02_app/marg_season_decade_emp.rds")

# # Remove large outliers !
# marg_season$Winter$transformed$Daroca <- marg_season$Winter$transformed$Daroca[
#   marg_season$Winter$transformed$Daroca[, 1] < 25,
# ]
# marg_season$Winter$transformed$Zamora <- marg_season$Winter$transformed$Zamora[
#   marg_season$Winter$transformed$Zamora[, 1] < 25,
# ]

# load map of continental Spain to use as background in plots
# TODO Save map so that we can pull with no internet!!!
# areas <- mapSpain::esp_get_munic_siane(moveCAN = TRUE) |>
areas <- read_sf("data/02_app/spain_shapefile.geojson") |>
  filter(
    !ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla")
  )

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# extract point location of each station for plotting on map
pts <- data %>%
  distinct(name = station_name, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas_ccaa))

station_names <- unique(data$station_name)


#### Threshold stability ####

# split data by location, convert to list of matrices
# data_lst <- data |>
#   select(name = station_name, temp, rain) |>
#   mutate(name = as.factor(name)) |>
#   group_split(name, .keep = TRUE)
# locs <- purrr::map_chr(data_lst, ~ as.character(.x$name[1]))
# names(data_lst) <- locs
# data_lst <- lapply(data_lst, \(x) {
#   matrix(
#     c(x$temp, x$rain),
#     ncol = 2,
#     dimnames = list(NULL, c("temp", "rain"))
#   )
# })

data_lst_season <- data |>
  mutate(season = factor(season, levels = seasons)) |>
  group_split(season, .keep = FALSE) |>
  lapply(\(x) {
    y <- x |>
      select(name = station_name, temp, rain) |>
      mutate(name = as.factor(name)) |>
      group_split(name, .keep = TRUE)
    locs <- purrr::map_chr(y, ~ as.character(.x$name[1]))
    names(y) <- locs
    y <- lapply(y, \(x) {
      matrix(
        c(x$temp, x$rain),
        ncol = 2,
        dimnames = list(NULL, c("temp", "rain"))
      )
    })
    y
  })
names(data_lst_season) <- unique(data$season)

# plot params at various quantiles for each loc, w/ bootstrapped samples
set.seed(123)
# quant_plots <- plot_boot_quant(
#   data_lst,
#   quantiles = c(0.6, 0.7, 0.8, 0.85, 0.88, 0.9, 0.95),
#   constrain = FALSE,
#   R = 100,
#   # R = 5,
#   # ncores = 1
#   ncores = parallel::detectCores() - 1,
#   # county_key_df = county_key_df
# )
# names(quant_plots) <- lapply(quant_plots, \(x) {
#   x$data$name[[1]]
# })
# quant_plots <- lapply(data_lst_season, \(data_lst) {
#   p <- plot_boot_quant(
#     data_lst,
#     quantiles = c(0.6, 0.7, 0.8, 0.85, 0.88, 0.9, 0.95),
#     constrain = FALSE,
#     R = 100,
#     ncores = parallel::detectCores() - 1,
#   )
#   names(p) <- lapply(p, \(x) {
#     x$data$name[[1]]
#   })
#   p
# })
#
# saveRDS(quant_plots, "data/02_app/boot_quantiles.RDS")
#
# quant_plots <- loadRDS("data/02_app/boot_quantiles.RDS")
#
# quant_plots$Summer$`Salamanca Aeropuerto`
#
# pdf(paste0(save_dir, "boot_quantiles_winter.pdf"), width = 14, height = 8)
# quant_plots$Winter
# dev.off()
# pdf(paste0(save_dir, "boot_quantiles_summer.pdf"), width = 14, height = 8)
# quant_plots$Winter
# dev.off()


# Conclusion|: For R = 30, 0.85 seems to give quite stable estimates, with
# bootstrap uncertainty increasing thereafter

#### Fit CE ####

# Fit CE model for each variable and season
dep_var_season <- lapply(dep_vars, \(spec_var) {
  lapply(marg_season, \(x) {
    cecl_dep(
      x,
      vars = c("temp", spec_var),
      cond_vars = c("temp", spec_var),
      cond_prob = cond_prob,
      nruns = 3,
      ncores = getOption("mc_cores", 2L)
    )
  })
})
# dep_var_season <- dep_var_season[[1]] # only looking at rain model
names(dep_var_season) <- dep_vars

# remove failed fits
dep_var_season[[1]] <- lapply(seq_along(dep_var_season[[1]]), \(i) {
  # object for editing
  dep_var_season_cp <- dep_var_season[[1]][[i]]
  # pull coefficients to check for failed fits
  coefs_var_season <- coef(dep_var_season_cp)

  rm_locs <- coefs_var_season |>
    filter(is.na(a)) |>
    pull(name)
  # remove if necessary
  if (length(rm_locs) > 0) {
    dep_var_season_cp$dependence[rm_locs] <- NULL
    dep_var_season_cp$transformed[rm_locs] <- NULL
    dep_var_season_cp$residual[rm_locs] <- NULL
  }

  dep_var_season_cp
})
names(dep_var_season[[1]]) <- seasons

# Function to plot CE diagnostics for each location and season
ce_plot <- \(obj, loc, season, dep_var) {
  p <- (
    plot_resid(obj[[season]], type = "ggplot", loc = loc, var = dep_var, cond_var = "temp") +
      plot_resid(obj[[season]], type = "ggplot", loc = loc, var = "temp", cond_var = dep_var)
  ) /
    (
      plot_quantile(obj[[season]], type = "ggplot", loc = loc, var = "temp", cond_var = dep_var) +
        plot_quantile(obj[[season]], type = "ggplot", loc = loc, var = "temp", cond_var = dep_var)
    ) +
    plot_layout() +
    plot_annotation(
      title = paste0(loc, ", ", season) # ,
      # theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
    )

  # p <- suppressWarnings(p &
  #   plot_layout() +
  #   ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5)))
  p <- suppressWarnings(p &
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5)))
  suppressWarnings(p)
}

dep_var_season_orig <- dep_var_season

# test out
# debugonce(plot_resid)
# debugonce(ce_plot)
ce_plot(dep_var_season_orig[[1]], "A Coruna", "Summer", dep_var = dep_vars[[1]])
ce_plot(dep_var_season_orig[[1]], "Valencia", "Winter", dep_var = dep_vars[[1]])

ce_plots <- lapply(seasons, \(season) {
  lapply(station_names, \(name) {
    tryCatch(
      {
        ce_plot(dep_var_season[[1]], name, season, dep_var = dep_vars[[1]])
      },
      error = \(e) {
        return(NA) # TODO Investigate failed fits!
      }
    )
  })
})
names(ce_plots) <- seasons

# save
# TODO Change to loop, really ugly
pdf(paste0(save_dir, "ce_diagnostics_winter.pdf"), width = 10, height = 8)
ce_plots$Winter
dev.off()
pdf(paste0(save_dir, "ce_diagnostics_spring.pdf"), width = 10, height = 8)
ce_plots$Spring
dev.off()
pdf(paste0(save_dir, "ce_diagnostics_summer.pdf"), width = 10, height = 8)
ce_plots$Summer
dev.off()
pdf(paste0(save_dir, "ce_diagnostics_autumn.pdf"), width = 10, height = 8)
ce_plots$Autumn
dev.off()

# flatten
dep_var_season_orig <- dep_var_season
dep_var_season <- purrr::list_flatten(
  dep_var_season,
  name_spec = "{outer}_{inner}"
)

# dataframe of list names for later
name_df <- tibble(name = names(dep_var_season)) |>
  tidyr::separate_wider_regex(
    name,
    c(
      var = ".*",
      "_",
      season = "[^_]+"
    )
  )
# name_df <- tibble(season = c("Summer", "Winter"))

# calculate alpha, beta estimates
ab_var_season <- lapply(dep_var_season, coef)

# plot results on map
ab_plots_season <- lapply(seq_along(dep_var_season), \(i) {
  map_plots <- plot_ab_map(
    ab_var_season[[i]],
    distinct(data, name = station_name, lon, lat),
    areas_ccaa,
    n_breaks = 8
  )
  map_plots[[1]] / map_plots[[2]] +
    plot_annotation(
      title = name_df$season[i],
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
})

ab_plots_season[[1]]
ab_plots_season[[2]]
# ab_plots_season[[3]]

# pdf("plots/02_app/evgam/04_ab_map.pdf", width = 12, height = 8)
pdf(paste0(save_dir, "04_ab_map.pdf"), width = 12, height = 8)
ab_plots_season[[1]]
ab_plots_season[[2]]
ab_plots_season[[3]]
ab_plots_season[[4]]
dev.off()


#### Cluster ####

if (is.null(names(marg_season))) {
  names(marg_season) <- seasons
}

# calculate distances
set.seed(123)
dist_var_season <- lapply(seq_along(dep_var_season), \(i) {
  cecl_dist(
    dep_var_season[[i]],
    marg_season[[name_df$season[[i]]]],
    ncores = getOption("mc_cores", 2L),
  )
})
names(dist_var_season) <- names(dep_var_season)

# find k via Elbow plot
twgss_vals <- bind_rows(lapply(seq_along(dist_var_season), \(i) {
  data.frame(
    "twgss"    = plot_scree(dist_var_season[[i]]$dist_mat, k_vals, "none"),
    "k"        = k_vals,
    "cond_var" = name_df$var[[i]],
    "season"   = name_df$season[[i]]
  )
}))

twgss_vals |>
  mutate(fac = paste(cond_var, " - ", season)) |>
  ggplot(aes(x = k, y = twgss)) +
  geom_point() +
  geom_line() +
  facet_wrap(~fac, scales = "free_y") +
  scale_x_continuous(breaks = c(k_vals)) +
  cecl_theme()

# For Summer, clear elbow at k = 3
# For Winter, less clear, maybe k = 4 though?
# k_vec <- ifelse(grepl("Year", names(dist_var_season)), 3, 2)

# try for k = 2, 3, 4
k_vec <- c(2, 3, 4)

# # Cluster for each variable and season
# clust_var_season <- lapply(seq_along(dep_var_season), \(i) {
#   cecl_clust(
#     dist_var_season[[i]],
#     cond_var = name_df$var[i],
#     k        = k_vec[i],
#     ncores   = getOption("mc_cores", 2L)
#   )
# })

# Cluster for each variable and season
clust_var_season_k <- lapply(k_vec, \(k) {
  print(paste0("k = ", k))
  ret <- lapply(seq_along(dep_var_season), \(i) {
    print(paste0("i = ", i))
    cecl_clust(
      dist_var_season[[i]],
      cond_var = name_df$var[i],
      k        = k,
      ncores   = getOption("mc_cores", 2L)
    )
  })
  names(ret) <- names(dep_var_season)
  ret
})
names(clust_var_season_k) <- k_vec

map_plots <- lapply(seq_along(k_vec), \(i) {
  mclapply(seq_along(dep_var_season), \(j) {
    plt_clust_map(
      pts,
      areas_ccaa,
      clust_var_season_k[[i]][[j]]$pam,
      plot_medoids = FALSE,
      pt_size = 6
    ) +
      ggtitle(paste0(name_df$var[j], " - ", name_df$season[j]))
  })
})

wrap_plots(map_plots[[1]]) # plots for k = 2
wrap_plots(map_plots[[2]])
wrap_plots(map_plots[[3]])

# Also plot heatmaps
heatmap_plots <- lapply(seq_along(map_plots), \(i) {
  lapply(seq_along(dep_var_season), \(j) {
    ggplot(clust_var_season_k[[i]][[j]], which = "image") +
      theme(
        legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(paste0(name_df$var[j], " - ", name_df$season[j]))
  })
})

wrap_plots(heatmap_plots[[1]])
wrap_plots(heatmap_plots[[2]])
wrap_plots(heatmap_plots[[3]])

# Seems to be strongest signal in heatmaps for Winter! And for three
# clusters it looks most sensible (along with elbow plot), but also fine for
# k = 2 and k = 4 (good separation in heatmap)

# save for k = 3
# ggsave("plots/02_app/evgam/04_clust_map.png", wrap_plots(map_plots[[2]]), width = 12, height = 8)
ggsave(paste0(save_dir, "04_clust_map.png"), wrap_plots(map_plots[[2]]), width = 12, height = 8)

# ggsave("plots/02_app/evgam/04_clust_heatmap.png", wrap_plots(heatmap_plots[[2]]), width = 16, height = 7)
ggsave(paste0(save_dir, "04_clust_heatmap.png"), wrap_plots(heatmap_plots[[2]]), width = 16, height = 7)


#### Cluster for individual dissimilarity matrices and plot ####

# find k via Elbow plot
twgss_vals_var <- bind_rows(lapply(seq_along(dist_var_season), \(i) {
  bind_rows(lapply(names(dist_var_season[[1]]$dist_mats), \(var) {
    data.frame(
      # "twgss"    = plot_scree(dist_var_season[[i]]$dist_mat, k_vals, "none"),
      "twgss"    = plot_scree(dist_var_season[[i]]$dist_mats[[var]], k_vals, "none"),
      "k"        = k_vals,
      # "cond_var" = name_df$var[[i]],
      "cond_var" = var,
      "season"   = name_df$season[[i]]
    )
  }))
}))

twgss_vals_var |>
  mutate(fac = paste(season, " - ", cond_var)) |>
  ggplot(aes(x = k, y = twgss)) +
  geom_point() +
  geom_line() +
  facet_wrap(~fac, scales = "free_y") +
  scale_x_continuous(breaks = c(k_vals)) +
  cecl_theme()
# all confidently say 3 clusters


# Cluster for each variable and season
# TODO Need to make this easier in package (would be very easy)
clust_var <- \(var = "rain") {
  lapply(k_vec, \(k) {
    print(paste0("k = ", k))
    ret <- mclapply(seq_along(dep_var_season), \(i) {
      dist_obj <- dist_var_season[[i]]
      dist_obj$dist_mat <- dist_obj$dist_mats[[var]]

      print(paste0("i = ", i))
      cecl_clust(
        # dist_var_season[[i]],
        dist_obj,
        k        = k,
        ncores   = getOption("mc_cores", 2L)
      )
    })
    names(ret) <- names(dep_var_season)
    ret
  })
}

clust_var_season_k_rain <- clust_var(dep_vars[[1]])
clust_var_season_k_temp <- clust_var("temp")

names(clust_var_season_k_rain) <- names(clust_var_season_k_temp) <- k_vec

clust_sets <- list(
  "Total" = clust_var_season_k,
  # "rain"  = clust_var_season_k_rain,
  "temp | rain" = clust_var_season_k_rain,
  "rain | temp" = clust_var_season_k_temp
)

# out_dir <- "plots/02_app/evgam/001_cluster_solutions_by_distance/"
out_dir <- paste0(save_dir, "001_cluster_solutions_by_distance/")
dir_create(out_dir)

plot_one_map <- \(clust_sets, dist_type, k, season, var = "rain", title = NULL) {
  clust_obj_spec <- clust_sets[[dist_type]][[as.character(k)]][[
    paste(var, season, sep = "_")
  ]]$pam

  p <- plt_clust_map(
    pts,
    areas_ccaa,
    # clust_sets[[dist_type]][[obj_name]]$pam,
    clust_obj_spec,
    plot_medoids = FALSE,
    pt_size = 6
  ) +
    cecl_theme()

  if (is.null(title)) {
    p <- p + ggtitle(dist_type)
  } else {
    p <- p + ggtitle(title)
  }

  p
}

plot_one_heatmap <- \(clust_sets, dist_type, k, season, var = "rain", title = NULL) {
  # ggplot(clust_sets[[dist_type]][[obj_name]], which = "image") +
  p <- ggplot(
    clust_sets[[dist_type]][[as.character(k)]][[paste(var, season, sep = "_")]],
    # fill_limits = c(0, 0.15),
    # fill_limits = c(0, 0.05),
    # fill_breaks = seq(0, 0.15, by = 0.03),
    # fill_breaks = seq(0, 0.05, length.out = 6),
    which = "image"
  ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, angle = 45, hjust = 1)
    )

  if (is.null(title)) {
    p <- p +
      ggtitle(dist_type)
  } else {
    p <- p +
      ggtitle(title)
  }
  p
}

make_season_page <- \(season, k, type = c("map", "heatmap"), title = NULL, var = "rain") {
  type <- match.arg(type)

  plots <- lapply(names(clust_sets), \(dist_type) {
    if (type == "map") {
      plot_one_map(clust_sets, dist_type, k, season, title, var = var)
    } else {
      plot_one_heatmap(clust_sets, dist_type, k, season, title, var = var)
    }
  })

  wrap_plots(plots, nrow = 1) +
    plot_annotation(title = paste0(type, ": k = ", k, ", ", season))
}

# debugonce(plot_one_map)
# debugonce(plot_one_heatmap)
for (k in k_vec) {
  pdf(glue("{out_dir}/cluster_maps_k{k}.pdf"), width = 14, height = 5)

  for (season in seasons) {
    print(k)
    print(season)
    print(make_season_page(season, k, type = "map", var = dep_vars[[1]]))
  }

  dev.off()

  pdf(glue("{out_dir}/cluster_heatmaps_k{k}.pdf"), width = 18, height = 8)

  for (season in seasons) {
    print(make_season_page(season, k, type = "heatmap", var = dep_vars[[1]]))
  }

  dev.off()
}


#### Cluster for each variable for different dependence quantiles ####

clust_plot_dqu <- \(dist_mat = "Total", dqu = 0.9) {
  # Fit CE model for each variable and season
  dep_var_season <- lapply(dep_vars, \(spec_var) {
    lapply(marg_season, \(x) {
      cecl_dep(
        x,
        vars = c("temp", spec_var),
        cond_vars = c("temp", spec_var),
        # cond_prob = cond_prob,
        cond_prob = dqu,
        nruns = 3,
        ncores = getOption("mc_cores", 2L)
      )
    })
  })
  names(dep_var_season) <- dep_vars
  # remove failed fits
  dep_var_season[[1]] <- lapply(seq_along(dep_var_season[[1]]), \(i) {
    # object for editing
    dep_var_season_cp <- dep_var_season[[1]][[i]]
    # pull coefficients to check for failed fits
    coefs_var_season <- coef(dep_var_season_cp)

    rm_locs <- coefs_var_season |>
      filter(is.na(a)) |>
      pull(name)
    # remove if necessary
    if (length(rm_locs) > 0) {
      dep_var_season_cp$dependence[rm_locs] <- NULL
      dep_var_season_cp$transformed[rm_locs] <- NULL
      dep_var_season_cp$residual[rm_locs] <- NULL
    }

    dep_var_season_cp
  })
  names(dep_var_season[[1]]) <- seasons

  # Produce CE diagnostic plots
  ce_plots <- lapply(seasons, \(season) {
    lapply(station_names, \(name) {
      tryCatch(
        {
          ce_plot(dep_var_season[[1]], name, season)
        },
        error = \(e) {
          return(NA) # TODO Investigate failed fits!
        }
      )
    })
  })
  names(ce_plots) <- seasons

  # flatten
  dep_var_season <- purrr::list_flatten(
    dep_var_season,
    name_spec = "{outer}_{inner}"
  )

  # dataframe of list names for later
  name_df <- tibble(name = names(dep_var_season)) |>
    tidyr::separate_wider_regex(
      name,
      c(
        var = ".*",
        "_",
        season = "[^_]+"
      )
    )
  # name_df <- tibble(season = c("Summer", "Winter"))


  # calculate distances
  set.seed(123)
  dist_var_season <- lapply(seq_along(dep_var_season), \(i) {
    cecl_dist(
      dep_var_season[[i]],
      marg_season[[name_df$season[[i]]]],
      ncores = getOption("mc_cores", 2L),
    )
  })

  # replace dist_mat with desired one, if specified
  if (dist_mat != "Total") {
    dist_var_season <- lapply(dist_var_season, \(x) {
      y <- x
      y$dist_mat <- y$dist_mats[[dist_mat]]
      y
    })
  }
  names(dist_var_season) <- names(dep_var_season)

  # Cluster each variable and season
  k <- 3
  clust_var_season <- mclapply(seq_along(dep_var_season), \(i) {
    cecl_clust(
      dist_var_season[[i]],
      cond_var = name_df$var[i],
      k = k,
      ncores = getOption("_cores", 2L)
    )
  })
  names(clust_var_season) <- names(dep_var_season)

  # Produce cluster maps and heatmaps
  map_plots <- mclapply(seq_along(dep_var_season), \(j) {
    plt_clust_map(
      pts,
      areas_ccaa,
      clust_var_season[[j]]$pam,
      plot_medoids = FALSE,
      pt_size = 6
    ) +
      ggtitle(paste0(name_df$var[j], " - ", name_df$season[j]))
  })

  map_plots <- wrap_plots(map_plots)

  # Also plot heatmaps
  heatmap_plots <- mclapply(seq_along(dep_var_season), \(j) {
    ggplot(clust_var_season[[j]], which = "image") +
      theme(
        legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(paste0(name_df$var[j], " - ", name_df$season[j]))
  })
  heatmap_plots <- wrap_plots(heatmap_plots)

  return(list(
    ce_plots = ce_plots,
    map_plots = map_plots,
    heatmap_plots = heatmap_plots
  ))
}

# remove large objects from previous runs
rm(ab_plots_season, ce_plots, clust_sets, dep_var_season, dist_var_season, clust_var_season_k, clust_var_season_k_rain, clust_var_season_k_temp)
gc()

# dqu_vals <- c(0.8, 0.85, 0.88, 0.9, 0.92, 0.95)
# TODO Very memory intensive, find way to make it not so bad!
# Code may just be badly written??
# dqu_vals <- c(0.8, 0.85, 0.88, 0.9)
dqu_vals <- c(0.75, 0.8, 0.85, 0.88, 0.9)
# dqu_vals <- dqu_vals[dqu_vals != cond_prob] # don't need to rerun

clust_dqu_plots_total <- lapply(dqu_vals, \(dqu) {
  print(paste0("dqu = ", dqu))
  clust_plot_dqu(dist_mat = "Total", dqu = dqu)
})

pdf(
  paste0(save_dir, "cluster_maps_dqu_total.pdf"),
  width = 12,
  height = 10
)

for (i in seq_along(dqu_vals)) {
  print(
    clust_dqu_plots_total[[i]]$map_plots +
      plot_annotation(
        title = paste0("DQU = ", dqu_vals[[i]])
      )
  )
}

dev.off()


pdf(
  paste0(save_dir, "heatmaps_dqu_total.pdf"),
  width = 12,
  height = 10
)

for (i in seq_along(dqu_vals)) {
  print(
    clust_dqu_plots_total[[i]]$heatmap_plots +
      plot_annotation(
        title = paste0("DQU = ", dqu_vals[[i]])
      )
  )
}

dev.off()

rm(clust_dqu_plots_total)
gc()

# var_dep and temp
clust_dqu_plots_var_dep <- lapply(dqu_vals, \(dqu) {
  clust_plot_dqu(dist_mat = dep_vars, dqu = dqu)
})
# pdf(paste0(save_dir, "cluster_maps_dqu_rain.pdf"), width = 12, height = 10)
pdf(paste0(save_dir, "cluster_maps_dqu_", dep_vars, ".pdf"), width = 12, height = 10)
for (i in seq_along(dqu_vals)) {
  print(
    clust_dqu_plots_var_dep[[i]]$map_plots +
      plot_annotation(
        title = paste0("DQU = ", dqu_vals[[i]])
      )
  )
}

rm(clust_dqu_plots_var_dep)
gc()

clust_dqu_plots_temp <- lapply(dqu_vals, \(dqu) {
  clust_plot_dqu(dist_mat = "temp", dqu = dqu)
})
pdf(paste0(save_dir, "cluster_maps_dqu_temp.pdf"), width = 12, height = 10)
for (i in seq_along(dqu_vals)) {
  print(
    clust_dqu_plots_temp[[i]]$map_plots +
      plot_annotation(
        title = paste0("DQU = ", dqu_vals[[i]])
      )
  )
}
dev.off()

rm(clust_dqu_plots_temp)
gc()


#### Fit & Cluster for each decade ####

# Fit CE model for each variable and season
dep_var_season_decade_deep <- lapply(dep_vars, \(spec_var) {
  lapply(marg_season_decade, \(x) {
    cecl_dep(
      x,
      vars = c("temp", spec_var),
      cond_vars = c("temp", spec_var),
      # cond_prob = 0.9,
      cond_prob = cond_prob,
      nruns = 3,
      ncores = getOption("mc_cores", 2L)
    )
  })
})
names(dep_var_season_decade_deep) <- dep_vars

# flatten deep list
dep_var_season_decade <- purrr::list_flatten(
  dep_var_season_decade_deep,
  name_spec = "{outer}_{inner}"
)

# dataframe of list names for later
name_df_decade <- tibble(name = names(dep_var_season_decade)) |>
  tidyr::separate_wider_regex(
    name,
    c(
      var = ".*",
      "_",
      season = "[^_]+"
    )
  ) |>
  tidyr::separate_wider_regex(
    var,
    c(
      var = ".*",
      "_",
      decade = "[^_]+"
    )
  )

# calculate distances
set.seed(123)
dist_var_season_decade <- mclapply(seq_along(dep_var_season_decade), \(i) {
  cecl_dist(
    dep_var_season_decade[[i]],
    marg_season_decade[[with(name_df_decade[i, ], paste0(decade, "_", season))]],
    ncores = getOption("mc_cores", 2L)
  )
})
names(dist_var_season_decade) <- names(dep_var_season_decade)

# Cluster for each variable and season
clust_var_season_decade_k <- lapply(k_vec, \(k) {
  print(paste0("k = ", k))
  ret <- mclapply(seq_along(dep_var_season_decade), \(i) {
    print(paste0("i = ", i))
    cecl_clust(
      dist_var_season_decade[[i]],
      k        = k,
      ncores   = getOption("mc_cores", 2L)
    )
  })
  names(ret) <- names(dep_var_season_decade)
  ret
})

names(clust_var_season_decade_k) <- k_vec

# flatten deep list
clust_var_season_decade <- purrr::list_flatten(
  clust_var_season_decade_k,
  name_spec = "{inner}_{outer}"
)

titles <- stringr::str_replace_all(names(clust_var_season_decade), "_", " ")
map_plots_decade <- mclapply(seq_along(clust_var_season_decade), \(i) {
  plt_clust_map(
    pts,
    areas_ccaa,
    # clust_var_season_k[[i]][[j]]$pam,
    clust_var_season_decade[[i]]$pam,
    plot_medoids = FALSE,
    pt_size = 6
  ) +
    # ggtitle(paste0(name_df$var[j], " - ", name_df$season[j]))
    ggtitle(titles[[i]])
})

# Also plot heatmaps
heatmap_plots_decade <- mclapply(seq_along(clust_var_season_decade), \(i) {
  ggplot(clust_var_season_decade[[i]], which = "image") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(titles[[i]])
})

names(map_plots_decade) <- names(heatmap_plots_decade) <- names(clust_var_season_decade)

# TODO Save and comment

#### Differences per decade ####

name_df_decade_k <- crossing(
  k = k_vec,
  name_df_decade
) |>
  relocate(k, .after = everything())

# double check everything matches up
all(
  with(name_df_decade_k, paste(var, decade, season, k, sep = "_")) ==
    names(clust_var_season_decade)
)

# pull clustering solutions for each
clust_sol_lst <- lapply(clust_var_season_decade, coef)

# join in information form name_df_decade_k
clust_sol_df <- bind_rows(lapply(seq_along(clust_sol_lst), \(i) {
  cbind(clust_sol_lst[[i]], name_df_decade_k[i, , drop = FALSE])
}))

# See if clustering has changed since 1960 and since previous decade
clust_sol_df_diff <- clust_sol_df %>%
  arrange(var, season, k, name, decade) %>%
  group_by(var, season, k, name) %>%
  mutate(
    cluster_1960 = cluster[decade == 1960][1],
    cluster_prev = lag(cluster),
    changed_since_1960 = (cluster != cluster_1960),
    changed_since_prev = (cluster != cluster_prev)
  ) %>%
  ungroup() |>
  mutate(changed_since_prev = ifelse(
    is.na(changed_since_prev), FALSE, changed_since_prev
  ))

# function to plot with
plot_diff_map <- \(x, pts, areas, decade, season, var, k) {
  false_pt_size <- ifelse(decade == 1960, 6, 4)

  p <- left_join(pts, x) |>
    filter(
      # decade == 1970,
      # season == "Winter",
      # var == "rain",
      # k == 3
      decade == !!decade,
      season == !!season,
      var == !!var,
      k == !!k
    ) |>
    ggplot() +
    geom_sf(
      data = areas,
      fill = NA,
      colour = "black"
    ) +
    coord_sf(expand = FALSE) +
    geom_sf(
      aes(
        fill = factor(cluster),
        shape = changed_since_1960,
        size = changed_since_prev
      ),
      alpha = 0.8,
      stroke = 1,
      colour = "black"
    ) +
    scale_fill_brewer(
      palette = "Set2",
      name = "Cluster"
    ) +
    scale_shape_manual(
      values = c(
        `FALSE` = 21,
        `TRUE` = 24
      ),
      labels = c(
        `FALSE` = "Same",
        `TRUE` = "Different"
      ),
      name = "vs 1960"
    ) +
    scale_size_manual(
      values = c(
        # `FALSE` = 4,
        `FALSE` = false_pt_size,
        `TRUE` = 7
      ),
      labels = c(
        `FALSE` = "No",
        `TRUE` = "Yes"
      ),
      name = "vs last decade"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          size = 5
        )
      ),
      shape = guide_legend(order = 2),
      size = guide_legend(order = 3)
    ) +
    guides(fill = "none") +
    cecl_theme()

  if (decade == 1960) {
    p <- p + guides(shape = "none", size = "none")
  }
  p
}


# loop over decades, plotting and saving individually for each season + k vals
# NOTE: Too memory intensive to run with mclapply!
diff_map_plots <- lapply(k_vec, \(k_i) {
  gc()
  lapply(decades, \(decade_i) {
    # lapply(c(1960, 1970), \(decade_i) {
    gc()
    season_plots <- lapply(seasons, \(season_i) {
      gc()
      p <- plot_diff_map(
        x = clust_sol_df_diff,
        pts = pts,
        areas = areas_ccaa,
        decade = decade_i,
        season = season_i,
        var = "rain",
        k = k_i
      ) +
        ggtitle(paste0("rain - ", season_i))
      if (season_i != "Winter") {
        p <- p + guides(shape = "none", size = "none")
      }
      p
    })


    wrap_plots(season_plots, nrow = 1) +
      plot_annotation(
        title = paste0("k = ", k_i, ", decade = ", decade_i),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
      )
  })
})

out_dir <- "plots/02_app/01_clust_diff_maps"
fs::dir_create(out_dir)

walk2(
  diff_map_plots,
  k_vec,
  \(plots_for_k, k_i) {
    walk2(
      plots_for_k,
      decades,
      \(p, decade_i) {
        ggsave(
          filename = glue::glue(
            "{out_dir}/diff_map_rain_k{k_i}_{decade_i}.png"
          ),
          plot = p,
          width = 12,
          height = 4,
          dpi = 300
        )
      }
    )
  }
)

# Now plot differences in dist matrices

# output folder
out_dir <- "plots/02_app/02_divergence_diff_heatmaps"
out_dir_1960 <- "plots/02_app/021_div_since_1960"
out_dir_prev <- "plots/02_app/022_div_since_previous"

fs::dir_create(out_dir_1960)
fs::dir_create(out_dir_prev)

# helper to get object names
make_dist_name <- \(var, decade, season) {
  paste(var, decade, season, sep = "_")
}

# function to plot difference from 1960
plot_divergence_diff <- \(dist_list, var, decade, season, k, base_year = 1960) {
  current_name <- make_dist_name(var, decade, season)
  base_name <- make_dist_name(var, base_year, season)

  current_obj <- dist_list[[current_name]]
  base_obj <- dist_list[[base_name]]

  diff_obj <- current_obj
  diff_obj$dist_mat <- current_obj$dist_mat - base_obj$dist_mat

  clust_obj <- cecl_clust(
    diff_obj,
    k = k
  )

  # p <- ggplot(
  #   clust_obj,
  #   which = "image"
  # )
  #
  # rng <- max(abs(p$data$Freq), na.rm = TRUE)
  #
  rng <- max(abs(diff_obj$dist_mat), na.rm = TRUE)

  ggplot(
    clust_obj,
    which = "image",
    show_xlab = TRUE,
    show_ylab = TRUE
  ) +
    ggtitle(paste0(var, " - ", season)) +
    scale_fill_gradient2(
      low = "blue3",
      mid = "white",
      high = "red3",
      midpoint = 0,
      limits = c(-rng, rng),
      na.value = "grey80"
    ) +
    theme(axis.text.x = element_blank())
}

plot_divergence_diff(dist_var_season_decade, "rain", 2000, "Winter", 3, 1960)

previous_decade <- \(decade) decade - 10

div_plots_since_1960 <- lapply(k_vec, \(k_i) {
  lapply(decades, \(decade_i) {
    season_plots <- lapply(seasons, \(season_i) {
      plot_divergence_diff(
        dist_list = dist_var_season_decade,
        var = "rain",
        decade = decade_i,
        season = season_i,
        k = k_i,
        base_year = 1960
      )
    })

    p <- wrap_plots(season_plots, nrow = 1) +
      plot_annotation(
        title = paste0("Change since 1960: k = ", k_i, ", decade = ", decade_i),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
      )

    ggsave(
      filename = glue("{out_dir_1960}/divergence_since_1960_rain_k{k_i}_{decade_i}.png"),
      plot = p,
      width = 14,
      height = 5,
      dpi = 300
    )

    p
  })
})

div_plots_since_prev <- lapply(k_vec, \(k_i) {
  lapply(decades[decades != 1960], \(decade_i) {
    base_i <- decade_i - 10

    season_plots <- lapply(seasons, \(season_i) {
      plot_divergence_diff(
        dist_list = dist_var_season_decade,
        var = "rain",
        decade = decade_i,
        season = season_i,
        k = k_i,
        base_year = base_i
      )
    })

    p <- wrap_plots(season_plots, nrow = 1) +
      plot_annotation(
        title = paste0("Change since ", base_i, ": k = ", k_i, ", decade = ", decade_i),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
      )

    ggsave(
      filename = glue("{out_dir_prev}/divergence_since_previous_rain_k{k_i}_{base_i}_to_{decade_i}.png"),
      plot = p,
      width = 14,
      height = 5,
      dpi = 300
    )

    p
  })
})

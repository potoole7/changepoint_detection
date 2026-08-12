#### Changepoint Detection for Spanish application ####

# Things to check from writing retreat!
# TODO - When does it happen that # years on LHS > RHS, or vice versa??
# TODO -

# Save plots from screening (done)
# TODO Write up summary of screening results
# TODO Run permutation tests with 100 permutations
# TODO Run permutation tests with 1000 permutations (??)
# TODO Code changes required to allow for multiple testing


# TODO Think about effect of different hyperparameters (n_perm, n_per_block, etc)

# How to deal with different number of observations for some locations?? (done)

# Maybe add progress in perm_test_fun? Can't really tell how things are# going currently (done)

# Look into Laplace sample being used across different comparisons, how
# to decide? Is the 80th quantile of the Laplace distribution a good shout? (done)

# Misc:
# TODO What about which model to use, how do I decide on that? Use total for
# now??

# TODO Move functions to own section

# Sliding window:
# Get working (done)
# Investigate errors for 10 years in August (done, not enough years)
# Investigate massive dissimilarities (done, occurs where not enough data)
# Identify peak(s) for each season (done)
# Try for different #s of min years (10, 15, 20 I guess? For each season) (done)
# Plot (done)

# Permutation test
# TODO Get working
# TODO Identify significant change points for each season

# TODO Test for 10, 15 and 20 year block sizes
# TODO Increase n_perm to 1000, ideally


#### libs ####

devtools::load_all("../CeCl")
library(grid)
library(lubridate)
library(RColorBrewer)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(purrr)
library(parallel)
library(evgam)
library(ggridges)
library(sf)
library(patchwork)

# source custom functions
source("src/00_functions.R")


#### Functions ####

# Function to perform screening for a single season and number of years
# Screen one season using complete season-year windows.
screen_one_setting <- \(
  season_name,
  n_years_per_block,
  verbose = TRUE,
  mc.cores = getOption("mc.cores", 1L)
) {
  # Prepare data
  season_data <- data_laplace_season[[season_name]] |>
    select(
      name,
      date,
      season_year,
      X1 = temp,
      X2 = drought_local
    ) |>
    filter(
      !is.na(name),
      !is.na(date),
      !is.na(season_year)
    ) |>
    arrange(name, date)

  available_years <- season_data |>
    distinct(season_year) |>
    arrange(season_year) |>
    pull(season_year)

  n_available_years <- length(available_years)

  if (n_available_years < 2L * n_years_per_block) {
    stop(
      "Season ", season_name,
      " has only ", n_available_years,
      " years, but ", 2L * n_years_per_block,
      " are required."
    )
  }

  candidate_positions <- seq.int(
    from = n_years_per_block,
    to = n_available_years - n_years_per_block
  )

  block_vec <- available_years[candidate_positions]

  # Try to obtain common starting values
  middle_candidate <- block_vec[
    ceiling(length(block_vec) / 2)
  ]

  start_fit <- tryCatch(
    single_run_explore(
      data = season_data,
      i = middle_candidate,
      n_years_per_block = n_years_per_block,
      n_vars = 2,
      laplace_trans = TRUE,
      skip_failed = TRUE,
      cond_val = dep_val,
      laplace_sample = laplace_sample,
      aLow = 0,
      nruns = 3,
      ret_dep = TRUE
    ),
    error = \(e) {
      NULL
    }
  )

  start_vals <- NULL

  if (
    is.list(start_fit) &&
      isTRUE(start_fit$success) &&
      !is.null(start_fit$dep)
  ) {
    candidate_start_vals <- tryCatch(
      lapply(start_fit$dep, coef),
      error = \(e) {
        NULL
      }
    )

    # valid_start_vals <- (
    #   !is.null(candidate_start_vals) &&
    #     length(candidate_start_vals) > 0L &&
    #     all(
    #       vapply(
    #         candidate_start_vals,
    #         \(x) {
    #           is.numeric(x) &&
    #             length(x) > 0L &&
    #             all(is.finite(x))
    #         },
    #         logical(1)
    #       )
    #     )
    # )
    valid_start_vals <- (
      !is.null(candidate_start_vals) &&
        length(candidate_start_vals) == 2L &&
        all(
          vapply(
            candidate_start_vals,
            \(x) {
              !is.null(x) && length(x) > 0L
            },
            logical(1)
          )
        )
    )

    if (valid_start_vals) {
      start_vals <- candidate_start_vals

      if (verbose) {
        message(
          "Season ", season_name,
          ", ", n_years_per_block,
          "-year windows: starting values obtained at candidate ",
          middle_candidate,
          "."
        )
      }
    }
  }

  if (is.null(start_vals)) {
    warning(
      "No valid common starting values for season ",
      season_name,
      " with ", n_years_per_block,
      "-year windows; screening candidates using default starts.",
      call. = FALSE
    )
  }

  # Screen individual candidates
  run_candidate <- \(candidate_index) {
    candidate_year <- block_vec[[candidate_index]]

    fit_args <- list(
      data = season_data,
      i = candidate_year,
      n_years_per_block = n_years_per_block,
      n_vars = 2,
      laplace_trans = TRUE,
      skip_failed = TRUE,
      cond_val = dep_val,
      laplace_sample = laplace_sample,
      aLow = 0
    )

    # Common starting values are optional.
    if (!is.null(start_vals)) {
      fit_args$start <- start_vals
    }

    result <- tryCatch(
      do.call(
        single_run_explore,
        fit_args
      ),
      error = \(e) {
        list(
          success = FALSE,
          block_info = list(
            split_after = candidate_year,
            split_before = available_years[
              candidate_positions[[candidate_index]] + 1L
            ],
            n_left = NA_integer_,
            n_right = NA_integer_
          ),
          frob = NA_real_,
          inf = NA_real_,
          spec = NA_real_,
          dep = NULL,
          error_stage = "single_run_explore",
          error_message = conditionMessage(e)
        )
      }
    )

    if (!is.list(result) || is.null(result$success)) {
      result <- list(
        success = FALSE,
        block_info = list(
          split_after = candidate_year,
          split_before = available_years[
            candidate_positions[[candidate_index]] + 1L
          ],
          n_left = NA_integer_,
          n_right = NA_integer_
        ),
        frob = NA_real_,
        inf = NA_real_,
        spec = NA_real_,
        dep = NULL,
        error_stage = "invalid_result",
        error_message = paste0(
          "`single_run_explore()` returned an invalid result at candidate ",
          candidate_year,
          "."
        )
      )
    }

    result
  }

  if (mc.cores > 1L && .Platform$OS.type != "windows") {
    screen_vals <- parallel::mclapply(
      seq_along(block_vec),
      run_candidate,
      mc.cores = mc.cores
    )
  } else {
    screen_vals <- lapply(
      seq_along(block_vec),
      run_candidate
    )
  }

  # Collate results
  screen_res <- bind_rows(
    lapply(
      seq_along(screen_vals),
      \(candidate_index) {
        result <- screen_vals[[candidate_index]]
        block_info <- result$block_info

        tibble(
          season = season_name,
          n_years_per_block = n_years_per_block,
          change_after_year = block_vec[[candidate_index]],
          change_before_year = if (
            is.list(block_info) &&
              !is.null(block_info$split_before)
          ) {
            as.integer(block_info$split_before)
          } else {
            available_years[
              candidate_positions[[candidate_index]] + 1L
            ]
          },
          n_left_dates = if (
            is.list(block_info) &&
              !is.null(block_info$n_left)
          ) {
            as.integer(block_info$n_left)
          } else {
            NA_integer_
          },
          n_right_dates = if (
            is.list(block_info) &&
              !is.null(block_info$n_right)
          ) {
            as.integer(block_info$n_right)
          } else {
            NA_integer_
          },
          success = isTRUE(result$success),
          frob = if (isTRUE(result$success)) {
            as.numeric(result$frob)
          } else {
            NA_real_
          },
          inf = if (isTRUE(result$success)) {
            as.numeric(result$inf)
          } else {
            NA_real_
          },
          spec = if (isTRUE(result$success)) {
            as.numeric(result$spec)
          } else {
            NA_real_
          },
          error_stage = if (
            is.null(result$error_stage)
          ) {
            NA_character_
          } else {
            as.character(result$error_stage)
          },
          error_message = if (
            is.null(result$error_message)
          ) {
            NA_character_
          } else {
            as.character(result$error_message)
          }
        )
      }
    )
  ) |>
    arrange(change_after_year)

  # Mark local Frobenius peaks
  screen_res <- screen_res |>
    mutate(
      local_peak_frob = (
        success &
          lag(success, default = FALSE) &
          lead(success, default = FALSE) &
          frob > lag(frob) &
          frob >= lead(frob)
      ),
      local_peak_frob = replace_na(
        local_peak_frob,
        FALSE
      )
    )

  failure_rate <- mean(!screen_res$success)

  if (verbose) {
    message(
      "Season ", season_name,
      ", ", n_years_per_block,
      "-year windows: ",
      sum(screen_res$success), " of ", nrow(screen_res),
      " candidates succeeded; failure rate = ",
      scales::percent(
        failure_rate,
        accuracy = 0.1
      ),
      "."
    )
  }

  screen_res
}

# Run the permutation scan for one season
run_season_permutation_scan <- \(
  season_name,
  n_years_per_block = 25L,
  n_perm = 100L,
  seed = 123L,
  use_start = TRUE,
  ret_dep = TRUE,
  verbose = TRUE
) {
  season_data <- data_laplace_season[[season_name]] |>
    select(
      name,
      date,
      season_year,
      X1 = temp,
      X2 = drought_local
    ) |>
    filter(
      !is.na(name),
      !is.na(date),
      !is.na(season_year)
    ) |>
    arrange(name, date)

  available_years <- season_data |>
    distinct(season_year) |>
    arrange(season_year) |>
    pull(season_year)

  n_available_years <- length(available_years)

  if (n_available_years < 2L * n_years_per_block) {
    stop(
      "Season ", season_name,
      " has only ", n_available_years,
      " seasonal years; ",
      2L * n_years_per_block,
      " are required."
    )
  }

  candidate_positions <- seq.int(
    from = n_years_per_block,
    to = n_available_years - n_years_per_block
  )

  candidate_years <- available_years[
    candidate_positions
  ]

  message(
    "Running ", n_perm,
    " permutations for season ", season_name,
    " over candidate years ",
    min(candidate_years),
    "–",
    max(candidate_years),
    " using ", n_years_per_block,
    " years per block."
  )

  set.seed(seed)

  permutation_results <- perm_test_fun(
    data = season_data,
    grid_vals = candidate_years,
    n_perm = n_perm,
    max_perm_attempts = 5L * n_perm,
    max_failure_rate = 0.10,
    n_years_per_block = n_years_per_block,
    n_vars = 2,
    laplace_trans = TRUE,
    use_start = use_start,
    ret_dep = ret_dep,
    use_dth = FALSE,
    verbose = verbose,
    cond_val = dep_val,
    laplace_sample = laplace_sample,
    n_mc = length(laplace_sample),
    aLow = 0,
    nruns = 1,
    validation_warnings = TRUE,
    permutation_validation_warnings = FALSE
  )

  list(
    season = season_name,
    n_years_per_block = n_years_per_block,
    n_perm = n_perm,
    candidate_years = candidate_years,
    results = permutation_results
  )
}



#### Metadata ####

dep_var <- c("drought_local_rev")
temp_var <- "temp_max"
decades <- seq(1960, 2010, by = 10)
seasons <- c("Winter", "Spring", "Summer", "Autumn")

seed <- 123 # random seed
# Conditional threshold and number of samples for Laplace sample used throughout
# dqu <- 0.8
dqu <- 0.85
n_samples <- 500
# grid_vals <- seq(950, 1050, by = 10) # TODO Add more values once code works

# run initially for just 100 permutations across full range
# TODO Change names, no need to call "screen", right??
# n_perm_screen <- 100L
# n_perm_screen <- 1000L
n_perm_screen <- 200L
n_years_per_block <- 25L # best choice, from screening round


#### Load Data ####

# data
data <- readr::read_csv(
  "data/02_app/ecad_clean.csv.gz"
) |>
  mutate(decade = factor(floor(year(date) / 10) * 10, levels = decades))

if (dep_var == "rain") {
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
  temp_var <- "temp"
}

# reverse drought_local variable to give positive alpha values, if desired
if (dep_var == "drought_local_rev") {
  data <- data |>
    mutate(drought_local = -drought_local)
  dep_var <- c("drought_local")
}

# check percentage of dates
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


marg_season <- readRDS("data/02_app/marg_season_roll_emp.rds")

# pull Laplace transformed data from marginal object, and join together
# data_laplace_season <- lapply(marg_season, \(x) {
#   y <- x$transformed
#
#   bind_rows(lapply(seq_along(y), \(i) {
#     as.data.frame(y[[i]]) |>
#       mutate(name = names(y)[[i]])
#   }))
# })
data_laplace_season <- lapply(marg_season, \(x) {
  bind_rows(Map(
    \(transformed, dates, station) {
      stopifnot(nrow(transformed) == length(dates))

      as.data.frame(transformed) |>
        mutate(
          date = as.Date(dates),
          name = station,
          .before = 1
        )
    },
    x$transformed,
    x$dates,
    names(x$transformed)
  ))
})

# confirm alignment by season between all locations
bind_rows(lapply(seq_along(data_laplace_season), \(i) {
  data_laplace_season[[i]] |>
    mutate(season = names(data_laplace_season)[[i]])
})) |>
  group_by(season, name) |>
  summarise(
    n_rows = n(),
    n_dates = n_distinct(date),
    first_date = min(date),
    last_date = max(date)
  ) |>
  distinct(season, n_rows, n_dates, first_date, last_date)

# check all stations have same dates
bind_rows(lapply(seq_along(data_laplace_season), \(i) {
  data_laplace_season[[i]] |>
    mutate(season = names(data_laplace_season)[[i]])
})) |>
  # group_by(season) |>
  distinct(name, date, season) |>
  group_by(season) |>
  count(date) |>
  count(n)

# add season_date
data_laplace_season <- lapply(
  names(data_laplace_season),
  \(s) {
    data_laplace_season[[s]] |>
      mutate(
        season = s,
        season_year = if_else(
          season == "Winter" & month(date) == 12L,
          year(date) + 1L,
          year(date)
        )
      )
  }
) |>
  setNames(names(marg_season))

# load map of continental Spain to use as background in plots
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


#### Precalculations ####

# use the same dependence threshold across all variables and locations

dep_val <- qlaplace(dqu) # for Laplace marginals

# Calculate cap for sampling from Laplace distribution
laplace_cap <- qlaplace(0.99)

# function to generate Laplace samples
rlaplace_trunc <- \(n, thresh_max = qlaplace(0.8), y_max = qlaplace(0.99)) {
  # get maximum point
  stopifnot(
    "y_max must be greater than thresh_max" = y_max > thresh_max
  )
  # get probability of being below this point from exponential CDF
  p_max <- 1 - exp(-(y_max - thresh_max))
  # sample from uniform distribution below this point
  U <- stats::runif(n, min = 0, max = p_max) # min=0 as we push up by thresh
  # inversion sampling from exponential distribution
  W <- -log(1 - U)

  # shift to the right by the threshold to get samples from truncated Laplace
  return(thresh_max + W)
}

set.seed(seed)
laplace_sample <- rlaplace_trunc(
  n = n_samples, # TODO Increase later ??
  thresh_max = dep_val,
  y_max = laplace_cap
)

# Also make a longer Laplace sample
# laplace_sample2 <- rlaplace_trunc(
#   n = 2000,
#   thresh_max = dep_val,
#   y_max = laplace_cap
# )

#### Screening ####

screen_setup_df <- tidyr::crossing(
  "season" = seasons,
  # 10 years too short (as is 15 actually) !!
  # n_years_per_block = c(10L, 15L, 20L, 25L)
  n_years_per_block = c(15L, 20L, 25L, 30L)
)

# sink("sink_screening_output.txt")
sink(paste0("sink_screening_output_dqu_", dqu, ".txt"))

screen_res_df <- bind_rows(lapply(seq_len(nrow(screen_setup_df)), \(i) {
  print(paste0(round(i / nrow(screen_setup_df) * 100, 2), "% of setups done"))
  print(paste0("season = ", screen_setup_df$season[[i]]))
  print(paste0("n_years_per_block = ", screen_setup_df$n_years_per_block[[i]]))
  if (i == 2) {
    # debugonce(single_run_explore)
  }
  with(
    screen_setup_df,
    screen_one_setting(season[[i]], n_years_per_block[[i]])
  )
}))

warnings()

# summarise failures by season and n_years_per_block
screen_failure_summary <- screen_res_df |>
  group_by(
    season,
    n_years_per_block
  ) |>
  summarise(
    n_candidates = n(),
    n_successful = sum(success),
    n_failed = sum(!success),
    failure_rate = mean(!success),
    .groups = "drop"
  )

screen_failure_summary

sink()

# save
# readr::write_csv(screen_res_df, "data/02_app/screen_res.csv.gz")
readr::write_csv(screen_res_df, paste0("data/02_app/screen_res_dqu_", dqu, ".csv.gz"))

# screen_res_df <- readr::read_csv("data/02_app/screen_res.csv.gz")
screen_res_df <- readr::read_csv(paste0("data/02_app/screen_res_dqu_", dqu, ".csv.gz"))


#### Plotting ####

screen_res_df_plt <- screen_res_df |>
  arrange(
    season,
    n_years_per_block,
    change_after_year
  ) |>
  group_by(
    season,
    n_years_per_block
  ) |>
  mutate(
    # Separate successful line segments at failed candidates.
    success_run = cumsum(!success),

    # Identify Frobenius local peaks.
    local_peak_frob = (
      success &
        lag(success, default = FALSE) &
        lead(success, default = FALSE) &
        frob > lag(frob) &
        frob >= lead(frob)
    ),
    local_peak_frob = replace_na(
      local_peak_frob,
      FALSE
    )
  ) |>
  ungroup() |>
  mutate(
    setting = paste0(
      season,
      ", ",
      n_years_per_block,
      " years"
    )
  )

# Preserve the desired facet ordering.
setting_levels <- screen_res_df_plt |>
  distinct(
    season,
    n_years_per_block,
    setting
  ) |>
  arrange(
    season,
    n_years_per_block
  ) |>
  pull(setting)

screen_res_df_plt <- screen_res_df_plt |>
  mutate(
    setting = factor(
      setting,
      levels = setting_levels
    )
  )


# summarise Frobenius peaks
candidate_peaks_frob <- screen_res_df_plt |>
  filter(
    success,
    local_peak_frob
  ) |>
  arrange(
    season,
    n_years_per_block,
    desc(frob)
  )

top_local_peaks_frob <- candidate_peaks_frob |>
  group_by(
    season,
    n_years_per_block
  ) |>
  slice_max(
    order_by = frob,
    n = 3,
    with_ties = FALSE
  ) |>
  arrange(
    season,
    n_years_per_block,
    desc(frob)
  ) |>
  ungroup()

top_local_peaks_frob


# plot Frobenius norm peaks
year_breaks <- seq(
  floor(
    min(
      screen_res_df_plt$change_after_year,
      na.rm = TRUE
    ) / 2
  ) * 2,
  ceiling(
    max(
      screen_res_df_plt$change_after_year,
      na.rm = TRUE
    ) / 2
  ) * 2,
  by = 2
  # by = 1
)

plot_frob <- \(df, spec_year = NULL) {
  df_plot <- df
  if (!is.null(spec_year)) {
    df_plot <- df_plot |>
      filter(n_years_per_block == spec_year)
  }
  df_plot |>
    ggplot(
      aes(
        x = change_after_year,
        y = frob
      )
    ) +
    geom_line(
      data = \(x) filter(x, success),
      aes(
        group = interaction(
          setting,
          success_run
        )
      ),
      colour = "black"
    ) +
    geom_point(
      data = \(x) filter(x, success),
      colour = "black"
    ) +
    geom_point(
      data = \(x) filter(x, local_peak_frob),
      colour = "red",
      size = 3
    ) +
    geom_rug(
      data = \(x) filter(x, !success),
      aes(x = change_after_year),
      inherit.aes = FALSE,
      sides = "b",
      colour = "grey50"
    ) +
    facet_wrap(
      ~setting,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = year_breaks
    ) +
    labs(
      # x = "Candidate change after seasonal year",
      x = "season year",
      y = "Frobenius discrepancy",
      # caption = paste(
      #   "Red points indicate local peaks.",
      #   "Grey axis marks indicate failed candidates."
      # )
    ) +
    cecl_theme() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

p_frob <- plot_frob(screen_res_df_plt)
p_frob_25 <- plot_frob(screen_res_df_plt, spec_year = 25)

# Convert all norms to long format
screen_res_long <- screen_res_df_plt |>
  pivot_longer(
    cols = c(
      frob,
      inf,
      spec
    ),
    names_to = "norm",
    values_to = "value"
  ) |>
  mutate(
    norm = recode(
      norm,
      frob = "Frobenius",
      inf = "Maximum",
      spec = "Spectral"
    ),
    norm = factor(
      norm,
      levels = c(
        "Frobenius",
        "Maximum",
        "Spectral"
      )
    )
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    change_after_year
  ) |>
  group_by(
    season,
    n_years_per_block,
    norm
  ) |>
  mutate(
    # Peaks are calculated separately for each norm.
    local_peak = (
      success &
        lag(success, default = FALSE) &
        lead(success, default = FALSE) &
        value > lag(value) &
        value >= lead(value)
    ),
    local_peak = replace_na(
      local_peak,
      FALSE
    )
  ) |>
  ungroup()


# Peak summaries for every norm
candidate_peaks_all <- screen_res_long |>
  filter(
    success,
    local_peak
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    desc(value)
  )

top_local_peaks_all <- candidate_peaks_all |>
  group_by(
    season,
    n_years_per_block,
    norm
  ) |>
  slice_max(
    order_by = value,
    n = 3,
    with_ties = FALSE
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    desc(value)
  ) |>
  ungroup()

# Combined plot for all norms
plot_all_norms <- \(df, spec_year = NULL) {
  df_plot <- df
  if (!is.null(spec_year)) {
    df_plot <- df_plot |>
      filter(n_years_per_block == spec_year)
  }
  df_plot |>
    # filter(n_years_per_block == 25) |>
    ggplot(
      aes(
        x = change_after_year,
        y = value,
        colour = norm
      )
    ) +
    geom_line(
      data = \(x) filter(x, success),
      aes(
        group = interaction(
          setting,
          norm,
          success_run
        )
      )
    ) +
    geom_point(
      data = \(x) filter(x, success)
    ) +
    geom_point(
      data = \(x) filter(x, local_peak),
      colour = "red",
      size = 3
    ) +
    geom_rug(
      data = \(x) filter(x, !success),
      aes(x = change_after_year),
      inherit.aes = FALSE,
      sides = "b",
      colour = "grey50"
    ) +
    facet_wrap(
      ~setting,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = year_breaks
    ) +
    scale_colour_brewer(
      palette = "Dark2"
    ) +
    labs(
      # x = "Candidate change after seasonal year",
      x = "season year",
      y = "Discrepancy",
      colour = "Norm",
      # caption = paste(
      #   "Red points indicate norm-specific local peaks.",
      #   "Grey axis marks indicate failed candidates."
      # )
    ) +
    cecl_theme() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )
}

p_all_norms <- plot_all_norms(screen_res_long)
p_all_norms_25 <- plot_all_norms(screen_res_long, spec_year = 25)

p_frob
p_frob_25

p_all_norms
p_all_norms_25

ggsave(paste0("plots/02_app/p_frob_dqu_", dqu, ".png"), p_frob, width = 12, height = 8)
ggsave(paste0("plots/02_app/p_frob_25_dqu_", dqu, ".png"), p_frob_25, width = 12, height = 8)
ggsave(paste0("plots/02_app/p_all_norms_dqu_", dqu, ".png"), p_all_norms, width = 12, height = 8)
ggsave(paste0("plots/02_app/p_all_norms_25_dqu_", dqu, ".png"), p_all_norms_25, width = 12, height = 8)


#### Permutation tests across full feasible range ####

# sink(paste0(
#   "sink_perm_output_nperm_",
#   n_perm_screen, "_n_years_",
#   n_years_per_block, ".txt"
# ))
sink(paste0(
  "sink_perm_output_nperm_",
  n_perm_screen, "_dqu_",
  dqu, ".txt"
))

source("src/00_functions.R")
permutation_scan_results <- setNames(
  lapply(
    seq_along(seasons),
    \(i) {
      message("Season = ", seasons[[i]])
      run_season_permutation_scan(
        season_name = seasons[[i]],
        n_years_per_block = n_years_per_block,
        n_perm = n_perm_screen,
        seed = seed + i - 1L,
        use_start = TRUE,
        ret_dep = TRUE,
        verbose = TRUE
      )
    }
  ),
  seasons
)

sink()

saveRDS(
  permutation_scan_results,
  paste0(
    "data/02_app/",
    "spain_perm_test_",
    "n_perm_",
    n_perm_screen,
    "_dqu_",
    dqu,
    ".rds"
  )
)

permutation_scan_results <- readRDS(
  paste0(
    "data/02_app/",
    "spain_perm_test_",
    "n_perm_",
    n_perm_screen,
    "_dqu_",
    dqu,
    ".rds"
  )
)


#### Plotting ####

# preprocess every season
permutation_scan_tidy <- lapply(
  permutation_scan_results,
  preprocess_permutation_scan
)

permutation_summary_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$summary
  )
)

permutation_values_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$permutations
  )
)

permutation_dependence_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$dependence
  )
)

# check failures
permutation_summary_df |>
  distinct(
    season,
    change_after_year,
    success,
    n_attempted,
    n_successful,
    n_failed,
    failure_rate,
    error_stage,
    error_message
  ) |>
  arrange(
    season,
    change_after_year
  )

# create plots for every season
permutation_scan_plots <- lapply(
  permutation_scan_tidy,
  plot_permutation_scan,
  plot_ce_parameters = TRUE,
  variable_names = c(
    X1 = "Maximum temperature",
    X2 = "Drought"
  )
)

# save
# plot_directory <- "plots/permutation_scan"
plot_directory <- paste0("plots/permutation_scan_n_perm_", n_perm_screen)

dir.create(
  plot_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

for (season_name in names(permutation_scan_plots)) {
  season_plots <-
    permutation_scan_plots[[season_name]]

  filename_prefix <- paste0(
    tolower(season_name),
    "_",
    n_years_per_block,
    "yr_",
    n_perm_screen,
    "perm"
  )

  ggsave(
    filename = file.path(
      plot_directory,
      paste0(
        filename_prefix,
        "_p_value_profile.png"
      )
    ),
    plot = season_plots$p_value_profile,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      plot_directory,
      paste0(
        filename_prefix,
        "_observed_vs_permuted.png"
      )
    ),
    plot =
      season_plots$permutation_boxplots,
    width = 11,
    height = 7,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      plot_directory,
      paste0(
        filename_prefix,
        "_permutation_histograms.png"
      )
    ),
    plot =
      season_plots$permutation_histograms,
    width = 12,
    height = 4 *
      ceiling(
        length(
          permutation_scan_results[[season_name]]$candidate_years
        ) / 3
      ),
    dpi = 300,
    limitsize = FALSE
  )
}

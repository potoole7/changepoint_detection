#### Simulation study for changepoint detection ####

# Simulate data (done)
# - Why is each date repeated 720 times, rather than 40, as in application? (done, problem in how dates were setup)
# Marginal transformation (done)
# TODO - Explain tanh linear predictor function
# TODO - Maybe move from using delta_z as argument to deriving from target rho
# Perform screening (done)
# TODO Change plotting code for screening
# TODO Perform changepoint detection algorithm for 1 iteration
# TODO Change plotting code for changepoints
# TODO Perform both for "many" iterations, with a given setup

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

#### metadata ####

# variables
# dep_var <- c("drought_local_rev")
dep_var <- c("drought_local")
# temp_var <- "temp_max"
temp_var <- "temp"

# number of locations, seasons and years to simulate for
n_locs <- 40
years <- 1960:2020
# seasons <- c("Winter", "Spring", "Summer", "Autumn")

seed <- 123 # random seed
# Conditional threshold and number of samples for Laplace sample used throughout
dqu <- 0.8
# dqu <- 0.85
n_samples <- 500

# run initially for just 100 permutations across full range
n_perm_screen <- 200L
n_years_per_block <- 25L # TODO Check this choice? Or just leave as best for app

# set minimum number of exceedances required for a successful fit
min_exceedances <- 15

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


#### Simulate Data ####

# Generate one simulated season with time-varying, site-specific
# bivariate t-copula dependence.
simulate_t_copula_season <- \(
  n_sites = 40L,
  n_years = 60L,
  n_obs_per_year = 12L,
  season_name = "Winter",
  start_year = 1961L,
  baseline_rho = c(
    rep(0.1, 14L),
    rep(0.5, 13L),
    rep(0.8, 13L)
  ),
  change_type = c("local", "global", "none"),
  affected_sites = 1L,
  change_start = 30L,
  change_end = 60L,
  delta_z = 0.5,
  df_t = 3,
  gpd_xi = c(-0.05, -0.05),
  gpd_sigma = c(1, 1),
  return_laplace = FALSE,
  seed = NULL
) {
  change_type <- match.arg(change_type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  stopifnot(
    n_sites >= 1L,
    n_years >= 1L,
    n_obs_per_year >= 1L,
    length(baseline_rho) == n_sites,
    all(abs(baseline_rho) < 1),
    length(gpd_xi) == 2L,
    length(gpd_sigma) == 2L,
    all(gpd_sigma > 0),
    change_start >= 1L,
    change_start <= n_years,
    change_end >= change_start,
    change_end <= n_years
  )

  site_ids <- sprintf(
    "site_%02d",
    seq_len(n_sites)
  )

  if (change_type == "none") {
    affected_index <- integer()
  } else if (change_type == "global") {
    affected_index <- seq_len(n_sites)
  } else {
    if (is.character(affected_sites)) {
      affected_index <- match(
        affected_sites,
        site_ids
      )

      if (anyNA(affected_index)) {
        stop(
          "At least one `affected_sites` name is invalid."
        )
      }
    } else {
      affected_index <- as.integer(
        affected_sites
      )

      if (
        anyNA(affected_index) ||
          any(!affected_index %in% seq_len(n_sites))
      ) {
        stop(
          "`affected_sites` contains invalid site indices."
        )
      }
    }

    affected_index <- unique(affected_index)
  }

  # Smooth ramp:
  #   0 before change_start
  #   linearly increases from 0 to 1
  #   1 from change_end onward
  year_index <- seq_len(n_years)

  change_progress <- pmin(
    1,
    pmax(
      0,
      (
        year_index - change_start
      ) /
        max(1, change_end - change_start)
    )
  )

  # Make the first affected year have a positive change when
  # change_start == change_end.
  if (change_start == change_end) {
    change_progress <- as.numeric(
      year_index >= change_start
    )
  }

  design <- tidyr::expand_grid(
    site_index = seq_len(n_sites),
    year_index = year_index
  ) |>
    dplyr::mutate(
      name = site_ids[site_index],
      season = season_name,
      season_year =
        start_year + year_index - 1L,
      baseline_rho =
        baseline_rho[site_index],
      cluster = dplyr::case_when(
        baseline_rho <= 0.1 ~ "low",
        baseline_rho <= 0.5 ~ "medium",
        TRUE ~ "high"
      ),
      affected =
        site_index %in% affected_index,
      change_progress =
        change_progress[year_index],
      eta = atanh(baseline_rho) +
        delta_z *
          change_progress *
          affected,
      rho = tanh(eta)
    )

  generated <- lapply(
    seq_len(nrow(design)),
    \(i) {
      cop <- copula::tCopula(
        param = design$rho[[i]],
        dim = 2L,
        df = df_t,
        df.fixed = TRUE,
        dispstr = "ex"
      )

      u <- copula::rCopula(
        n_obs_per_year,
        cop
      )

      tibble::tibble(
        within_year_index =
          seq_len(n_obs_per_year),
        U1 = u[, 1],
        U2 = u[, 2]
      )
    }
  )

  simulated <- dplyr::bind_cols(
    design,
    tibble::tibble(draws = generated)
  ) |>
    tidyr::unnest(draws) |>
    dplyr::mutate(
      X1 = qgpd(
        U1,
        xi = gpd_xi[[1]],
        sigma = gpd_sigma[[1]],
        u = 0
      ),
      X2 = qgpd(
        U2,
        xi = gpd_xi[[2]],
        sigma = gpd_sigma[[2]],
        u = 0
      )
    )

  if (return_laplace) {
    # Because the true marginal CDF values are known in simulation,
    # transform U directly to standard Laplace margins.
    simulated <- simulated |>
      dplyr::mutate(
        X1_laplace = qlaplace(U1),
        X2_laplace = qlaplace(U2)
      )
  }

  simulated |>
    dplyr::select(
      name,
      site_index,
      cluster,
      season,
      season_year,
      year_index,
      within_year_index,
      affected,
      change_progress,
      baseline_rho,
      rho,
      X1,
      X2,
      dplyr::any_of(
        c("X1_laplace", "X2_laplace")
      )
    )
}

sim_local <- simulate_t_copula_season(
  n_sites = 40L,
  n_years = 60L,
  baseline_rho = c(
    rep(0.2, 14L),
    rep(0.5, 13L),
    rep(0.8, 13L)
  ),
  # change_type = "global",
  change_type = "local",
  # affected_sites = 1:5,
  # affected_sites = 1:30,
  # affected_sites = 1:40,
  affected_sites = 1:14,
  change_start = 30L, # starts on middle year
  change_end = 60L,
  # TODO maybe parametrise by desired rho_t??
  delta_z = 0.3, # controls strength of increase
  # delta_z = 0.7, # controls strength of increase
  return_laplace = TRUE,
  seed = 123L
) |>
  group_by(name) |>
  dplyr::mutate(
    date = as.Date(
      paste0(season_year, "-01-01")
    ) +
      7L * (within_year_index - 1L)
  ) |>
  ungroup()

sim_local |>
  filter(name == "site_01") |>
  distinct(rho, baseline_rho, season_year)

# check everything's working properly
sim_local |>
  distinct(
    name,
    season_year,
    baseline_rho,
    rho,
    affected
  ) |>
  filter(
    name %in% c(
      "site_01",
      "site_15",
      "site_35"
    )
  ) |>
  ggplot(
    aes(
      x = season_year,
      y = rho,
      colour = name
    )
  ) +
  geom_line(linewidth = 1) +
  geom_hline(
    aes(yintercept = baseline_rho),
    linetype = "dashed",
    alpha = 0.4
  ) +
  labs(
    x = "Season-year",
    y = expression(rho[t]),
    colour = "Site"
  ) +
  theme_bw()


#### Marginal transformation ####

# convert to cecl_marg object
data_marg <- sim_local |>
  mutate(name = factor(name)) |>
  group_split(name) |>
  lapply(\(x) {
    ret1 <- x |>
      arrange(date) |>
      select(X1_laplace, X2_laplace) |>
      as.matrix()

    # colnames(ret1) <- c("X1", "X2")
    colnames(ret1) <- c(temp_var, dep_var)
    ret1
  })
names(data_marg) <- unique(sim_local$name)

marg <- as_cecl_marg(data_marg)

# add other metadata to the marg object
sites <- unique(sim_local$name)
# marg$dates <- unique(sim_local$date) # TODO Why do dates not go up to end of 2020?
marg$dates <- lapply(sites, \(x) as.character(unique(sim_local$date)))
names(marg$dates) <- sites

marg <- list("Winter" = marg) # make dummy list for different "seasons"


#### Screening ####

screen_setup_df <- tidyr::crossing(
  # "season" = seasons,
  "season"            = "Winter",
  "n_years_per_block" = c(15L, 20L, 25L, 30L)
)

# pull Laplace transformed data from marginal object, and join together
data_laplace_season <- lapply(marg, \(x) {
  bind_rows(Map(
    \(transformed, dates, station) {
      # print(nrow(transformed))
      # print(length(dates))
      # TODO Why is this not the case??
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
  setNames(names(marg))


# sink(paste0("sink_screening_output_dqu_", dqu, ".txt"))

screen_res_df <- bind_rows(lapply(seq_len(nrow(screen_setup_df)), \(i) {
  print(paste0(round(i / nrow(screen_setup_df) * 100, 2), "% of setups done"))
  print(paste0("season = ", screen_setup_df$season[[i]]))
  print(paste0("n_years_per_block = ", screen_setup_df$n_years_per_block[[i]]))
  with(
    screen_setup_df,
    screen_one_setting(
      data_laplace_season,
      season_name = "Winter",
      n_years_per_block[[i]],
      min_exceedances = min_exceedances
    )
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

# sink()

# save
# readr::write_csv(screen_res_df, "data/02_app/screen_res.csv.gz")
readr::write_csv(screen_res_df, paste0("data/01_sim/screen_res_dqu_", dqu, ".csv.gz"))

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

# plot_frob <- \(df, spec_year = NULL) {
#   df_plot <- df
#   if (!is.null(spec_year)) {
#     df_plot <- df_plot |>
#       filter(n_years_per_block == spec_year)
#   }
#   df_plot |>
#     ggplot(
#       aes(
#         x = change_after_year,
#         y = frob
#       )
#     ) +
#     geom_line(
#       data = \(x) filter(x, success),
#       aes(
#         group = interaction(
#           setting,
#           success_run
#         )
#       ),
#       colour = "black"
#     ) +
#     geom_point(
#       data = \(x) filter(x, success),
#       colour = "black"
#     ) +
#     geom_point(
#       data = \(x) filter(x, local_peak_frob),
#       colour = "red",
#       size = 3
#     ) +
#     geom_rug(
#       data = \(x) filter(x, !success),
#       aes(x = change_after_year),
#       inherit.aes = FALSE,
#       sides = "b",
#       colour = "grey50"
#     ) +
#     facet_wrap(
#       ~setting,
#       scales = "free_y"
#     ) +
#     scale_x_continuous(
#       breaks = year_breaks
#     ) +
#     labs(
#       # x = "Candidate change after seasonal year",
#       x = "season year",
#       y = "Frobenius discrepancy",
#       # caption = paste(
#       #   "Red points indicate local peaks.",
#       #   "Grey axis marks indicate failed candidates."
#       # )
#     ) +
#     cecl_theme() +
#     theme(
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1
#       )
#     )
# }

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
# plot_all_norms <- \(df, spec_year = NULL) {
#   df_plot <- df
#   if (!is.null(spec_year)) {
#     df_plot <- df_plot |>
#       filter(n_years_per_block == spec_year)
#   }
#   df_plot |>
#     # filter(n_years_per_block == 25) |>
#     ggplot(
#       aes(
#         x = change_after_year,
#         y = value,
#         colour = norm
#       )
#     ) +
#     geom_line(
#       data = \(x) filter(x, success),
#       aes(
#         group = interaction(
#           setting,
#           norm,
#           success_run
#         )
#       )
#     ) +
#     geom_point(
#       data = \(x) filter(x, success)
#     ) +
#     geom_point(
#       data = \(x) filter(x, local_peak),
#       colour = "red",
#       size = 3
#     ) +
#     geom_rug(
#       data = \(x) filter(x, !success),
#       aes(x = change_after_year),
#       inherit.aes = FALSE,
#       sides = "b",
#       colour = "grey50"
#     ) +
#     facet_wrap(
#       ~setting,
#       scales = "free_y"
#     ) +
#     scale_x_continuous(
#       breaks = year_breaks
#     ) +
#     scale_colour_brewer(
#       palette = "Dark2"
#     ) +
#     labs(
#       # x = "Candidate change after seasonal year",
#       x = "season year",
#       y = "Discrepancy",
#       colour = "Norm",
#       # caption = paste(
#       #   "Red points indicate norm-specific local peaks.",
#       #   "Grey axis marks indicate failed candidates."
#       # )
#     ) +
#     cecl_theme() +
#     theme(
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1
#       ),
#       legend.position = "bottom"
#     )
# }
plot_all_norms <- \(df, spec_year = NULL) {
  df_plot <- df
  if (!is.null(spec_year)) {
    df_plot <- df_plot |>
      filter(n_years_per_block == spec_year)
  }
  p <- df_plot |>
    group_by(norm, n_years_per_block) |>
    # scale between 0 and 1
    # mutate(value = scale(value, center = TRUE, scale = TRUE)) |>
    mutate(value = boot::inv.logit(scale(value))) |>
    filter(norm != "Spectral") |>
    ggplot(
      aes(
        x = change_after_year,
        y = value,
        # colour = norm
        colour = season
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
      ),
      show.legend = FALSE
    ) +
    geom_point(
      data = \(x) filter(x, success)
    ) +
    geom_point(
      data = \(x) filter(x, local_peak),
      # colour = "red",
      colour = "black",
      shape = 4,
      size = 5
    ) +
    # facet_wrap(
    #   ~setting,
    #   scales = "free_y"
    # ) +
    facet_wrap(
      ~norm,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = year_breaks
    ) +
    scale_colour_brewer(
      palette = "Dark2"
    ) +
    labs(
      x = "Season Year",
      y = expression(D),
      # colour = "Norm",
      colour = "Season",
    ) +
    cecl_theme() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 6)))
}


p_all_norms <- plot_all_norms(screen_res_long)
p_all_norms_25 <- plot_all_norms(screen_res_long, spec_year = 25)

# p_frob
# p_frob_25
#
# p_all_norms
p_all_norms_25


#### Changepoint Detection ####

seasons <- "Winter"
permutation_scan_results <- setNames(
  lapply(
    seq_along(seasons),
    \(i) {
      message("Season = ", seasons[[i]])
      # message("Season = Spring")
      run_season_permutation_scan(
        data_laplace_season,
        season_name = seasons[[i]],
        n_years_per_block = n_years_per_block,
        n_perm = n_perm_screen,
        min_exceedances = min_exceedances,
        seed = seed + i - 1L,
        use_start = TRUE,
        ret_dep = TRUE,
        verbose = TRUE,
        permutation_validation_warnings = TRUE
      )
    }
  ),
  seasons
)

saveRDS(
  permutation_scan_results,
  paste0(
    "data/01_sim/",
    "sim_perm_test_",
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
  ) |>
  print(n = Inf)

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

# Only one season!
permutation_scan_plots <- permutation_scan_plots[[1]]

# save
plot_directory <- paste0(
  "plots/02_sim/permutation_scan_n_perm_",
  n_perm_screen,
  "_dqu_",
  dqu
)

dir.create(
  plot_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# for (season_name in names(permutation_scan_plots)) {
# season_plots <- permutation_scan_plots[[season_name]]
season_plots <- permutation_scan_plots

filename_prefix <- paste0(
  # tolower(season_name),
  # "_",
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
  plot = season_plots$permutation_histograms,
  width = 12,
  height = 4 *
    ceiling(
      length(
        # permutation_scan_results[[season_name]]$candidate_years
        permutation_scan_results[[1]]$candidate_years
      ) / 3
    ),
  dpi = 300,
  limitsize = FALSE
)
# }

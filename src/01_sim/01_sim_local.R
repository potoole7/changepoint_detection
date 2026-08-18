#### Simulation study for changepoint detection ####

# Simulate data (done)
# TODO Marginal transformation
# TODO - Explain tanh linear predictor function
# TODO Perform changepoint detection algorithm for 1 iteration
# TODO Perform for "many" iterations

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
dep_var <- c("drought_local_rev")
temp_var <- "temp_max"

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
  change_start = 31L,
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
    rep(0.1, 14L),
    rep(0.5, 13L),
    rep(0.8, 13L)
  ),
  change_type = "local",
  # affected_sites = 1:5,
  # affected_sites = 1:30,
  affected_sites = 1:40,
  change_start = 31L,
  change_end = 60L,
  delta_z = 0.5, # controls strength of increase
  return_laplace = TRUE,
  seed = 123L
) |>
  dplyr::mutate(
    date = as.Date(
      paste0(season_year, "-01-01")
    ) +
      7L * (within_year_index - 1L)
  )

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

#### Screening ####


#### Changepoint Detection ####


#### Plotting ####

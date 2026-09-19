#### Plot results of simulation study ####

# TODO Change labels for local scenarios; seem a bit unprofessional

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

global_dir <- "data/01_sim/changepoints/global_first/"
# global_dir <- "data/01_sim/changepoints/global_second/"
# global_dir <- "data/01_sim/changepoints/global_third/"

# number of times to repeat simulations
nreps <- 100

# variables
dep_var <- c("drought_local")
temp_var <- "temp"

# number of locations, seasons and years to simulate for
n_locs <- 40
years <- 1960:2020
# seasons <- c("Winter", "Spring", "Summer", "Autumn")

seed <- 123 # random seed
# Conditional threshold and number of samples for Laplace sample used throughout
dqu <- 0.8

# run initially for just 100 permutations across full range
n_perm_screen <- 100L
n_years_per_block <- 25L

# Year when gradual change begins
cp_year <- 1990

#### Functions ####

# function to add Wilson confidence intervals to a data frame for a given rate and sample size
add_wilson_interval <- \(data, rate, n) {
  z <- qnorm(0.975)
  z2 <- z^2

  data |>
    mutate(
      .rate = {{ rate }},
      .n = {{ n }},
      ci_lower = (
        .rate +
          z2 / (2 * .n) -
          z * sqrt(
            .rate * (1 - .rate) / .n +
              z2 / (4 * .n^2)
          )
      ) / (
        1 + z2 / .n
      ),
      ci_upper = (
        .rate +
          z2 / (2 * .n) +
          z * sqrt(
            .rate * (1 - .rate) / .n +
              z2 / (4 * .n^2)
          )
      ) / (
        1 + z2 / .n
      )
    ) |>
    select(
      -.rate,
      -.n
    )
}

#### Screening ####

# # load data
# types <- c("none", "global", "local")
# screen_res_df <- bind_rows(
#   lapply(
#     types,
#     \(x) {
#       list.files(
#         "data/01_sim/screening/",
#         pattern = sprintf("%s_\\d{3}\\.csv.gz$", x),
#         full.names = TRUE
#       ) |>
#         lapply(read.csv) |>
#         bind_rows(.id = "simulation_id") |>
#         mutate(
#           simulation_id = as.integer(simulation_id),
#           scenario = x
#         )
#     }
#   )
# ) |>
#   mutate(
#     scenario = factor(
#       scenario,
#       levels = c("none", "local", "global"),
#       labels = c("No change", "Local change", "Global change")
#     )
#   )

screen_sources <- tibble::tribble(
  ~scenario, ~directory, ~file_prefix,
  "No change", "data/01_sim/screening", "none",
  "Local (medium -> high)", "data/01_sim/screening/local_med_high", "local",
  "Local (low -> high)", "data/01_sim/screening/local_low_high", "local",
  "Global change", "data/01_sim/screening/global_first", "global"
)

screen_res_df <- purrr::pmap_dfr(
  screen_sources,
  function(scenario, directory, file_prefix) {
    files <- list.files(
      directory,
      pattern = sprintf("^screen_%s_[0-9]{3}\\.csv\\.gz$", file_prefix),
      full.names = TRUE
    )

    if (length(files) == 0L) {
      stop("No screening files found in: ", directory)
    }

    names(files) <- sub(
      ".*_([0-9]{3})\\.csv\\.gz$",
      "\\1",
      basename(files)
    )

    bind_rows(lapply(files, read.csv), .id = "simulation_id") |>
      mutate(
        simulation_id = as.integer(simulation_id),
        scenario = scenario
      )
  }
) |>
  mutate(
    scenario = factor(scenario, levels = screen_sources$scenario)
  )

# screen_res_long <- screen_res_df |>
#   pivot_longer(
#     cols = c(frob, inf, inf2),
#     names_to = "norm",
#     values_to = "statistic"
#   ) |>
#   mutate(
#     norm = recode(
#       norm,
#       frob = "Frobenius",
#       inf  = "Maximum",  # verify these two mappings
#       inf2 = "Infinity"
#     )
#   )

screen_res_long <- screen_res_df |>
  pivot_longer(
    cols = c(frob, inf, inf2),
    names_to = "norm",
    values_to = "statistic"
  ) |>
  mutate(
    norm = recode(
      norm,
      frob = "Frobenius",
      inf  = "Maximum", # verify these two mappings
      inf2 = "Infinity"
    )
  )

# TODO Plot
alpha <- 0.05

null_critical_values <- screen_res_long |>
  filter(scenario == "No change") |>
  group_by(n_years_per_block, norm) |>
  summarise(
    critical_value = quantile(
      statistic,
      probs = 1 - alpha,
      na.rm = TRUE,
      type = 8
    ),
    .groups = "drop"
  )

screen_performance <- screen_res_long |>
  left_join(
    null_critical_values,
    by = c("n_years_per_block", "norm")
  ) |>
  group_by(
    scenario,
    n_years_per_block,
    norm
  ) |>
  summarise(
    n = sum(!is.na(statistic)),
    rejections = sum(
      statistic > critical_value,
      na.rm = TRUE
    ),
    rejection_rate = rejections / n,
    .groups = "drop"
  ) |>
  add_wilson_interval(
    rate = rejection_rate,
    n = n
  )

ggplot(
  screen_performance,
  aes(
    x = n_years_per_block,
    y = rejection_rate,
    colour = scenario,
    group = scenario
  )
) +
  geom_hline(
    yintercept = alpha,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_ribbon(
    aes(
      ymin = ci_lower,
      ymax = ci_upper,
      fill = scenario
    ),
    colour = NA,
    alpha = 0.12,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~norm) +
  scale_x_continuous(
    breaks = sort(
      unique(screen_performance$n_years_per_block)
    )
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Number of years per block",
    y = "Screening rejection rate",
    colour = "Scenario"
  ) +
  cecl_theme()

null_reference <- screen_res_long |>
  filter(scenario == "No change") |>
  group_by(
    n_years_per_block,
    # season_year,
    change_after_year,
    norm
  ) |>
  summarise(
    null_mean = mean(statistic, na.rm = TRUE),
    null_sd = sd(statistic, na.rm = TRUE),
    .groups = "drop"
  )

screen_standardised <- screen_res_long |>
  left_join(null_reference) |>
  mutate(
    z_statistic = (statistic - null_mean) / null_sd
  )

screen_profile <- screen_standardised |>
  # filter(n_years_per_block == 25L) |>
  filter(
    scenario != "No change"
  ) |>
  group_by(
    scenario,
    norm,
    # season_year
    change_after_year,
    n_years_per_block
  ) |>
  summarise(
    mean_z = mean(z_statistic, na.rm = TRUE),
    se_z = sd(z_statistic, na.rm = TRUE) /
      sqrt(sum(!is.na(z_statistic))),
    lower = mean_z - 1.96 * se_z,
    upper = mean_z + 1.96 * se_z,
    .groups = "drop"
  )

# TODO Give each plot it's own facets
p_screen <- screen_profile |>
  filter(norm != "Maximum") |>
  mutate(ind = paste(norm, " - ", scenario)) |>
  ggplot(
    aes(
      x = change_after_year,
      y = mean_z,
      colour = factor(n_years_per_block),
      group = factor(n_years_per_block)
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_vline(
    xintercept = 1990,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = factor(n_years_per_block)
    ),
    alpha = 0.2
  ) +
  geom_line(
    linewidth = 0.9, show.legend = FALSE
  ) +
  # facet_grid(
  #   norm ~ scenario,
  #   scales = "free_y"
  # ) +
  facet_wrap(
    # vars(norm, scenario),
    ~ ind,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    x    = "Candidate season-year",
    y    = "Null-standardised screening statistic",
    fill = "# years per block"
  ) +
  cecl_theme() +
  scale_x_continuous(limits = c(1970, 2010)) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1)
    )
  )

p_screen

ggsave(
  # "plots/01_sim/sim_screening_plot.png",
  "latex/plots/sim_screening_plot.png",
  plot = p_screen,
  width = 10,
  height = 8
)


# screen_peaks <- screen_standardised |>
#   filter(
#     #   # scenario != "No change",
#     n_years_per_block == 25L
#   ) |>
#   group_by(
#     simulation_id,
#     scenario,
#     n_years_per_block,
#     norm
#   ) |>
#   slice_max(
#     z_statistic,
#     n = 1,
#     with_ties = FALSE
#   ) |>
#   ungroup() |>
#   mutate(
#     absolute_error = abs(change_after_year - 1990),
#     inside_transition = between(
#       change_after_year,
#       1980,
#       2000
#     )
#   )
#
# screen_peaks |>
#   # filter(
#   #   scenario != "No change"
#   # ) |>
#   ggplot(
#     aes(x = change_after_year)
#   ) +
#   # annotate(
#   #   "rect",
#   #   xmin = 1980,
#   #   xmax = 2000,
#   #   ymin = -Inf,
#   #   ymax = Inf,
#   #   fill = "grey70",
#   #   alpha = 0.25
#   # ) +
#   geom_vline(
#     xintercept = 1990,
#     linetype = "dotted"
#   ) +
#   geom_histogram(
#     binwidth = 1,
#     boundary = 0.5,
#     fill = "#0072B2",
#     colour = "white"
#   ) +
#   facet_grid(norm ~ scenario) +
#   labs(
#     x = "Year of maximum screening statistic",
#     y = "Number of simulations"
#   ) +
#   theme_bw()
#
# peak_prob <- screen_peaks |>
#   count(scenario, norm, change_after_year, name = "n") |>
#   group_by(scenario, norm) |>
#   mutate(probability = n / sum(n)) |>
#   ungroup() |>
#   complete(
#     scenario,
#     norm,
#     change_after_year,
#     fill = list(
#       n = 0,
#       probability = 0
#     )
#   )
#
# null_peak_prob <- peak_prob |>
#   filter(scenario == "No change") |>
#   select(
#     norm,
#     change_after_year,
#     null_probability = probability
#   )
#
# peak_excess <- peak_prob |>
#   filter(scenario != "No change") |>
#   left_join(
#     null_peak_prob,
#     by = c("norm", "change_after_year")
#   ) |>
#   mutate(
#     excess_probability =
#       probability - null_probability
#   )
#
# ggplot(
#   peak_excess,
#   aes(
#     x = change_after_year,
#     y = excess_probability
#   )
# ) +
#   # annotate(
#   #   "rect",
#   #   xmin = 1980,
#   #   xmax = 2000,
#   #   ymin = -Inf,
#   #   ymax = Inf,
#   #   fill = "grey70",
#   #   alpha = 0.2
#   # ) +
#   geom_hline(
#     yintercept = 0,
#     colour = "grey40"
#   ) +
#   geom_vline(
#     xintercept = 1990,
#     linetype = "dotted"
#   ) +
#   geom_col(fill = "#0072B2") +
#   facet_grid(norm ~ scenario) +
#   scale_y_continuous(
#     labels = scales::label_percent()
#   ) +
#   labs(
#     x = "Year of maximum screening statistic",
#     y = "Excess peak probability relative to no change"
#   ) +
#   theme_bw()
#
#
# ## compare n_years_per_block
#
# screen_profile_blocks <- screen_standardised |>
#   group_by(
#     scenario,
#     norm,
#     n_years_per_block,
#     change_after_year
#   ) |>
#   summarise(
#     mean_z = mean(z_statistic, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   mutate(
#     n_years_per_block = factor(
#       n_years_per_block,
#       levels = c(15L, 20L, 25L, 30L),
#       labels = c(
#         "15 years",
#         "20 years",
#         "25 years",
#         "30 years"
#       )
#     )
#   )
#
# ggplot(
#   filter(
#     screen_profile_blocks,
#     scenario != "No change"
#   ),
#   aes(
#     x = change_after_year,
#     y = mean_z,
#     colour = n_years_per_block,
#     group = n_years_per_block
#   )
# ) +
#   annotate(
#     "rect",
#     xmin = 1980,
#     xmax = 2000,
#     ymin = -Inf,
#     ymax = Inf,
#     fill = "grey70",
#     alpha = 0.18
#   ) +
#   geom_hline(
#     yintercept = 0,
#     colour = "grey50",
#     linetype = "dashed"
#   ) +
#   geom_vline(
#     xintercept = 1990,
#     colour = "black",
#     linetype = "dotted"
#   ) +
#   geom_line(linewidth = 0.9) +
#   # facet_grid(
#   #   norm ~ scenario,
#   #   scales = "free_y"
#   # ) +
#   facet_wrap(
#     scenario ~ norm,
#     scales = "free_y"
#   ) +
#   labs(
#     x = "Candidate season-year",
#     y = "Mean null-standardised screening statistic"
#   ) +
#   # theme_bw() +
#   cecl_theme()


#### Changepoints: Plot Type 1 Error for null case ####

# load res
changepoint_res_none <- list.files(
  "data/01_sim/changepoints",
  pattern = sprintf("%s_\\d{3}\\.RDS$", "none"),
  full.names = TRUE
) |>
  lapply(readRDS)

# pull through summary dataframes
changepoint_df_none <- bind_rows(
  lapply(
    changepoint_res_none,
    `[[`,
    "summary"
  ),
  .id = "simulation_id"
) |>
  mutate(
    # label each simulation with an integer id
    simulation_id = as.integer(
      simulation_id
    ),
    # change year to date for plotting
    change_after_year = as.Date(paste0(change_after_year, "-01-01")),
    norm = recode(
      norm,
      "Infinity"  = "Maximum",
      "Infinity2" = "Infinity"
    )
  )

# Calculate global Type I error rates
global_error_1_sim_df <- changepoint_df_none |>
  group_by(
    simulation_id,
    norm
  ) |>
  summarise(
    n_valid_candidates = sum(
      !is.na(p_value)
    ),
    min_p_value = if (
      n_valid_candidates > 0
    ) {
      min(p_value, na.rm = TRUE)
    } else {
      NA_real_
    },
    reject_any_05 = if (
      n_valid_candidates > 0
    ) {
      any(p_value < 0.05, na.rm = TRUE)
    } else {
      NA
    },
    reject_any_10 = if (
      n_valid_candidates > 0
    ) {
      any(p_value < 0.10, na.rm = TRUE)
    } else {
      NA
    },
    .groups = "drop"
  )

# Calculate familywise Type I error rates
global_error_1_df <- global_error_1_sim_df |>
  group_by(norm) |>
  summarise(
    familywise_type1_05 = mean(
      reject_any_05,
      na.rm = TRUE
    ),
    familywise_type1_10 = mean(
      reject_any_10,
      na.rm = TRUE
    ),
    n_simulations = sum(
      !is.na(reject_any_05)
    ),
    .groups = "drop"
  )

# calculate pointwise Type I error (rejection rate at each candidate year)
pointwise_error_1_df <- changepoint_df_none |>
  group_by(
    change_after_year,
    norm
  ) |>
  summarise(
    type1_05 = mean(
      p_value < 0.05,
      na.rm = TRUE
    ),
    type1_10 = mean(
      p_value < 0.10,
      na.rm = TRUE
    ),
    n_valid = sum(!is.na(p_value)),
    .groups = "drop"
  )
pointwise_error_1_df

# Plot pointwise Type I error rates with Wilson confidence intervals
pointwise_error_1_plt <- pointwise_error_1_df |>
  pivot_longer(
    cols = c(
      type1_05,
      type1_10
    ),
    names_to = "threshold",
    values_to = "type1_error"
  ) |>
  mutate(
    nominal_alpha = case_when(
      threshold == "type1_05" ~ 0.05,
      threshold == "type1_10" ~ 0.10
    ),
    # threshold = factor(
    #   threshold,
    #   levels = c(
    #     "type1_05",
    #     "type1_10"
    #   ),
    #   labels = c(
    #     expression(alpha == 0.05),
    #     expression(alpha == 0.10)
    #   )
    # )
    threshold = factor(
      threshold,
      levels = c(
        "type1_05",
        "type1_10"
      ),
      labels = c(
        "p < 0.05",
        "p < 0.10"
      )
    )
  ) |>
  add_wilson_interval(
    rate = type1_error,
    n = n_valid
  )

p_cp1 <- pointwise_error_1_plt |>
  filter(norm != "Spectral") |>
  ggplot(
    aes(
      x = change_after_year,
      y = type1_error,
      colour = threshold,
      group = threshold
    )
  ) +
  geom_hline(
    data = pointwise_error_1_plt |>
      distinct(
        threshold,
        nominal_alpha
      ),
    aes(
      yintercept = nominal_alpha,
      colour = threshold
    ),
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.15,
    alpha = 0.55
  ) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    vars(norm)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(
      0,
      NA
    )
  ) +
  labs(
    x = "Season year",
    y = "Pointwise Type I error",
    colour = "P-value threshold"
  ) +
  cecl_theme() +
  scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_cp1

ggsave(
  "plots/01_sim/changepoint_pointwise_type1_error.png",
  plot = p_cp1,
  width = 10,
  height = 8
)

global_error_1_plt <- global_error_1_df |>
  pivot_longer(
    cols = c(
      familywise_type1_05,
      familywise_type1_10
    ),
    names_to = "threshold",
    values_to = "familywise_type1"
  ) |>
  mutate(
    nominal_alpha = case_when(
      threshold == "familywise_type1_05" ~ 0.05,
      threshold == "familywise_type1_10" ~ 0.10
    ),
    threshold = factor(
      threshold,
      levels = c(
        "familywise_type1_05",
        "familywise_type1_10"
      ),
      labels = c(
        "p < 0.05",
        "p < 0.10"
      )
    )
  ) |>
  add_wilson_interval(
    rate = familywise_type1,
    n = n_simulations
  )

p_cp2 <- global_error_1_plt |>
  filter(norm != "Spectral") |>
  ggplot(
    aes(
      x = norm,
      y = familywise_type1,
      fill = norm
    )
  ) +
  geom_hline(
    data = global_error_1_plt |>
      distinct(
        threshold,
        nominal_alpha
      ),
    aes(
      yintercept = nominal_alpha
    ),
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_col(
    width = 0.65,
    alpha = 0.75
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.15
  ) +
  facet_wrap(
    vars(threshold),
    labeller = label_parsed
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(
      0,
      NA
    )
  ) +
  labs(
    x = "Norm",
    # y = paste(
    #   "Simulations with at least",
    #   "one rejection"
    # ),
    y = "Null simulations with at least one rejection",
    fill = "Norm"
  ) +
  cecl_theme() +
  theme(
    legend.position = "none"
  )

p_cp2

ggsave(
  # "plots/01_sim/changepoint_global_type1_error.png",
  "latex/plots/sim_changepoint_global_type1_error.png",
  plot = p_cp2,
  width = 10,
  height = 8
)


#### Type II Error/Power (global change) ####

# load res
changepoint_res_global <- list.files(
  # "data/01_sim/changepoints",
  global_dir,
  pattern = sprintf("%s_\\d{3}\\.RDS$", "global"),
  full.names = TRUE
) |>
  lapply(readRDS)

# pull through summary dataframes
changepoint_df_global <- bind_rows(
  lapply(
    changepoint_res_global,
    `[[`,
    "summary"
  ),
  .id = "simulation_id"
) |>
  mutate(
    # label each simulation with an integer id
    simulation_id = as.integer(
      simulation_id
    ),
    # change year to date for plotting
    change_after_year = as.Date(paste0(change_after_year, "-01-01")),
    norm = recode(
      norm,
      "Infinity"  = "Maximum",
      "Infinity2" = "Infinity"
    )
  )

# Combine null and global-change simulations
changepoint_df_global_all <- bind_rows(
  changepoint_df_none |>
    mutate(
      scenario = "No change"
    ),
  changepoint_df_global |>
    mutate(
      scenario = "Global change"
    )
) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "No change",
        "Global change"
      )
    )
  )

# Calculate power at each year
global_pointwise_power_df <- changepoint_df_global_all |>
  group_by(
    scenario,
    change_after_year,
    norm
  ) |>
  summarise(
    power_05 = mean(
      p_value < 0.05,
      na.rm = TRUE
    ),
    power_10 = mean(
      p_value < 0.10,
      na.rm = TRUE
    ),
    n_valid = sum(
      !is.na(p_value)
    ),
    .groups = "drop"
  ) |>
  mutate(
    type2_05 = 1 - power_05,
    type2_10 = 1 - power_10
  )

# plot power
pointwise_power_plt <- global_pointwise_power_df |>
  select(
    scenario,
    change_after_year,
    norm,
    n_valid,
    power_05,
    power_10
  ) |>
  pivot_longer(
    cols = c(
      power_05,
      power_10
    ),
    names_to = "threshold",
    values_to = "rejection_rate"
  ) |>
  mutate(
    nominal_alpha = case_when(
      threshold == "power_05" ~ 0.05,
      threshold == "power_10" ~ 0.10
    ),
    threshold = factor(
      threshold,
      levels = c(
        "power_05",
        "power_10"
      ),
      labels = c(
        "p < 0.05",
        "p < 0.10"
      )
    )
  ) |>
  add_wilson_interval(
    rate = rejection_rate,
    n = n_valid
  )

# Calculate maximum power for each norm and threshold combination
pointwise_power_plt |>
  group_by(norm, threshold) |>
  filter(rejection_rate == max(rejection_rate, na.rm = TRUE)) |>
  arrange(
    desc(threshold),
    rejection_rate
  )

p_global_power <- pointwise_power_plt |>
  filter(
    norm != "Spectral",
    norm != "Maximum"
  ) |>
  mutate(ind = paste(norm, " - ", scenario)) |>
  ggplot(
    aes(
      x = change_after_year,
      y = rejection_rate,
      colour = threshold,
      group = threshold
    )
  ) +
  geom_hline(
    # data = global_pointwise_power_df |>
    data = pointwise_power_plt |>
      distinct(
        threshold,
        nominal_alpha
      ),
    aes(
      yintercept = nominal_alpha,
      colour = threshold
    ),
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = as.Date(
      paste0(cp_year, "-01-01")
    ),
    linetype = "dotted",
    colour = "grey30",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 60,
    alpha = 0.45
  ) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.8
  ) +
  geom_point(
    size = 2
  ) +
  # facet_grid(
  #   rows = vars(norm),
  #   cols = vars(scenario)
  # ) +
  facet_wrap(~ ind, scales = "fixed") +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2000-01-01")),
    date_labels = "%Y",
    date_breaks = "2 years"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(0, NA),
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    x = "Season year",
    y = "Pointwise rejection rate",
    colour = "P-value threshold"
  ) +
  cecl_theme() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_global_power

ggsave(
  # "plots/01_sim/changepoint_global_change_pointwise_power.png",
  "latex/plots/sim_changepoint_global_change_pointwise_power.png",
  # "latex/plots/sim_changepoint_global_change_pointwise_power_curve.png",
  plot = p_global_power,
  width = 10,
  height = 8
)


#### Type I Error ####

# # load res
# local_vals <- c(0.3, 0.6)
# changepoint_df_local <- bind_rows(lapply(local_vals, \(x) {
#   changepoint_res_spec <- list.files(
#     paste0("data/01_sim/changepoints/local_", x),
#     full.names = TRUE
#   ) |>
#     lapply(readRDS)
#
#   bind_rows(
#     lapply(
#       changepoint_res_spec,
#       `[[`,
#       "summary"
#     ),
#     .id = "simulation_id"
#   ) |>
#     mutate(
#       simulation_id = as.integer(simulation_id),
#       change_after_year = as.Date(paste0(change_after_year, "-01-01"))
#     )
# }), .id = "cp_size") |>
#   mutate(
#     cp_size = case_when(
#       cp_size == 1 ~ local_vals[[1]],
#       cp_size == 2 ~ local_vals[[2]],
#       TRUE ~ NA
#     )
#   )

# load res
# changepoint_res_local <- list.files(
#   "data/01_sim/changepoints",
#   pattern = sprintf("%s_\\d{3}\\.RDS$", "local"),
#   full.names = TRUE
# ) |>
#   lapply(readRDS)

# load res
# changepoint_res_local <- list.files(
#   "data/01_sim/changepoints",
#   pattern = sprintf("%s_\\d{3}\\.RDS$", "local"),
#   full.names = TRUE
# ) |>
#   lapply(readRDS)

dirs <- list.dirs(
  "data/01_sim/changepoints",
  recursive = FALSE,
  full.names = TRUE
)
dirs <- dirs[grepl("local", dirs)]

changepoint_res_local <- unlist(lapply(dirs, \(dir) {
  # load changepoint results for specific run (and setting)
  out <- list.files(
    dir,
    # pattern = sprintf("%s_\\d{3}\\.RDS$", "local"),
    full.names = TRUE
  ) |>
    lapply(readRDS)
  # label each simulation with the scenario name (i.e. local change magnitude)
  lapply(out, \(x) {
    # x$scenario <- basename(dir)
    x$summary$scenario <- basename(dir)
    x
  })
}), recursive = FALSE)

# pull through summary dataframes
# changepoint_df_local <- bind_rows(
#   lapply(
#     changepoint_res_local,
#     `[[`,
#     "summary"
#   ),
#   .id = "simulation_id"
# )

n_sim <- length(changepoint_res_local) / length(dirs)
changepoint_df_local <- bind_rows(
  lapply(
    changepoint_res_local,
    `[[`,
    "summary"
  ),
  .id = "simulation_id"
) |>
  mutate(
    # label each simulation with an integer id
    simulation_id = as.integer(
      simulation_id
    ),
    # change year to date for plotting
    change_after_year = as.Date(paste0(change_after_year, "-01-01")),
    norm = recode(
      norm,
      "Infinity"  = "Maximum",
      "Infinity2" = "Infinity"
    )
  )

# correct simulation_id for multiple local change scenarios
if (length(dirs) > 1) {
  changepoint_df_local <- changepoint_df_local |>
    mutate(
      simulation_id = ifelse(
        simulation_id > n_sim,
        simulation_id - n_sim,
        simulation_id
      )
    )
} else {
  changepoint_df_local <- changepoint_df_local |>
    mutate(scenario = "Local change")
}

# also add null case to this
changepoint_df_local_all <- bind_rows(
  changepoint_df_none |>
    mutate(scenario = "No change"),
  changepoint_df_local
) |>
  # mutate(scenario = factor(scenario, levels = c("No change", "Local change")))
  mutate(
    scenario = case_when(
      scenario == "local_med_high" ~ "Local (medium -> high)",
      scenario == "local_low_high" ~ "Local (low -> high)",
      TRUE ~ scenario
    ),
    scenario = factor(
      scenario,
      levels = c("No change", "Local (medium -> high)", "Local (low -> high)")
    )
  )

# calculate power for each value
local_pointwise_power_df <- changepoint_df_local_all |>
  group_by(
    # cp_size,
    scenario,
    change_after_year,
    norm
  ) |>
  summarise(
    power_05 = mean(
      p_value < 0.05,
      na.rm = TRUE
    ),
    power_10 = mean(
      p_value < 0.10,
      na.rm = TRUE
    ),
    n_valid = sum(
      !is.na(p_value)
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(
      power_05,
      power_10
    ),
    names_to = "threshold",
    values_to = "rejection_rate"
  ) |>
  mutate(
    nominal_alpha = case_when(
      threshold == "power_05" ~ 0.05,
      threshold == "power_10" ~ 0.10
    ),
    threshold = factor(
      threshold,
      levels = c(
        "power_05",
        "power_10"
      ),
      labels = c(
        "p < 0.05",
        "p < 0.10"
      )
    ),
    # scenario = factor(
    #   cp_size,
    #   levels = c(0, local_vals),
    #   # TODO Change
    #   labels = c(
    #     "No change",
    #     # "Correlation increase = 0.3",
    #     "Local change (0.3)",
    #     "Local change (0.6)"
    #   )
    # )
  ) |>
  add_wilson_interval(
    rate = rejection_rate,
    n = n_valid
  )

p_local_power <- local_pointwise_power_df |>
  filter(
    # norm != "Sectral"
    norm != "Maximum"
  ) |>
  # TODO Make ind a factor!
  mutate(
    ind = paste(norm, "-", scenario),
    ind = factor(
      ind,
      levels = crossing(
        "norm" = c("Frobenius", "Infinity"),
        "scenario" = c("No change", "Local (medium -> high)", "Local (low -> high)")
      ) |>
        arrange(norm, desc(scenario)) |>
        mutate(ind = paste0(norm, " - ", scenario)) |>
        pull(ind)
    )
  ) |>
  ggplot(
    aes(
      x = change_after_year,
      y = rejection_rate,
      colour = threshold,
      group = threshold
    )
  ) +
  geom_hline(
    data = local_pointwise_power_df |>
      distinct(
        threshold,
        nominal_alpha
      ),
    aes(
      yintercept = nominal_alpha,
      colour = threshold
    ),
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = as.Date(
      paste0(cp_year, "-01-01")
    ),
    linetype = "dotted",
    colour = "grey30",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    # Width is measured in days because x is a Date
    width = 60,
    alpha = 0.45,
    position = position_dodge(
      width = 50
    )
  ) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.8
  ) +
  geom_point(
    size = 2,
    position = position_dodge(
      width = 50
    )
  ) +
  # facet_grid(
  #   rows = vars(norm),
  #   cols = vars(scenario)
  # ) +
  # facet_wrap(
  #   norm ~ scenario, scales = "free"
  # ) +
  facet_wrap(~ ind, scales = "fixed") +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2000-01-01")),
    date_labels = "%Y",
    date_breaks = "2 years"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(0, NA),
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    x = "Season year",
    y = "Pointwise rejection rate",
    colour = "P-value threshold"
  ) +
  cecl_theme() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_local_power

ggsave(
  # "plots/01_sim/changepoint_local_change_pointwise_power_curve.png",
  # "latex/plots/sim_changepoint_local_change_pointwise_power_curve.png",
  "latex/plots/sim_changepoint_local_change_pointwise_power.png",
  plot = p_local_power,
  width = 10,
  height = 8
)

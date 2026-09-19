#### Examine correlations in application data ####

# Struggling to find a good simulation setup which has decent power

# To aid in this, we want to see what correlations look like in the application
# data, and how they change over time (and before/after changepoints).
# This might help us design something similar
# Choosing Summer to look at as it had by far the lowest p-values for quite
# a number of years

# TODO Could also do decades, rather than yearly?? Would be less noisy
# TODO
# TODO Summarise

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


#### Metadata ####

dep_var <- c("drought_local_rev")
temp_var <- "temp_max"
decades <- seq(1960, 2010, by = 10)
seasons <- c("Winter", "Spring", "Summer", "Autumn")

# identified changepoint years
change_year_winter <- 1998
change_year_summer <- 1990
change_year_autumn <- 1986 # change for autumn at 10% level

change_year_df <- data.frame(
  season = c("Winter", "Summer", "Autumn"),
  change_year = c(change_year_winter, change_year_summer, change_year_autumn)
)

use_laplace <- FALSE
use_laplace <- TRUE # Use Laplace data instead of original scale


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

station_names <- unique(data$station_name)

# # marginal fits
marg_season <- readRDS("data/02_app/marg_season_roll_emp.rds")
marg_season_summer <- marg_season$Summer

# also load in clusterings before and after changepoints
clust_sol <- readRDS("data/02_app/clust_cp_plt_all_dqu_0.8.rds")
clust_sol_summer <- clust_sol[grepl("Summer", names(clust_sol))]

# to access clustering solution, do e.g.
clust_sol_summer$`Summer - 1991-2020`$pam$clustering

# pull through clustering solutions
clust_sol_df <- bind_rows(lapply(seq_along(clust_sol_summer), \(i) {
  clust <- clust_sol_summer[[i]]$pam$clustering

  ret <- data.frame(
    "name" = names(clust),
    "cluster" = unname(clust),
    "season" = "Summer",
    "period" = names(clust_sol_summer)[i]
  )
  ret$period <- stringr::str_remove(ret$period, "Summer - ")
  ret
}))

#### Pull through Laplace estimates ####

station_names_marg <- names(marg_season_summer$transformed)
data_laplace_summer <- bind_rows(lapply(station_names_marg, \(spec_loc) {
  trans <- data.frame(marg_season_summer$transformed[[spec_loc]]) |>
    mutate(
      station_name = spec_loc,
      date = marg_season_summer$dates[[spec_loc]]
    )
})) |>
  # TODO Pull year (as season_year)
  mutate(
    season_year = as.numeric(substr(date, 0, 4))
  ) |>
  relocate(station_name, date, season_year)


#### Summer Kendall dependence ####

summer_change_year <- change_year_summer

# These reproduce the two 25-season-year blocks used by the test
# summer_pre_years <- (summer_change_year - 24L):summer_change_year
summer_pre_years <-  1960:summer_change_year
# summer_post_years <- (summer_change_year + 1L):
# (summer_change_year + 25L)
summer_post_years <- (summer_change_year + 1):2024

# function to calculate Kendall's tau and the corresponding Pearson correlation coefficient
kendall_dependence <- \(x, y) {
  keep <- complete.cases(x, y) # ensure that we have no missing values

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 3L ||
      dplyr::n_distinct(x) < 2L ||
      dplyr::n_distinct(y) < 2L
  ) {
    return(
      tibble(
        n = length(x),
        tau = NA_real_,
        rho = NA_real_
      )
    )
  }

  tau <- cor(
    x,
    y,
    method = "kendall"
  )

  # convert from tau to eliptical copula correlation coefficient (rho)
  tibble(
    n = length(x),
    tau = tau,
    rho = sin(pi * tau / 2)
  )
}

summer_data <- data |>
  filter(season == "Summer")
# replace original scale data with Laplace scale, to be more consistent with
# simulations, which are done on Laplace scale
if (use_laplace == TRUE) {
  summer_data <- data_laplace_summer
}


#### Annual estimates ####

# First, calculate kendalls tau for every year
summer_dependence_yearly <- summer_data |>
  group_by(
    # station_id,
    station_name,
    season_year
  ) |>
  group_modify(
    \(dat, key) {
      kendall_dependence(
        dat[[temp_var]],
        dat[[dep_var]]
      )
    }
  ) |>
  ungroup()

summer_dependence_yearly


#### Pre- and post-changepoint estimates ####

# Next, calculate kendalls tau for the two 25-year blocks before and after
# the changepoint (hopefully less noisy)
summer_dependence_period <- summer_data |>
  mutate(
    period = case_when(
      season_year %in% summer_pre_years ~ "Before",
      season_year %in% summer_post_years ~ "After",
      TRUE ~ NA_character_
    ),
    period = factor(
      period,
      levels = c("Before", "After")
    )
  ) |>
  filter(!is.na(period)) |>
  group_by(
    # station_id,
    station_name,
    period
  ) |>
  group_modify(
    \(dat, key) {
      kendall_dependence(
        dat[[temp_var]],
        dat[[dep_var]]
      )
    }
  ) |>
  ungroup()

summer_dependence_period

# tabulate in wide form
summer_dependence_change <- summer_dependence_period |>
  select(
    # station_id,
    station_name,
    period,
    n,
    tau,
    rho
  ) |>
  pivot_wider(
    names_from = period,
    values_from = c(n, tau, rho),
    names_glue = "{.value}_{tolower(period)}"
  ) |>
  mutate(
    delta_tau = tau_after - tau_before,
    delta_rho = rho_after - rho_before,
    abs_delta_rho = abs(delta_rho)
  ) |>
  arrange(desc(abs_delta_rho))

summer_dependence_change


#### Add Summer cluster membership ####

# process cluster memberships for Summer before and after changepoints,
summer_clusters <- clust_sol_df |>
  filter(season == "Summer") |>
  mutate(
    cluster_period = case_when(
      period == "1960-1990" ~ "cluster_before",
      period == "1991-2020" ~ "cluster_after", # TODO Should be 2024 lol
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(cluster_period)) |>
  select(
    station_name = name,
    cluster_period,
    cluster
  ) |>
  pivot_wider(
    names_from = cluster_period,
    values_from = cluster
  ) |>
  # also check where sites have switched clusters
  mutate(
    cluster_before = factor(cluster_before),
    cluster_after = factor(cluster_after),
    switched = cluster_before != cluster_after,
    transition = paste0(
      cluster_before,
      " -> ",
      cluster_after
    )
  )

summer_clusters

# find original size of clusters
summer_clusters |>
  group_by(cluster_before) |>
  summarise(n(), .groups = "drop")
# 18, 8, 14

# find final size of clusters
summer_clusters |>
  group_by(cluster_after) |>
  summarise(n(), .groups = "drop")
# 18. 6, 16

# count the number of switching clusters for each cluster
summer_clusters |>
  group_by(cluster_before) |>
  summarise(sum(switched))
# 8, 5, 6

# join cluster membership information in
summer_dependence_change <- summer_dependence_change |>
  left_join(
    summer_clusters,
    by = "station_name"
  )

summer_dependence_yearly <- summer_dependence_yearly |>
  left_join(
    summer_clusters,
    by = "station_name"
  )


#### Plotting ####

# cluster_colours <- c(
#   "1" = "#0072B2",
#   "2" = "#D55E00",
#   "3" = "#009E73"
# )

# plot correlations before and after changepoint
p1 <- ggplot(
  summer_dependence_change,
  aes(
    x = rho_before,
    y = rho_after,
    colour = cluster_before, # colour by original cluster
    shape = switched # shape by whether the station switched clusters
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  cecl_theme(legend.position = "bottom") +
  scale_shape_manual(
    values = c(
      "FALSE" = 16,
      "TRUE" = 17
    ),
    labels = c(
      "FALSE" = "Did not switch",
      "TRUE" = "Switched"
    ),
    name = "Cluster switching"
  ) +
  coord_equal() +
  labs(
    x = "Correlation before changepoint",
    y = "Correlation after changepoint",
    colour = "Cluster before changepoint"
  )
p1

plt_name1 <- ifelse(isTRUE(use_laplace), "plot1_laplace.png", "plot1.png")
ggsave(plt_name1, plot = p1, width = 12, height = 10)

# also facet by cluster before and after changepoint
p1a <- ggplot(
  summer_dependence_change,
  aes(
    x = rho_before,
    y = rho_after,
    colour = cluster_before, # colour by original cluster
    shape = switched # shape by whether the station switched clusters
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  facet_grid(cluster_before ~ cluster_after) +
  cecl_theme(legend.position = "bottom") +
  scale_shape_manual(
    values = c(
      "FALSE" = 16,
      "TRUE" = 17
    ),
    labels = c(
      "FALSE" = "Did not switch",
      "TRUE" = "Switched"
    ),
    name = "Cluster switching"
  ) +
  coord_equal() +
  labs(
    x = "Correlation before changepoint",
    y = "Correlation after changepoint",
    colour = "Cluster before changepoint"
  )
p1a

plt_name1a <- ifelse(isTRUE(use_laplace), "plot1a_laplace.png", "plot1a.png")
ggsave(plt_name1a, plot = p1a, width = 12, height = 10)

# find the largest changes in correlation, to label on the plot
largest_changes <- summer_dependence_change |>
  slice_max(
    abs_delta_rho,
    n = 8L,
    with_ties = FALSE
  )

p2 <- p_before_after +
  geom_text(
    data = largest_changes,
    aes(label = station_name),
    size = 3,
    nudge_y = 0.015,
    show.legend = FALSE,
    check_overlap = TRUE
  )
p2

plt_name2 <- ifelse(isTRUE(use_laplace), "plot2_laplace.png", "plot2.png")
ggsave(plt_name2, plot = p2, width = 12, height = 10)

# also plot yearly
p3 <- summer_dependence_yearly |>
  arrange(cluster_before, station_name, season_year) |>
  mutate(
    station_name = factor(station_name, levels = unique(station_name))
  ) |>
  ggplot(
    aes(
      x = season_year,
      y = rho,
      # y = tau,
      colour = cluster_before,
      # group = station_name,
      linetype = switched
    )
  ) +
  # facet_wrap(~cluster_before) +
  facet_wrap(~station_name) +
  geom_line(
    alpha = 0.3
  ) +
  geom_smooth() +
  cecl_theme() +
  labs(
    x = "Year",
    y = "Correlation",
    colour = "Cluster before changepoint",
    linetype = "Cluster switching"
  ) +
  # remove grey backgrounds from legend
  guides(
    colour = guide_legend(override.aes = list(fill = NA)),
    linetype = guide_legend(override.aes = list(fill = NA, linewidth = 4, alpha = 1))
  )
p3

plt_name3 <- ifelse(isTRUE(use_laplace), "plot3_laplace.png", "plot3.png")
ggsave(plt_name3, plot = p3, width = 12, height = 10)

p4 <- summer_dependence_yearly |>
  arrange(cluster_before, station_name, season_year) |>
  mutate(
    station_name = factor(station_name, levels = unique(station_name))
  ) |>
  ggplot(
    aes(
      x = season_year,
      y = rho,
      # y = tau,
      colour = cluster_before,
      # group = station_name,
      # linetype = switched
    )
  ) +
  # facet_wrap(~cluster_before) +
  # facet_wrap(~cluster_before) +
  facet_grid(cluster_before ~ cluster_after) +
  geom_line(
    aes(group = station_name),
    alpha = 0.3
  ) +
  geom_smooth() +
  cecl_theme() +
  labs(
    x = "Year",
    y = "Correlation",
    colour = "Cluster before changepoint"
  ) +
  guides(
    colour = guide_legend(override.aes = list(fill = NA, linewidth = 3))
  )
p4

plt_name4 <- ifelse(isTRUE(use_laplace), "plot4_laplace.png", "plot4.png")
ggsave(plt_name4, plot = p4, width = 12, height = 10)


#### Cluster-correlation bars ####

cluster_rho_order_before <- summer_dependence_change |>
  group_by(cluster_before) |>
  summarise(
    n_sites = n(),
    median_rho = median(
      rho_before,
      na.rm = TRUE
    ),
    rho_q25 = quantile(rho_before, 0.25, na.rm = TRUE),
    rho_q75 = quantile(rho_before, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(median_rho) |>
  mutate(
    simulated_regime = c(
      "low",
      "medium",
      "high"
    )
  )

# ensure simulated regimes are mapped correctly afterwards
regime_mapping <- cluster_rho_order_before |>
  select(cluster_before, simulated_regime)

cluster_rho_order_after <- summer_dependence_change |>
  group_by(cluster_before) |>
  summarise(
    n_sites = n(),
    median_rho = median(rho_after, na.rm = TRUE),
    rho_q25 = quantile(rho_after, 0.25, na.rm = TRUE),
    rho_q75 = quantile(rho_after, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    regime_mapping,
    by = "cluster_before"
  )

# plot
cluster_rho_summary <- bind_rows(
  cluster_rho_order_before |>
    rename(cluster = cluster_before) |>
    mutate(period = "Before"),
  cluster_rho_order_after |>
    rename(cluster = cluster_before) |>
    mutate(period = "After")
) |>
  mutate(
    cluster = factor(
      cluster,
      levels = c("1", "2", "3")
    ),
    period = factor(
      period,
      levels = c("Before", "After")
    )
  ) |>
  arrange(period, cluster)

sink_name <- ifelse(isTRUE(use_laplace), "bar_laplace.txt", "bar.txt")
sink(sink_name)
cluster_rho_summary
sink()

p_bar <- ggplot(
  cluster_rho_summary,
  aes(
    y = cluster,
    x = median_rho,
    fill = cluster
  )
) +
  geom_col(
    width = 0.65,
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(
      xmin = rho_q25,
      xmax = rho_q75
    ),
    width = 0.2,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~period
  ) +
  labs(
    x = "Median implied copula correlation",
    y = "Cluster",
    fill = "Cluster"
  ) +
  cecl_theme(
    legend.position = "bottom"
  )
p_bar

plt_name_bar <- ifelse(isTRUE(use_laplace), "plot_bar_laplace.png", "plot_bar.png")
ggsave(plt_name_bar, plot = p_bar, width = 12, height = 10)

#### Cluster-switching correlation bars ####

# Do the same, but inspect how switching effects things
cluster_rho_order_before_switch <- summer_dependence_change |>
  group_by(cluster_before, cluster_after) |>
  summarise(
    n_sites = n(),
    median_rho = median(
      rho_before,
      na.rm = TRUE
    ),
    rho_q25 = quantile(rho_before, 0.25, na.rm = TRUE),
    rho_q75 = quantile(rho_before, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(median_rho) |>
  # mutate(
  #   simulated_regime = c(
  #     "low",
  #     "medium",
  #     "high"
  #   )
  # ) |>
  identity()

cluster_rho_order_after_switch <- summer_dependence_change |>
  group_by(cluster_before, cluster_after) |>
  summarise(
    n_sites = n(),
    median_rho = median(rho_after, na.rm = TRUE),
    rho_q25 = quantile(rho_after, 0.25, na.rm = TRUE),
    rho_q75 = quantile(rho_after, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# plot
cluster_rho_summary_switch <- bind_rows(
  cluster_rho_order_before_switch |>
    # rename(cluster = cluster_before) |>
    mutate(period = "Before"),
  cluster_rho_order_after_switch |>
    # rename(cluster = cluster_before) |>
    mutate(period = "After")
) |>
  mutate(
    # cluster = factor(
    #   cluster,
    #   levels = c("1", "2", "3")
    # ),
    # cluster_before = factor(
    #   cluster_before,
    #   levels = c("1", "2", "3")
    # ),
    # cluster_after = factor(
    #   cluster_before,
    #   levels = c("1", "2", "3")
    # ),
    period = factor(
      period,
      levels = c("Before", "After")
    )
  ) |>
  arrange(period, cluster_before, cluster_after)

sink_name_switch <- ifelse(isTRUE(use_laplace), "bar_switch_laplace.txt", "bar_switch.txt")
sink(sink_name_switch)
cluster_rho_summary_switch
sink()

p_bar_switch <- ggplot(
  cluster_rho_summary_switch,
  aes(
    # y = cluster,
    y = period,
    x = median_rho,
    # fill = cluster
    fill = period
  )
) +
  geom_col(
    width = 0.65,
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(
      xmin = rho_q25,
      xmax = rho_q75
    ),
    width = 0.2,
    linewidth = 0.7
  ) +
  # facet_wrap(~period) +
  facet_grid(cluster_before~cluster_after) +
  labs(
    x = "Median implied copula correlation",
    y = "Cluster",
    fill = "Cluster"
  ) +
  cecl_theme(
    legend.position = "bottom"
  )
p_bar_switch

plt_name_bar_switch <- ifelse(isTRUE(use_laplace), "plot_bar_switch_laplace.png", "plot_bar_switch.png")
ggsave(plt_name_bar_switch, plot = p_bar_switch, width = 12, height = 10)

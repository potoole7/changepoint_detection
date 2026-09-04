#### Simulation study for changepoint detection, for 100 repititions ####

# TODO
# - Also run for 15 years per block (and not 30) for 15 years per block

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

# devtools::load_all("../CeCl")
library(CeCl)
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

# number of times to repeat simulations
# nreps <- 50 # start with 50, then move on to 100
nreps <- 100

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

mc_cores <- parallel::detectCores() - 1L

# run initially for just 100 permutations across full range
n_perm_screen <- 200L
# n_perm_screen <- 100L
n_years_per_block <- 25L # TODO Check this choice? Or just leave as best for app

# set minimum number of exceedances required for a successful fit
min_exceedances <- 15

cp_type <- "global"
# cp_type <- "local"
# cp_type <- "none"

baseline_rho <- final_rho <- c(
  low = 0.1,
  medium = 0.5,
  high = 0.8
)

if (cp_type == "global") {
  final_rho <- c(
    low = 0.4,
    medium = 0.8,
    high = 0.9
  )
} else if (cp_type == "local") {
  final_rho <- c(
    low = 0.5,
    medium = 0.5,
    high = 0.8
  )
}

# start and end years for the changepoint (for global change, all sites affected; for local change, only some sites affected)
start_year <- 1980L
end_year <- 2000L

save_dir <- "data/01_sim"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


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


#### Loop ####

screen_setup_df <- tidyr::crossing(
  "season"            = "Winter",
  # TODO Change back to include 15
  "n_years_per_block" = c(15L, 20L, 25L, 30L)
  # "n_years_per_block" = c(20L, 25L, 30L)
)

# arguments for simulate_t_copula_season function
sim_args <- list(
  # chosen to emulate the application data
  "n_sites" = 40L,
  "n_years" = 60L,
  # well separated "clusters" of sites, with different dependence structures
  # "baseline_rho"   = c(
  #   rep(0.1, 14L),
  #   rep(0.5, 13L),
  #   rep(0.8, 13L)
  # ),
  "baseline_rho" = baseline_rho,
  "final_rho" = final_rho,
  "change_type" = cp_type,
  # TODO maybe parametrise by desired rho_t??
  # "delta_z"        = 0.45,  # approx gives 0.5 for first cluster (atan(0.5) - atain(0.1) ~= 0.45)
  # "delta_z"        = 0.21,  # approx gives 0.3
  # "delta_z"        = 0.593, # approx gives 0.6
  # TODO Change to 25:35, see if that helps? More centred around 25 season years
  # "change_start"   = 30L, # starts on middle year
  # "change_end"     = 60L,
  "change_start_year" = start_year,
  "change_end_year" = end_year,
  "affected_sites" = ifelse(cp_type == "local", 1:5, NA_integer_),
  "return_laplace" = TRUE
)

res <- lapply(seq_len(nreps), \(k) {
  system(sprintf(
    'echo "\n%s\n"',
    paste0(k, " of ", nreps, " repititions")
  ))

  # skip iteration if you've already saved the associated file
  file <- sprintf(
    "data/01_sim/changepoints/changepoint_%s_%03d.RDS",
    cp_type,
    k
  )

  if (file.exists(file)) {
    return(NULL)
  }

  ## Simulate Data ##
  sim_local <- do.call(
    simulate_t_copula_season_line,
    c(sim_args, list("seed" = seed + k))
  ) |>
    group_by(name) |>
    dplyr::mutate(
      date = as.Date(
        paste0(season_year, "-01-01")
      ) +
        7L * (within_year_index - 1L)
    ) |>
    ungroup()

  ## Marginal transformation ##

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
  marg$dates <- lapply(sites, \(x) as.character(unique(sim_local$date)))
  names(marg$dates) <- sites

  marg <- list("Winter" = marg) # make dummy list for different "seasons"

  # pull Laplace transformed data from marginal object, and join together
  data_laplace_season <- lapply(marg, \(x) {
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


  ## Screening ##

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
        min_exceedances = min_exceedances,
        mc.cores = mc_cores
      )
    )
  }))

  ## Changepoint Detection ##

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
          seed = seed + k - 1L,
          use_start = TRUE,
          ret_dep = TRUE,
          verbose = TRUE,
          permutation_validation_warnings = TRUE,
          mc_cores = mc_cores
        )
      }
    ),
    seasons
  )


  # preprocess every season
  permutation_scan_tidy <- lapply(
    permutation_scan_results,
    preprocess_permutation_scan
  )[[1]] # only using one season!


  # outputs
  # return(list(
  #   "screen_res"      = screen_res_df,
  #   "changepoint_res" = permutation_scan_tidy
  # ))

  # save here so we don't have to keep in memory
  readr::write_csv(
    screen_res_df,
    sprintf(
      "data/01_sim/screening/screen_%s_%03d.csv.gz",
      cp_type,
      k
    )
  )

  saveRDS(object = permutation_scan_tidy, file)

  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(k / nreps, 3) * 100, "% completed", collapse = "")
  ))

  NULL
})


#### Join outputs ####

# TODO Change for changepoints since we output a list
screen_res_all <- list.files(
  "data/01_sim/screening",
  pattern = sprintf("%s_\\d{3}\\.csv\\.gz$", cp_type),
  full.names = TRUE
) |>
  lapply(readr::read_csv, show_col_types = FALSE) |>
  bind_rows()

readr::write_csv(
  screen_res_all,
  paste0("data/01_sim/screen_res_all_", cp_type, ".csv.gz")
)

# changepoint_res_all <- list.files(
#   "data/01_sim/changepoints",
#   pattern = sprintf("%s_\\d{3}\\.RDS$", cp_type),
#   full.names = TRUE
# ) |>
#   lapply(readRDS)


# #### Screening Plots ####
#
# # plot screening results
# screen_res_all_plt <- screen_res_all |>
#   # add rep number
#   group_by(n_years_per_block, change_after_year) |>
#   mutate(rep = row_number()) |>
#   ungroup() |> # stack metrics into one column
#   pivot_longer(c(frob, inf, spec)) |>
#   # scale between 0 and 1
#   group_by(name, n_years_per_block) |>
#   mutate(value = boot::inv.logit(scale(value))) |>
#   ungroup() |>
#   mutate(
#     name = case_when(
#       name == "frob" ~ "Frobenius",
#       name == "inf" ~ "Infinity",
#       TRUE ~ "Spectral"
#     )
#   )
#
# p_screen1 <- screen_res_all_plt |>
#   ggplot(
#     aes(x = change_after_year, y = value, colour = factor(n_years_per_block))
#   ) +
#   geom_smooth(show.legend = FALSE) +
#   # geom_line(aes(group = interaction(name, rep)), alpha = 0.2) +
#   geom_line(aes(group = interaction(n_years_per_block, rep)), alpha = 0.2) +
#   # facet_wrap(~n_years_per_block) +
#   facet_wrap(~name) +
#   labs(colour = "Years per block") +
#   cecl_theme() +
#   scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(colour = guide_legend(override.aes = list(linewidth = 6, alpha = 1)))
#
# p_screen1
#
# ggsave(
#   filename = paste0("plots/01_sim/screen_facet_norm_", cp_type, ".png"),
#   plot = p_screen1,
#   width = 8,
#   height = 6
# )
#
# p_screen2 <- screen_res_all_plt |>
#   filter(name != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = value, colour = name)) +
#   geom_smooth(show.legend = FALSE) +
#   geom_line(aes(group = interaction(n_years_per_block, rep)), alpha = 0.2) +
#   facet_wrap(~n_years_per_block) +
#   labs(colour = "Norm Type") +
#   cecl_theme() +
#   scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(colour = guide_legend(override.aes = list(linewidth = 6, alpha = 1)))
#
# p_screen2
#
# ggsave(
#   filename = paste0("plots/01_sim/screen_facet_years_", cp_type, ".png"),
#   plot = p_screen2,
#   width = 8,
#   height = 6
# )
#
#
# #### Changepoint Plotting ####
#
# changepoint_df <- bind_rows(
#   lapply(
#     changepoint_res_all,
#     `[[`,
#     "summary"
#   ),
#   .id = "simulation_id"
# ) |>
#   mutate(
#     simulation_id = as.integer(
#       simulation_id
#     )
#   )
#
# # For each norm and change_after_year, calculate times p-value < 0.05 and < 0.1
# rejection_df <- changepoint_df |>
#   group_by(change_after_year, norm) |>
#   summarise(
#     pval_0.05 = sum(p_value < 0.05) / n(),
#     pval_0.10 = sum(p_value < 0.1) / n(),
#     n = n(),
#     .groups = "drop"
#   )
#
# rejection_df_plt <- rejection_df |>
#   # filter(norm != "Spectral") |>
#   pivot_longer(cols = c(pval_0.05, pval_0.10), names_to = "pval_threshold", values_to = "rejection_rate") |>
#   mutate(pval_threshold = recode(pval_threshold, pval_0.05 = "p < 0.05", pval_0.10 = "p < 0.1")) |>
#   # convert numeric year to date object
#   mutate(change_after_year = as.Date(paste0(change_after_year, "-01-01")))
#
# # find average rejection rate for each norm and pval_threshold (ignores multiple comparisons)
# rejection_df_plt |>
#   group_by(norm, pval_threshold) |>
#   summarise(avg_rejection_rate = mean(rejection_rate), .groups = "drop") |>
#   arrange(norm, pval_threshold)
#
# # find rejection rate for years before and after 1990 (when change occurs)
# rejection_df_plt |>
#   filter(change_after_year >= as.Date("1990-01-01")) |>
#   group_by(norm, pval_threshold) |>
#   summarise(avg_rejection_rate = mean(rejection_rate), .groups = "drop") |>
#   arrange(norm, pval_threshold)
#
# rejection_df_plt |>
#   filter(change_after_year < as.Date("1990-01-01")) |>
#   group_by(norm, pval_threshold) |>
#   summarise(avg_rejection_rate = mean(rejection_rate), .groups = "drop") |>
#   arrange(norm, pval_threshold)
#
# # plot
# rejection_df_plt |>
#   filter(norm != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = rejection_rate, color = pval_threshold)) +
#   geom_point(size = 2) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
#   facet_wrap(~norm) +
#   labs(
#     x = "Season Year",
#     y = "Rejection Rate",
#     colour = "P-value Threshold"
#   ) +
#   cecl_theme() +
#   scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# # barplot of number of changepoints detected for each change_after_year,
# # for each norm type
# p_cp1 <- rejection_df_plt |>
#   filter(norm != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = rejection_rate, fill = pval_threshold)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_abline(aes(slope = 0, intercept = 0.05), linewidth = 1, linetype = "dashed", color = "black") +
#   geom_abline(aes(slope = 0, intercept = 0.1), linewidth = 1, linetype = "dashed", color = "black") +
#   facet_wrap(~norm) +
#   labs(x = "Season year", y = "Rejection Rate", fill = "P-value Threshold") +
#   cecl_theme() +
#   scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p_cp1
#
# # save
# ggsave(
#   filename = paste0("plots/01_sim/changepoint_rejection_rate_", cp_type, ".png"),
#   plot = p_cp1,
#   width = 8,
#   height = 6
# )
#
# # calculate pointwise Type I error (rejection rate at each candidate year)
# pointwise_error_1_df <- changepoint_df |>
#   group_by(
#     change_after_year,
#     norm
#   ) |>
#   summarise(
#     type1_05 = mean(
#       p_value < 0.05,
#       na.rm = TRUE
#     ),
#     type1_10 = mean(
#       p_value < 0.10,
#       na.rm = TRUE
#     ),
#     n_valid = sum(!is.na(p_value)),
#     .groups = "drop"
#   )
# pointwise_error_1_df
#
# # calculate familywise Type I error (rejection rate across all candidate years)
# global_error_1_df <- changepoint_df |>
#   group_by(
#     simulation_id,
#     norm
#   ) |>
#   summarise(
#     min_p_value = min(
#       p_value,
#       na.rm = TRUE
#     ),
#     reject_any_05 = any(
#       p_value < 0.05,
#       na.rm = TRUE
#     ),
#     reject_any_10 = any(
#       p_value < 0.10,
#       na.rm = TRUE
#     ),
#     .groups = "drop"
#   ) |>
#   group_by(norm) |>
#   summarise(
#     familywise_type1_05 =
#       mean(reject_any_05),
#     familywise_type1_10 =
#       mean(reject_any_10),
#     n_simulations = n(),
#     .groups = "drop"
#   )
#
# global_error_1_df
#
#
# #### Changepoint Plotting (old) ####
#
# # add rep number
# changepoint_df_plt <- changepoint_df |>
#   group_by(norm, change_after_year) |>
#   mutate(rep = row_number()) |>
#   ungroup()
#
#
# # add rep number
# changepoint_df_plt <- changepoint_df |>
#   group_by(norm, change_after_year) |>
#   mutate(rep = row_number()) |>
#   ungroup()
#
# # plot p-values
# p_cp1 <- changepoint_df_plt |>
#   filter(norm != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = p_value, colour = norm)) +
#   geom_smooth(show.legend = FALSE) +
#   labs(x = "Season year", y = "p-value", colour = "Norm type") +
#   cecl_theme() +
#   scale_x_continuous(breaks = seq(1985, 1995, by = 1)) +
#   scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(colour = guide_legend(override.aes = list(linewidth = 6, alpha = 1)))
# p_cp1
#
# ggsave(
#   filename = paste0("plots/01_sim/changepoint_pval_", cp_type, ".png"),
#   plot = p_cp1,
#   width = 8,
#   height = 6
# )
#
# p_cp2 <- changepoint_df_plt |>
#   filter(norm != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = p_value, colour = norm)) +
#   geom_smooth(show.legend = FALSE) +
#   geom_line(aes(group = interaction(norm, rep)), alpha = 0.2) +
#   geom_hline(
#     aes(yintercept = 0.05),
#     linetype = "dashed",
#     alpha = 1,
#     colour = "grey10"
#   ) +
#   geom_hline(
#     aes(yintercept = 0.1),
#     linetype = "dashed",
#     alpha = 1,
#     colour = "darkgreen"
#   ) +
#   facet_wrap(~norm) +
#   labs(x = "Season year", y = "p-value", colour = "Norm type") +
#   cecl_theme() +
#   scale_x_continuous(breaks = seq(1985, 1995, by = 1)) +
#   scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(colour = guide_legend(override.aes = list(linewidth = 6, alpha = 1)))
# p_cp2
#
# ggsave(
#   filename = paste0("plots/01_sim/changepoint_pval_facet_", cp_type, ".png"),
#   plot = p_cp2,
#   width = 8,
#   height = 6
# )
#
# # plot changepoint locations (for alpha = 0.05 and alpha = 0.1)
# changepoint_df_bar_plt <- changepoint_df_plt |>
#   mutate(
#     cp_0.05 = ifelse(p_value <= 0.05, 1, 0),
#     cp_0.10 = ifelse(p_value <= 0.10, 1, 0)
#   ) |>
#   pivot_longer(c(cp_0.05, cp_0.10), names_to = "alpha", values_to = "cp") |>
#   mutate(
#     alpha = ifelse(alpha == "cp_0.05", "5% significance", "10% significance"),
#     alpha = factor(alpha, levels = c("5% significance", "10% significance"))
#   ) |>
#   group_by(norm, change_after_year, alpha) |>
#   summarise(cp = sum(cp), .groups = "drop")
#
# lims <- c(0, max(changepoint_df_bar_plt$cp) + 1)
#
# # barplot of number of changepoints detected for each change_after_year,
# # for each norm type
# p_cp3 <- changepoint_df_bar_plt |>
#   filter(norm != "Spectral") |>
#   ggplot(aes(x = change_after_year, y = as.integer(cp), fill = norm)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~alpha) +
#   labs(x = "Season year", y = "# detected changepoints", fill = "Norm type") +
#   cecl_theme() +
#   scale_y_continuous(
#     limits = lims, breaks = seq(0, max(changepoint_df_bar_plt$cp) + 1, by = 1)
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p_cp3
#
# ggsave(
#   filename = paste0("plots/01_sim/changepoint_bar_", cp_type, ".png"),
#   plot = p_cp3,
#   width = 8,
#   height = 6
# )
#
# # barplot with percentages
# p_cp4 <- changepoint_df_bar_plt |>
#   filter(norm != "Spectral") |>
#   group_by(norm, change_after_year, alpha) |>
#   mutate(cp_perc = cp / length(changepoint_res_all)) |>
#   ungroup() |>
#   ggplot(aes(x = change_after_year, y = cp_perc, fill = norm)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~alpha) +
#   labs(
#     x = "Season year",
#     y = "Percentage of detected changepoints",
#     fill = "Norm type"
#   ) +
#   cecl_theme() +
#   scale_y_continuous(labels = scales::percent) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p_cp4
#
# ggsave(
#   filename = paste0("plots/01_sim/changepoint_bar_perc_", cp_type, ".png"),
#   plot = p_cp4,
#   width = 8,
#   height = 6
# )

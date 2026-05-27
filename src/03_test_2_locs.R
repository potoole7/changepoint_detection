#### Run changepoint detection process for two locations, one model ####

# For each difference in correlation:
# 1. Apply screening approach to find candidate changepoints
# - For 100 simulations, calculate and create line and boxplots
# - Determine how many times "peak" is at known changepoint
# 2. Use permutation test to test points around peak
# - For 100 simulations, create line and boxplots
# - Also calculate the Type II Error & Power of test at various rho vals (done)

# TODO Run for diff = 1.5 and 2.5, as power seems to scale up considerably between
# 0.2 and 0.3 (and less so between 0.1 and 0.2)

# TODO Read changepoint literature and see what other analyses they do.
# Some ideas (from Chat) are:
# - heatmap of p-values for different diff values (done)
# - Global type I error over all times (done)
# - Changepoint estimation accuracy (Bias, MAE, RMSE, coverage)
# - Sharpness of peak for screening: contrast/ratio
# - Sensitivity analysis on threshold(s), as well as diff values
#   - Also for window size and permutations?
# - Histogram of p-values under no difference: Should be uniform! (also KS test?)
# - Estimation of computational cost vs grid size and # permutations

# Also make sure you're using producing the same number of exceedances from
# both blocks when thresholding! (done)

# TODO Test both methods for performing screening test to see which works better
# across all simulations!
# - Method 1 using same threshold for both blocks, method 2 using same
# number of exceedances for both blocks (so different thresholds)

# TODO Increase number of simulations if possible (get on Nimbus!)

# TODO Use nicer, more sensible functions and code!

# TODO Allow last CE param values to be used in next step
# Will break embarassingly parallel code, but fine for these simulations
# as we're not running the tests in parallel, just the simulations!



#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(grid) # for unit()
library(RColorBrewer)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(parallel)
library(evgam)
library(ggridges)
library(patchwork)
source("src/00_functions.R")

#### metadata ####

seed <- 123 # not great but passable
# seed <- 12345 # worse
# seed <- 321
# seed <- 16366303
n_vars <- 2 # number of variables
n_locs <- 2
n_clust <- 2
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per loc times n_locs per group
n_per_block <- n * n_clust # number of obs per block (two blocks here)
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.9) # (initial) correlation parameters for each group

# laplace_q <- 0.9
laplace_q <- 0.95

# diff_vec <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5) # differences in correlation for second block
diff_vec <- c(0, 0.05, 0.1, 0.2, 0.3) # differences in correlation for second block
mc_cores <- detectCores() - 1 # number of cores for parallel processing
cond_var <- "X2" # variable on RHS of conditioning event

# number of simulations to test
nsims <- 100

# vector of block positions to explore (i.e., size of first block)
block_vec <- seq(500, 1500)

# vector of block positions to calculate permutation test for
# grid_vals <- seq(950, 1050, by = 10) # TODO Add more values once code works
grid_vals <- seq(800, 1200, by = 20)

#### Functions ####

# function to generate data from
gen_dat <- \(cor_t_vec, diff, laplace_trans = TRUE) {
  data_loc <- bind_rows(lapply(seq_along(diff), \(i) {
    cor_spec <- cor_t_vec
    # increase first correlation value by diff
    cor_spec[[1]] <- cor_spec[[1]] + diff[[i]]
    print(cor_spec)
    # generate data, with optional Laplace transformation
    ret <- do.call(
      rbind,
      lapply(cor_spec, gen_t, n = n, laplace_trans = laplace_trans)
    )
    ret |>
      mutate(
        name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
        block = i # time block
      )
  }))
  return(data_loc)
}

# set.seed(seed)
# data_loc <- gen_dat(
#   cor_t_vec,
#   diff = c(0, diff_vec)
#   # diff = c(0, last(diff_vec)), # TODO Temp, just do for largest difference
#   laplace_trans = TRUE
# )


#### Pre-calculations ####

# use the same dependence threshold across all variables and locations
# Choose larger quantile so it's always (?) larger than block-specific thresholds!
dep_val <- qlaplace(0.8) # for Laplace marginals
# dep_val <- qlaplace(0.95) # for Laplace marginals

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
  n = 500, # TODO Increase later
  thresh_max = dep_val,
  y_max = laplace_cap
)


#### Screening ####

set.seed(seed)
norm_vals_lst <- mclapply(seq_len(nsims), \(j) {
  # norm_vals_lst <- mclapply(1:3, \(j) { # test
  # norm_vals_lst <- lapply(1:3, \(j) { # test

  # generate data
  data_loc <- gen_dat(
    cor_t_vec,
    diff = c(0, diff_vec) # control + all diff values
  )

  # create list where each diff_vec value is separated into a separate entry
  data_loc_lst <- lapply(unique(data_loc$block)[-1], \(i) {
    data_loc |>
      filter(block %in% c(1, i)) |>
      mutate(block = ifelse(block != 1, 2, block))
  })

  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(j / nsims, 3) * 100, "% completed", collapse = "")
  ))

  # lopp through diff values
  norm_vals_lst <- lapply(diff_vec, \(diff) {
    # pull data specific to current diff_vec value
    # print(paste("diff =", diff))
    data_loc <- data_loc_lst[[which(diff_vec == diff)]]

    # Calculate start values
    norm1 <- single_run_explore(
      data_loc,
      block_vec[(length(block_vec) / 2) + 1],
      cond_val = dep_val,
      laplace_sample = laplace_sample,
      aLow = 0,
      nruns = 3,
      ret_dep = TRUE
    )
    dep <- norm1$dep
    start_vals <- lapply(dep, coef)

    # at each block_vec value, calculate the test statistics (i.e. "screen")
    norm_vals <- mclapply(block_vec, \(i) {
      x <- single_run_explore(
        data_loc,
        i,
        cond_var = "X2",
        cond_prob = 0.8,
        use_dth = TRUE,
        laplace_sample = laplace_sample,
        aLow = 0,
        nruns = 2,
        start = start_vals
      )
      # x <- single_run_explore(
      #   data_loc,
      #   i,
      #   cond_val = dep_val,
      #   laplace_sample = laplace_sample,
      #   aLow = 0,
      #   nruns = 3,
      #   start = start_vals
      # )
      x
    })
  })
})

# save output
saveRDS(
  norm_vals_lst,
  file = "data/norm_vals_lst.rds"
)

# load output
norm_vals_lst <- readRDS("data/norm_vals_lst.rds")

#### Plot ####

# extract data
norm_vals_frob_df <- bind_rows(lapply(norm_vals_lst, \(dat) {
  bind_rows(lapply(seq_along(dat), \(i) {
    norm_vals <- dat[[i]]
    data.frame(
      diff = diff_vec[i],
      time = block_vec,
      D    = sapply(norm_vals, \(x) x$frob)
    )
  }))
}), .id = "iteration")

# line plots and boxplots
p_line1 <- norm_vals_frob_df |>
  ggplot(aes(x = time, y = D, colour = as.factor(diff))) +
  # geom_point(alpha = 0.8) +
  geom_line(aes(group = iteration), alpha = 0.1) +
  geom_smooth(alpha = 1, linewidth = 1) +
  facet_wrap(~diff) +
  labs(
    x = "time",
    y = "D",
    colour = "Correlation increase"
  ) +
  cecl_theme() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8, alpha = 1))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(plot = p_line1, "screen_plot_line.png", height = 8, width = 10)

# boxplot
# TODO Thin times! Or bin later??
p_box1 <- norm_vals_frob_df |>
  # cut in 10s so boxplots aren't too large!
  filter(time %% 50 == 0) |>
  # filter(time %% 50 == 0) |>
  # ggplot(aes(x = time, y = D, colour = as.factor(diff), group = iteration)) +
  ggplot(aes(x = as.factor(time), y = D, colour = as.factor(diff))) +
  geom_boxplot() +
  geom_vline(
    xintercept = as.factor(1000),
    # colour = "black",
    # linewidth = 1,
    linetype = "dashed"
  ) +
  facet_wrap(~diff) +
  labs(
    x = "time",
    y = "D",
    colour = "Correlation increase"
  ) +
  cecl_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(plot = p_box1, "screen_plot_box.png", height = 8, width = 13)

# how many times does "peak" correspond to changepoint in iterations?
# For different "sub-windows"
# TODO Also do within windows, i.e. 1-10, 11-20, etc

find_peak_val <- \(df, min_time = NULL, max_time = NULL) {
  df_ret <- df
  if (!is.null(min_time)) {
    df_ret <- filter(df_ret, time >= min_time)
  }
  if (!is.null(max_time)) {
    df_ret <- filter(df_ret, time <= max_time)
  }
  df_ret |>
    group_by(iteration, diff) |>
    filter(D == max(D)) |>
    ungroup() |>
    count(diff, time) |>
    mutate(perc = n / nsims) |>
    arrange(diff, desc(perc)) |>
    group_by(diff) |>
    slice(1) |>
    ungroup()
}

# default min and max are 500 and 1500, respectively
find_peak_val(norm_vals_frob_df)
find_peak_val(norm_vals_frob_df, 700, 1300)

# norm_vals_frob_df |>
#   group_by(iteration, diff) |>
#   filter(D == max(D)) |>
#   ungroup() |>
#   count(diff, time) |>
#   mutate(perc = n / nsims) |>
#   arrange(diff, desc(perc)) |>
#   group_by(diff) |>
#   slice(1)


#### Permutation test ####

# Apply permutation test over timepoints around peak
set.seed(seed)
perm_test_res_lst <- mclapply(seq_len(nsims), \(j) {
  # perm_test_res_lst <- mclapply(1:3, \(j) { # test
  # perm_test_res_lst <- lapply(1:3, \(j) { # test

  # re-generate data, too memory intensive to store from last step
  data_loc <- gen_dat(
    cor_t_vec,
    diff = c(0, diff_vec)
  )
  data_loc_lst <- lapply(unique(data_loc$block)[-1], \(i) {
    data_loc |>
      filter(block %in% c(1, i)) |>
      mutate(block = ifelse(block != 1, 2, block))
  })

  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(j / nsims, 3) * 100, "% completed", collapse = "")
  ))

  # perform permutation test for each diff and grid value
  perm_res <- lapply(data_loc_lst, \(data_loc) {
    perm_test_fun(
      data = select(data_loc, X1, X2, name, block),
      grid_vals = grid_vals,
      n_per_block = min(grid_vals), # use max data available
      n_perm = 100,
      laplace_trans = TRUE,
      ret_dep = TRUE,
      aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
      cond_prob = 0.8, # use block-specific threshs to ensure same exceedances
      laplace_sample = laplace_sample # control Laplace sample
    )
  })
})

# save output
saveRDS(
  perm_test_res_lst,
  file = "data/perm_test_res_lst.rds"
)

# load output
perm_test_res_lst <- readRDS("data/perm_test_res_lst.rds")

#### Plot ####

# collate and plot using histogram
preprocess_fun <- \(x, grid_vals = NULL, n_per_block = NULL) {
  bind_rows(lapply(seq_along(x), \(j) {
    ret <- with(
      x[[j]],
      data.frame(
        norm_value = c(perm_norms_inf, perm_norms_frob),
        norm_orig = c(
          rep(norm_orig_inf, length(perm_norms_inf)),
          rep(norm_orig_frob, length(perm_norms_frob))
        ),
        p_value = c(
          rep(p_value_inf, length(perm_norms_inf)),
          rep(p_value_frob, length(perm_norms_frob))
        ),
        norm_type = rep(c("Infinity", "Frobenius"), each = length(perm_norms_inf)) # ,
        # method = names(perm_test_res)[i],
        # grid_val = grid_vals[[j]]
      )
    )
    if (!is.null(grid_vals)) {
      ret$grid_val <- grid_vals[[j]]
    }
    if (!is.null(n_per_block)) {
      ret$n_per_block <- n_per_block[[j]]
    }
    ret
  }))
}

# extract
perm_plot_df <- bind_rows(lapply(perm_test_res_lst, \(x) {
  lapply(x, preprocess_fun, grid_vals = grid_vals) |>
    bind_rows(.id = "iteration") |>
    mutate(iteration = as.numeric(iteration)) |>
    filter(norm_type == "Frobenius") |>
    # replace iteration column with correct diff values
    left_join(
      norm_vals_frob_df |>
        distinct(diff) |>
        mutate(iteration = row_number()),
      by = "iteration"
    ) |>
    select(-iteration)
}), .id = "iteration")

# extract just p_values
p_value_df <- perm_plot_df |>
  group_by(diff, iteration, norm_type, grid_val) |>
  slice(1) |>
  select(-norm_value) |>
  ungroup()

# boxplot of p_values
p_box_perm <- p_value_df |>
  ggplot(aes(x = as.factor(grid_val), y = p_value, colour = as.factor(diff))) +
  geom_boxplot() +
  facet_wrap(~diff) +
  labs(x = "time", y = "p", colour = "Correlation Increase") +
  cecl_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(plot = p_box_perm, "perm_plot_box.png", height = 8, width = 13)

# heatmap of p_values
p_heat <- p_value_df |>
  ggplot(aes(x = grid_val, y = factor(diff), fill = p_value)) +
  geom_tile() +
  labs(x = "time", y = "Changepoint Magnitude", fill = "p") +
  cecl_theme(nejm_pal = FALSE, legend.position = "right") +
  scale_fill_viridis_c(
    option = "A",
    name = "p-value",
  ) +
  theme(
    legend.title = element_text(size = 14, face = "bold")
  )
p_heat
ggsave(plot = p_heat, "perm_plot_heat.png", height = 6, width = 8)

# calculate type II error & power and plot
# TODO Could also get type I error for other grid values??
type2_err_df <- p_value_df |>
  # true changepoint
  filter(grid_val == 1000) |>
  mutate(
    accept_null = p_value > 0.05,
    reject_null = p_value <= 0.05
  ) |>
  # aggregate across simulations
  group_by(diff) |>
  summarise(
    n = n(),
    type2_err = mean(accept_null),
    power = mean(reject_null),
    se_type2_err = sqrt(type2_err * (1 - type2_err) / n),
    se_power = sqrt(power * (1 - power) / n),
    ci_low_type2 = type2_err - 1.96 * se_type2_err,
    ci_high_type2 = type2_err + 1.96 * se_type2_err,
    ci_low_power = power - 1.96 * se_power,
    ci_high_power = power + 1.96 * se_power,
    .groups = "drop"
  )

p_type2 <- type2_err_df %>%
  # clamp error rates
  mutate(
    ymin = pmax(0, ci_low_type2),
    ymax = pmin(1, ci_high_type2)
  ) %>%
  ggplot(aes(x = diff, y = type2_err)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.02) +
  labs(
    x = "Correlation Increase",
    y = "Type II error"
  ) +
  cecl_theme()
p_type2

p_power <- type2_err_df %>%
  # clamp error rates
  mutate(
    ymin = pmax(0, ci_low_power),
    ymax = pmin(1, ci_high_power)
  ) %>%
  ggplot(aes(x = diff, y = power)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.02) +
  labs(
    # x = "Correlation Increase",
    x = "Changepoint Magnitude",
    y = "Power"
  ) +
  cecl_theme()
p_power
ggsave(plot = p_power, "perm_plot_power.png", height = 6, width = 8)

# Now calculate type I error for all other grid points across diff values!
type1_err_df <- p_value_df |>
  # remove true changepoint
  filter(grid_val != 1000) |>
  mutate(
    accept_null = p_value > 0.05,
    reject_null = p_value <= 0.05
  ) |>
  # aggregate across simulations
  group_by(diff, grid_val) |>
  summarise(
    n = n(),
    type1_err = mean(reject_null),
    se_type1 = sqrt(type1_err * (1 - type1_err) / n),
    ci_low = pmax(0, type1_err - 1.96 * se_type1),
    ci_high = pmin(1, type1_err + 1.96 * se_type1),
    .groups = "drop"
  )

p_type1 <- type1_err_df %>%
  ggplot(aes(x = grid_val, y = type1_err, colour = as.factor(diff))) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(
    xintercept = 1000,
    colour = "black",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.1
  ) +
  facet_wrap(~diff) +
  labs(
    x = "time",
    y = "Type I error",
    colour = "Changepoint Magnitude"
  ) +
  cecl_theme(legend.position = "right") +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.95, 0.15),
    legend.justification = c(1, 0),
    # legend.background = element_rect(fill = "white", colour = "black")
  )
p_type1

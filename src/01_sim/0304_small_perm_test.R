#### Simple permutation test for two sites ####

# TODO But what about the Laplace sample? What happens there when we change
# the threshold used for each block s.t. the # exceedances is the same in both?

# TODO Idea: Could maybe use only 500 or so points to left and right of block
# position, that way the effect of having more points in one block compared
# to the other is reduced!

# TODO Write conclusions from this!

# TODO Why are p-values for Frob and Inf norms the same??

# TODO Output distribution of alpha, beta values to make boxplots with!
# (i.e. output for each permutation)

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
source("src/00_functions.R")

#### Metadata ####

seed <- 123
n_vars <- 2 # number of variables
# n_locs <- 12 # number of locations
n_locs <- 2
# n_locs <- 60
# n_clust <- 3 # number of clusters
n_clust <- 2
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per loc times n_locs per group
n_per_block <- n * n_clust # number of obs per block (two blocks here)
df_t <- 3 # degrees of freedom for t-copula
# cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group
cor_t_vec <- c(0.1, 0.9) # (initial) correlation parameters for each group

# diff <- 0 # try for no difference in correlation (but stochastically different data!
# diff <- 0.05
diff <- 0.1
# diff <- 0.1
# diff <- 0.2
# diff <- 0.4
# diff <- 0.3
# use_test_case <- TRUE # try for absolutely no difference (i.e. same exact data in both blocks)
use_test_case <- FALSE
mc_cores <- detectCores() - 1 # number of cores for parallel processing
# mc_cores <- 14
# cond_var <- NULL
cond_var <- "X2" # variable on RHS of conditioning event
grid_vals <- seq(950, 1050, by = 10) # TODO Add more values once code works

#### Generate data ####

gen_dat <- \(cor_t_vec, increments, diff) {
  data_loc <- bind_rows(lapply(seq_along(increments), \(i) {
    cor_spec <- cor_t_vec
    cor_spec[[1]] <- cor_spec[[1]] + (diff * increments[[i]])
    print(cor_spec)
    # generate data
    ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n))
    ret |>
      mutate(
        name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
        block = i # time block
      )
  }))
  return(data_loc)
}


increments <- c(0, 1)
set.seed(seed)
data_loc <- gen_dat(cor_t_vec, increments, diff)

# use same data for both blocks if test case
data_loc_test <- data_loc
if (use_test_case && diff == 0) {
  data_loc_test <- data_loc |>
    filter(block == 1) |>
    bind_rows(
      data_loc |>
        filter(block == 1) |>
        mutate(block = 2)
    )
}

# check data is correct
table(data_loc_test$block)
table(data_loc_test$name)

# should look the same
data_loc_test |>
  ggplot() +
  geom_point(aes(X1, X2)) +
  facet_wrap(name ~ block) +
  cecl_theme()

# Should be the same if test case, but different if not
data_loc_test |>
  group_by(name, block) |>
  summarise(cor = cor(X1, X2), .groups = "drop")



#### Precalculations ####

# use the same dependence threshold across all variables and locations
dep_val <- qlaplace(0.8) # for Laplace marginals

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

# Also make a longer Laplace sample
# laplace_sample2 <- rlaplace_trunc(
#   n = 2000,
#   thresh_max = dep_val,
#   y_max = laplace_cap
# )


#### Run test ####

# Double check that function works properly the whole way through (no bugs) (debugonce)

# # n_perm <- 100 # TODO Increase
# n_perm <- 1000
# # n_perm <- 1
#
# # source("src/00_functions.R")
# set.seed(123)
# perm_test_res <- perm_test_fun(
#   data = select(data_loc_test, X1, X2, name, block),
#   grid_vals = grid_vals,
#   n_perm = n_perm,
#   laplace_trans = TRUE,
#   ret_dep = TRUE,
#   cond_var = cond_var,
#   aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#   cond_val = dep_val, laplace_sample = laplace_sample # control Laplace sample
# )
# #
# # # save
# # saveRDS(perm_test_res, paste0("small_perm_test_diff_", diff, ".rds"))
#
# # load
# # perm_test_res <- readRDS(paste0("small_perm_test_diff_", diff, ".rds"))

# function to pull out important information from permutation test results
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

# function to produce plots of interest
plot_fun <- \(perm_test_res, grid_vals, plot_ce_params = TRUE) {
  # pull out important information
  perm_test_res_df <- preprocess_fun(perm_test_res, grid_vals = grid_vals)

  # points of original norm values
  cols <- colorRampPalette(brewer.pal(12, "Set3"))
  myPal <- cols(length(unique(perm_test_res_df$grid_val)))
  p1 <- perm_test_res_df |>
    mutate(grid_val_facet = paste0(grid_val, " (p = ", p_value, ")")) |>
    mutate(
      grid_val_facet = factor(grid_val_facet, levels = unique(grid_val_facet))
    ) |>
    # ggplot(aes(x = norm_value, fill = factor(grid_val))) +
    ggplot(aes(x = norm_value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      fill = "red",
      position = "identity",
      alpha = 0.5,
      bins = 30
    ) +
    geom_vline(
      # aes(xintercept = norm_orig, colour = factor(grid_val)),
      aes(xintercept = norm_orig),
      colour = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(x = "Norm", y = "Density") +
    # scale_fill_manual(values = myPal) +
    # scale_colour_manual(values = myPal) +
    # facet_wrap(norm_type ~ grid_val_facet, scales = "free") +
    facet_wrap(grid_val_facet ~ norm_type, scales = "free") +
    cecl_theme(nejm_pal = FALSE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
    )

  # boxplots
  p2 <- perm_test_res_df |>
    ggplot(aes(x = grid_val, y = norm_value)) +
    geom_boxplot(aes(group = grid_val), width = 2, outliers = TRUE, key_glyph = "rect") +
    # geom_hline(
    #   data = perm_test_res_df |>
    #     distinct(grid_val, norm_orig, norm_type),
    #   aes(yintercept = norm_orig, group = grid_val),
    #   linetype = "dashed",
    #   size = 1
    # ) +
    geom_point(aes(y = norm_orig), colour = "blue", size = 3) +
    # geom_boxplot(aes(group = grid_val)) +
    # geom_smooth() +
    facet_wrap(~norm_type, scales = "free") +
    cecl_theme()

  ret <- list(
    hist_plot = p1,
    box_plot = p2
  )

  if (!plot_ce_params) {
    return(ret)
  }

  # plot dependence values
  perm_test_dep_df <- bind_rows(lapply(seq_along(perm_test_res), \(j) {
    perm_test_res[[j]]$dep_vals |>
      as.data.frame() |>
      mutate(grid_val = grid_vals[[j]])
  }))

  p3 <- perm_test_dep_df |>
    pivot_longer(cols = c("a", "b", "m", "s", "ll"), names_to = "param", values_to = "value") |>
    mutate(param = ifelse(param == "ll", "LogLik", param)) |>
    arrange(name, block, param) |>
    mutate(param = paste0(param, ", ", name, ", block ", block)) |>
    ggplot(aes(x = grid_val, y = value, colour = name)) +
    geom_point() +
    geom_line() +
    # facet_wrap(param ~ name + block, scales = "free", ncol = 4) +
    facet_wrap(~param, scales = "free", ncol = 4) +
    cecl_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  return(
    list(
      hist_plot = p1,
      box_plot = p2,
      param_plot = p3
    )
  )
}

# plots <- plot_fun(perm_test_res, grid_vals)
#
# pdf(paste0("hist_diff_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("hist_diff_", diff, "_not_test.pdf"), width = 10, height = 8)
# plots[[1]] + ggtitle(paste0("Diff = ", diff))
# # plots[[1]] + ggtitle(paste0("Diff = ", diff, " (not test case)"))
# dev.off()
#
# pdf(paste0("box_diff_", diff, ".pdf"), width = 10, height = 6)
# # pdf(paste0("box_diff_", diff, "not_test.pdf"), width = 10, height = 6)
# plots[[2]] + ggtitle(paste0("Diff = ", diff))
# # plots[[2]] + ggtitle(paste0("Diff = ", diff, " (not test case)"))
# dev.off()
#
# pdf(paste0("ce_params_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("ce_params_", diff, "not_test.pdf"), width = 10, height = 8)
# plots[[3]] + ggtitle(paste0("Diff = ", diff))
# # plots[[3]] + ggtitle(paste0("Diff = ", diff, " (not test case)"))
# dev.off()

# conclusion:
# For diff = 0; very little difference between permutation distributions,
# but original norm values can differ significantly, which gives much
# different p-values
# Hence, significant p-values at grid values ???
# when no change point is actually present

# For diff = 0.4;

# For diff = 0.2;

#### Test with equal points in each block ####

# Run again with equal number of points in each block,
# to check that results are not just due to more points in one block than the
# other (as this can affect the CE estimation)

# (n_per_block_perm <- min(grid_vals))
#
# set.seed(123)
# perm_test_res_n_per_block <- perm_test_fun(
#   data = select(data_loc_test, X1, X2, name, block),
#   grid_vals = grid_vals,
#   # n_perm = n_perm,
#   n_perm = n_perm,
#   n_per_block = n_per_block_perm,
#   laplace_trans = TRUE,
#   cond_var = cond_var,
#   ret_dep = TRUE,
#   aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#   cond_val = dep_val, laplace_sample = laplace_sample # control Laplace sample
# )
#
# saveRDS(perm_test_res_n_per_block, paste0("small_perm_test_diff_", diff, "_n_per_block_", n_per_block_perm, "not_test.rds"))
#
# # perm_test_res_n_per_block <- readRDS(paste0("small_perm_test_diff_", diff, "_n_per_block_", n_per_block_perm, ".rds"))
#
# plots2 <- plot_fun(perm_test_res_n_per_block, grid_vals)
#
# pdf(paste0("hist2_diff_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("hist2_diff_", diff, "_not_test.pdf"), width = 10, height = 8)
# plots2[[1]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[1]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()
#
# pdf(paste0("box2_diff_", diff, ".pdf"), width = 10, height = 6)
# # pdf(paste0("box2_diff_", diff, "_not_test.pdf"), width = 10, height = 6)
# plots2[[2]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[2]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()
#
# pdf(paste0("ce_params2_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("ce_params2_", diff, "_not_test.pdf"), width = 10, height = 8)
# plots2[[3]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[3]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()

#### Test w/ equal points in each, ENSURING N EXCEEDANCES IN EACH BLOCK ####

# TODO Ensure that threshold taken for both blocks is the same (and max value)
# (see `use_dth` argument)

# TODO check # thresholded obs in each block is actually the same (and equal to 1 - cond_prob)

# diff <- 0 # try for no difference in correlation (but stochastically different data!
# diff <- 0.05
# diff <- 0.1
# diff <- 0.1
# diff <- 0.2
# diff <- 0.4
# diff <- 0.3

# use_test_case <- FALSE # try for absolutely no difference (i.e. same exact data in both blocks)
#
# data_loc <- gen_dat(cor_t_vec, increments, diff)
#
# # use same data for both blocks if test case
# data_loc_test <- data_loc
# if (diff == 0 && use_test_case) {
#   data_loc_test <- data_loc |>
#     filter(block == 1) |>
#     bind_rows(
#       data_loc |>
#         filter(block == 1) |>
#         mutate(block = 2)
#     )
# }
#
# (n_per_block_perm <- min(grid_vals))
# n_perm <- 1000
# # n_perm <- 100
#
# set.seed(123)
# perm_test_res_thresh_block <- perm_test_fun(
#   data = select(data_loc_test, X1, X2, name, block),
#   grid_vals = grid_vals,
#   n_perm = n_perm,
#   n_per_block = n_per_block_perm,
#   laplace_trans = TRUE,
#   cond_var = cond_var,
#   # define thresh as m-th largest order statistic s.t. P(X > thresh) = cond_prob, to ensure n exceedances in each block
#   cond_prob = 0.8,
#   ret_dep = TRUE,
#   aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#   use_dth = TRUE, # Use different threshold for each block
#   # cond_val = dep_val,
#   laplace_sample = laplace_sample # control Laplace sample
# )
#
# saveRDS(perm_test_res_thresh_block, paste0("small_perm_test_diff_", diff, "_n_per_block_", n_per_block_perm, "_thresh_blocks.rds"))
# # # perm_test_res_thresh_block <- readRDS(paste0("small_perm_test_diff_", diff, "_n_per_block_", n_per_block_perm, ".rds"))
# plots2 <- plot_fun(perm_test_res_thresh_block, grid_vals)
#
# pdf(paste0("hist3_diff_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("hist2_diff_", diff, "_not_test.pdf"), width = 10, height = 8)
# plots2[[1]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[1]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()
#
# pdf(paste0("box3_diff_", diff, ".pdf"), width = 10, height = 6)
# # pdf(paste0("box2_diff_", diff, "_not_test.pdf"), width = 10, height = 6)
# plots2[[2]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[2]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()
#
# pdf(paste0("ce_params3_", diff, ".pdf"), width = 10, height = 8)
# # pdf(paste0("ce_params2_", diff, "_not_test.pdf"), width = 10, height = 8)
# plots2[[3]] + ggtitle(paste0("Diff = ", diff, ", npoints = ", n_per_block_perm))
# # plots2[[3]] + ggtitle(paste0("Diff = ", diff, " (not test case), npoints = ", n_per_block_perm))
# dev.off()

#### Test with different block sizes ####

# run with different values of n_per_perm, for the same

# n_per_block_vals <- seq(200, 900, by = 100)
n_per_block_vals <- seq(100, 1000, by = 50)

# set.seed(123)
# perm_test_res_block <- mclapply(n_per_block_vals, \(n_per_block_perm) {
#   perm_test_fun(
#     data = select(data_loc_test, X1, X2, name, block),
#     grid_vals = 1000,
#     n_perm = n_perm,
#     n_per_block = n_per_block_perm,
#     laplace_trans = TRUE,
#     cond_var = cond_var,
#     ret_dep = TRUE,
#     aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#     cond_val = dep_val, laplace_sample = laplace_sample # control Laplace sample
#   )
# })
#
# saveRDS(perm_test_res_block, paste0("small_perm_test_diff_", diff, "_n_per_block_vals.rds"))


#### Line plots for all values of diff ####

diff_vals <- c(0, 0.05, 0.075, 0.1, 0.2, 0.3)

## experiment 1 (different grid values) ##
# files <- list.files(pattern = "small_perm_test_diff", full.names = TRUE)
# files <- files[!grepl("n_per_block", files)]
# # extract numeric diff from filename
# file_diffs <- as.numeric(
#   sub(".*small_perm_test_diff_([0-9.]+)\\.rds$", "\\1", files)
# )
# # reorder files to match diff_vals
# files <- files[match(diff_vals, file_diffs)]
#
# perm_lst <- lapply(files, readRDS)
# perm_df <- bind_rows(lapply(seq_along(perm_lst), \(i) {
#   preprocess_fun(perm_lst[[i]], grid_vals = grid_vals) |>
#     mutate(diff = diff_vals[[i]])
# }))
#
# perm_df |>
#   filter(norm_type == "Frobenius") |>
#   ggplot(aes(x = grid_val, y = p_value, colour = as.factor(diff))) +
#   geom_point() +
#   geom_line() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Grid value",
#     y = "p-value",
#     colour = "Correlation increase"
#   ) +
#   cecl_theme() # TODO Do for second experiment
#
#
# ## experiment 2 (different grid values, same block sizes) ##
# files2 <- list.files(pattern = "small_perm_test_diff_[0-9.]+_n_per_block_[0-9]+\\.rds", full.names = TRUE)
# file_diffs2 <- as.numeric(
#   sub(".*small_perm_test_diff_([0-9.]+)_n_per_block_[0-9]+\\.rds$", "\\1", files2)
# )
# files2 <- files2[match(diff_vals, file_diffs2)]
#
# perm_block_lst <- lapply(files2, readRDS)
# perm_block_df <- bind_rows(lapply(seq_along(perm_block_lst), \(i) {
#   preprocess_fun(perm_block_lst[[i]], grid_vals = grid_vals) |>
#     mutate(diff = diff_vals[[i]])
# }))
#
# perm_block_df |>
#   filter(norm_type == "Frobenius") |>
#   ggplot(aes(x = grid_val, y = p_value, colour = as.factor(diff))) +
#   geom_point() +
#   geom_line() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Grid value",
#     y = "p-value",
#     colour = "Correlation increase"
#   ) +
#   cecl_theme()
#
# ## experiment 3 (same grid value, vary (same) block size) ##
#
# # list files
# files3 <- list.files(pattern = "n_per_block_vals.rds", full.names = TRUE)
# file_diffs3 <- as.numeric(
#   sub(".*small_perm_test_diff_([0-9.]+)_n_per_block_vals\\.rds$", "\\1", files3)
# )
# files3 <- files3[match(diff_vals, file_diffs3)]
#
# perm_lst_block3 <- lapply(files3, readRDS)
#
# # pull out important information and plot (flatten list as only one grid val)
# perm_lst_block3 <- lapply(perm_lst_block3, \(x) {
#   preprocess_fun(purrr::flatten(x), n_per_block = n_per_block_vals)
# })
#
# perm_block_df3 <- bind_rows(lapply(seq_along(perm_lst_block3), \(i) {
#   mutate(perm_lst_block3[[i]], diff = diff_vals[[i]])
# }))
#
# # plot
# perm_block_df3 |>
#   filter(norm_type == "Frobenius") |>
#   ggplot(aes(x = n_per_block, y = p_value, colour = as.factor(diff))) +
#   geom_point() +
#   geom_line() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value",
#     colour = "Correlation increase"
#   ) +
#   ggtitle("Permutation test p-values for different block sizes") +
#   scale_x_continuous(breaks = n_per_block_vals) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   cecl_theme()
#
# ## experiment 4 (no difference, different seed numbers)
#
# seed_vals <- 123:132
#
# files4 <- list.files(pattern = "diff_0_seed_", full.names = TRUE)
#
# perm_lst_block4 <- lapply(files4, readRDS)
#
# # pull out important information and plot (flatten list as only one grid val)
# perm_lst_block4 <- lapply(perm_lst_block4, \(x) {
#   preprocess_fun(purrr::flatten(x), n_per_block = n_per_block_vals)
# })
#
# perm_block_df4 <- bind_rows(lapply(seq_along(perm_lst_block4), \(i) {
#   # mutate(perm_lst_block4[[i]], diff = diff_vals[[i]])
#   mutate(perm_lst_block4[[i]], diff = 0, seed = seed_vals[[i]])
# }))
#
# # plot
# perm_block_df4 |>
#   filter(norm_type == "Frobenius") |>
#   # ggplot(aes(x = n_per_block, y = p_value, colour = as.factor(seed))) +
#   ggplot(aes(x = n_per_block, y = p_value)) +
#   # geom_point(colour = "red", alpha = 0.3) +
#   geom_line(aes(group = as.factor(seed)), colour = "red", alpha = 0.3) +
#   geom_smooth() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value"
#   ) +
#   ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
#   scale_x_continuous(breaks = n_per_block_vals) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   cecl_theme(nejm_pal = FALSE)


#### Run sliding window test case many times ####

# # n_per_block_vals <- seq(200, 900, by = 100)
# n_per_block_vals <- seq(100, 1000, by = 50)
# diff_test <- 0
#
# set.seed(123)
# perm_test_res_loop <- mclapply(1:100, \(i) {
#   system(sprintf("echo %s", paste("repetition", i)))
#   data_loc <- gen_dat(cor_t_vec, increments, diff_test)
#
#   # use same data for both blocks if test case
#   data_loc_test <- data_loc
#   if (diff_test == 0) {
#     data_loc_test <- data_loc |>
#       filter(block == 1) |>
#       bind_rows(
#         data_loc |>
#           filter(block == 1) |>
#           mutate(block = 2)
#       )
#   }
#
#   lapply(n_per_block_vals, \(n_per_block_perm) {
#     perm_test_fun(
#       data = select(data_loc_test, X1, X2, name, block),
#       grid_vals = 1000,
#       n_perm = n_perm,
#       n_per_block = n_per_block_perm,
#       laplace_trans = TRUE,
#       cond_var = cond_var,
#       ret_dep = TRUE,
#       aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#       cond_val = dep_val, laplace_sample = laplace_sample # control Laplace sample
#     )
#   })
# })
#
# saveRDS(perm_test_res_loop, "perm_test_res_loop.rds")

#### Plot results of sliding window ####

# load
# files_slide <- list.files(pattern = "perm_test_res_loop_", full.names = TRUE)
# files_slide <- files_slide[!grepl("loop2", files_slide)]
# file_not_test <- files_slide[grepl("not_test", files_slide)]
# files_slide <- files_slide[files_slide != file_not_test]
# diffs_slide <- as.numeric(
#   sub(".*perm_test_res_loop_diff_([0-9.]+)\\.rds$", "\\1", files_slide)
# )
# files_slide <- files_slide[match(diff_vals, diffs_slide)]
# files_slide <- c(file_not_test, files_slide)
#
# data_slide <- lapply(files_slide, readRDS)
#
# # pull out important information and plot (flatten list as only one grid val)
# # perm_lst_slide <- lapply(data_slide, \(x) {
# #   preprocess_fun(purrr::flatten(x), n_per_block = n_per_block_vals)
# # })
# # names(data_slide) <- c("not_test", diff_vals)
# names(data_slide) <- c("0", "Exact Same", diff_vals[2:length(diff_vals)])
# diff_vals_name <- names(data_slide)
# # x <- data_slide[[1]]
# data_slide_df <- bind_rows(lapply(seq_along(data_slide), \(i) {
#   x <- data_slide[[i]]
#   bind_rows(lapply(x, \(y) {
#     preprocess_fun(purrr::flatten(y), n_per_block = n_per_block_vals)
#   }), .id = "rep") |>
#     mutate(diff = diff_vals_name[[i]])
# })) |>
#   relocate(rep, .after = n_per_block) |>
#   mutate(diff = factor(diff, levels = c("Exact Same", "0", as.character(diff_vals[2:length(diff_vals)]))))
#
# readr::write_csv(
#   data_slide_df,
#   "data_slide_df.csv"
# )
#
# # plot
# p <- data_slide_df |>
#   filter(norm_type == "Frobenius") |>
#   ggplot(aes(x = n_per_block, y = p_value)) +
#   geom_line(aes(group = as.factor(rep)), colour = "red", alpha = 0.3) +
#   geom_smooth() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value"
#   ) +
#   ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
#   # scale_x_continuous(breaks = n_per_block_vals) +
#   scale_x_continuous(breaks = n_per_block_vals[seq(1, length(n_per_block_vals), by = 2)]) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   cecl_theme(nejm_pal = FALSE) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   facet_wrap(~diff)
# ggsave(plot = p, filename = "plot_not_test.png", height = 10, width = 12)
#
# # also plot smooth lines for each, without red lines in the background
# p2 <- data_slide_df |>
#   # filter(norm_type == "Frobenius") |>
#   # group_by(diff, n_per_block) |>
#   # slice(1:10) |>
#   # ungroup() |>
#   ggplot(aes(x = n_per_block, y = p_value, colour = diff, group = diff)) +
#   geom_smooth() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value",
#     colour = "Correlation increase"
#   ) +
#   scale_x_continuous(breaks = n_per_block_vals) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
#   cecl_theme() +
#   guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# ggsave(plot = p2, filename = "plot2_not_test.png", height = 10, width = 12)
#
# rm(p1, p2, data_slide)
# gc()


### Run sliding window test case many times (WITH SAME DATA!!) ####

# # n_per_block_vals <- seq(200, 900, by = 100)
# n_per_block_vals <- seq(100, 1000, by = 50)
# # n_per_block_vals <- seq(950, 1000, by = 50)
# # diff_test <- 0
# # diff_test <- 0.05
# # diff_test <- 0.075
# diff_test <- 0.1
# # diff_test <- 0.2
#
# set.seed(123)
# data_loc <- gen_dat(cor_t_vec, increments, diff_test)
#
# # use same data for both blocks if test case
# data_loc_test <- data_loc
# if (diff_test == 0) {
#   data_loc_test <- data_loc |>
#     filter(block == 1) |>
#     bind_rows(
#       data_loc |>
#         filter(block == 1) |>
#         mutate(block = 2)
#     )
# }
#
# # perm_test_res_loop_same <- mclapply(1:100, \(i) {
# perm_test_res_loop_same <- mclapply(1:2, \(i) {
#   system(sprintf("echo %s", paste("repetition", i)))
#
#   lapply(n_per_block_vals, \(n_per_block_perm) {
#     perm_test_fun(
#       data = select(data_loc_test, X1, X2, name, block),
#       grid_vals = 1000,
#       n_perm = n_perm,
#       n_per_block = n_per_block_perm,
#       laplace_trans = TRUE,
#       cond_var = cond_var,
#       ret_dep = TRUE,
#       aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
#       cond_val = dep_val, laplace_sample = laplace_sample # control Laplace sample
#     )
#   })
# })
#
# # saveRDS(perm_test_res_loop, "perm_test_res_loop.rds")
# saveRDS(perm_test_res_loop_same, paste0("perm_test_res_loop2_diff_", diff_test, ".rds"))

# check that they've run properly!!
# perm_test_res_loop_same_test <- readRDS(
#   # "perm_test_res_loop2_diff_0.05.rds" # fine
#   "perm_test_res_loop2_diff_0.075.rds"  # error
# )

#### Plot ####

# # load
# files_slide2 <- list.files(pattern = "perm_test_res_loop2", full.names = TRUE)
# file_not_test2 <- files_slide2[grepl("not_test", files_slide2)]
# files_slide2 <- files_slide2[files_slide2 != file_not_test2]
# # TODO temp, remove afterwards (also fix, not running!)
# # files_slide2 <- files_slide2[grepl("diff_0.rds", files_slide2)]
# # files_slide2 <- files_slide2[!grepl("0.05", files_slide2)]
# files_slide2 <- files_slide2[!grepl("0.075", files_slide2)]
# diffs_slide2 <- as.numeric(
#   sub(".*perm_test_res_loop2_diff_([0-9.]+)\\.rds$", "\\1", files_slide2)
# )
# matches <- match(diff_vals, diffs_slide2)
# files_slide2 <- files_slide2[matches]
# # use regular expressions to extract number
# # (diff_vals_spec <- diff_vals[sort(matches[!is.na(matches)])])
# diff_vals_spec <- as.numeric(
#   sub(".*perm_test_res_loop2_diff_([0-9.]+)\\.rds$", "\\1", files_slide2)
# )
# (diff_vals_spec <- diff_vals_spec[!is.na(diff_vals_spec)])
#
# # add back in non-test file
# files_slide2 <- c(file_not_test2, files_slide2)
# (files_slide2 <- files_slide2[!is.na(files_slide2)]) # rm diffs not ran
#
# # load data
# data_slide2 <- lapply(files_slide2, readRDS)
#
# # pull out important information and plot (flatten list as only one grid val)
# # perm_lst_slide <- lapply(data_slide2, \(x) {
# #   preprocess_fun(purrr::flatten(x), n_per_block = n_per_block_vals)
# # })
# # names(data_slide2) <- c("not_test", diff_vals)
# # TODO Temp, remove afterwards (also fix, not running!)
# # names(data_slide2) <- c("0", "Exact Same", diff_vals_spec[2:length(diff_vals_spec)])
# names(data_slide2) <- c("0", "Exact Same", diff_vals_spec[2:length(diff_vals_spec)])
# # names(data_slide2) <- c("0", "Exact Same")
# (diff_vals_name2 <- names(data_slide2))
# # x <- data_slide2[[1]]
#
# data_slide2_df <- bind_rows(lapply(seq_along(data_slide2), \(i) {
#   print(i)
#   print(diff_vals_name2[[i]])
#   x <- data_slide2[[i]]
#   bind_rows(lapply(x, \(y) {
#     preprocess_fun(purrr::flatten(y), n_per_block = n_per_block_vals)
#   }), .id = "rep") |>
#     mutate(diff = diff_vals_name2[[i]])
# })) |>
#   relocate(rep, .after = n_per_block) |>
#   mutate(diff = factor(diff, levels = c("Exact Same", "0", as.character(diff_vals[2:length(diff_vals)]))))
# readr::write_csv(
#   data_slide2_df,
#   "data_slide2_df.csv"
# )
#
# # plot
# p <- data_slide2_df |>
#   filter(norm_type == "Frobenius") |>
#   ggplot(aes(x = n_per_block, y = p_value)) +
#   geom_line(aes(group = as.factor(rep)), colour = "red", alpha = 0.3) +
#   geom_smooth() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value"
#   ) +
#   ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
#   # scale_x_continuous(breaks = n_per_block_vals) +
#   scale_x_continuous(breaks = n_per_block_vals[seq(1, length(n_per_block_vals), by = 2)]) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   cecl_theme(nejm_pal = FALSE) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   facet_wrap(~diff)
# ggsave(plot = p, filename = "plot_same.png", height = 10, width = 12)
#
# # also plot smooth lines for each, without red lines in the background
# p2 <- data_slide2_df |>
#   # filter(norm_type == "Frobenius") |>
#   # group_by(diff, n_per_block) |>
#   # slice(1:10) |>
#   # ungroup() |>
#   ggplot(aes(x = n_per_block, y = p_value, colour = diff, group = diff)) +
#   geom_smooth() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#   labs(
#     x = "Number of points per block",
#     y = "p-value",
#     colour = "Correlation increase"
#   ) +
#   scale_x_continuous(breaks = n_per_block_vals) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
#   cecl_theme() +
#   guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# ggsave(plot = p2, filename = "plot_same2.png", height = 10, width = 12)

#### Run sliding window test case many times (W/ SAME DATA, SAME EXCEEDANCES) ####

n_per_block_vals <- seq(100, 1000, by = 50)
diff_test <- 0
diff_test <- 0.05
diff_test <- 0.075
diff_test <- 0.1
diff_test <- 0.2
diff_test <- 0.3
use_test_case <- FALSE

set.seed(123)
perm_test_res_loop_thresh_blocks <- mclapply(1:100, \(i) {
  system(sprintf("echo %s", paste("repetition", i)))
  data_loc <- gen_dat(cor_t_vec, increments, diff_test)

  # use same data for both blocks if test case
  data_loc_test <- data_loc
  if (diff_test == 0 && use_test_case == TRUE) {
    data_loc_test <- data_loc |>
      filter(block == 1) |>
      bind_rows(
        data_loc |>
          filter(block == 1) |>
          mutate(block = 2)
      )
  }

  lapply(n_per_block_vals, \(n_per_block_perm) {
    perm_test_fun(
      data = select(data_loc_test, X1, X2, name, block),
      grid_vals = 1000,
      n_perm = n_perm,
      n_per_block = n_per_block_perm,
      laplace_trans = TRUE,
      cond_var = cond_var,
      cond_prob = 0.8, # define thresh as m-th largest order statistic s.t. P(X > thresh) = cond_prob, to ensure n exceedances in each block
      ret_dep = TRUE,
      use_dth = TRUE,
      aLow = 0, nruns = 2, use_start = TRUE, # optimal specification of CE
      # cond_val = dep_val,
      laplace_sample = laplace_sample # control Laplace sample
    )
  })
})

# saveRDS(perm_test_res_loop_thresh_blocks, "perm_test_res_loop_thresh_blocks.rds")
saveRDS(
  perm_test_res_loop_thresh_blocks,
  # "perm_test_res_loop_thresh_blocks.rds"
  paste0("perm_test_res_loop3_thresh_blocks_diff_", diff_test, ".rds")
)

#### Plot ####

# TODO Functionalise, repeated throughout!

# load
files_slide3 <- list.files(pattern = "perm_test_res_loop3", full.names = TRUE)
file_not_test3 <- files_slide3[grepl("not_test", files_slide3)]
if (length(file_not_test3) > 0) {
  files_slide3 <- files_slide3[files_slide3 != file_not_test3]
}
# TODO temp, remove afterwards (also fix, not running!)
# files_slide3 <- files_slide3[grepl("diff_0.rds", files_slide3)]
# files_slide3 <- files_slide3[!grepl("0.05", files_slide3)]
# files_slide3 <- files_slide3[!grepl("0.075", files_slide3)]
# files_slide3 <- files_slide3[!grepl("0.rds", files_slide3)]
diffs_slide3 <- as.numeric(
  # sub(".*perm_test_res_loop3_diff_([0-9.]+)\\.rds$", "\\1", files_slide3)
  sub(".*perm_test_res_loop3_thresh_blocks_diff_([0-9.]+)\\.rds$", "\\1", files_slide3)
)
matches <- match(diff_vals, diffs_slide3)
files_slide3 <- files_slide3[matches]
# use regular expressions to extract number
# (diff_vals_spec <- diff_vals[sort(matches[!is.na(matches)])])
diff_vals_spec <- as.numeric(
  # sub(".*perm_test_res_loop3_diff_([0-9.]+)\\.rds$", "\\1", files_slide3)
  sub(".*perm_test_res_loop3_thresh_blocks_diff_([0-9.]+)\\.rds$", "\\1", files_slide3)
)
(diff_vals_spec <- diff_vals_spec[!is.na(diff_vals_spec)])

# add back in non-test file
files_slide3 <- c(file_not_test3, files_slide3)
(files_slide3 <- files_slide3[!is.na(files_slide3)]) # rm diffs not ran

# load data
data_slide3 <- lapply(files_slide3, readRDS)

# pull out important information and plot (flatten list as only one grid val)
# perm_lst_slide <- lapply(data_slide3, \(x) {
#   preprocess_fun(purrr::flatten(x), n_per_block = n_per_block_vals)
# })
# names(data_slide3) <- c("not_test", diff_vals)
# TODO Temp, remove afterwards (also fix, not running!)
# names(data_slide3) <- c("0", "Exact Same", diff_vals_spec[3:length(diff_vals_spec)])
# names(data_slide3) <- c("0", "Exact Same", diff_vals_spec[3:length(diff_vals_spec)])
# names(data_slide3) <- c("0", diff_vals_spec[3:length(diff_vals_spec)])
names(data_slide3) <- diff_vals_spec
# names(data_slide3) <- c("0", "Exact Same")
(diff_vals_name3 <- names(data_slide3))
# x <- data_slide3[[1]]

data_slide3_df <- bind_rows(lapply(seq_along(data_slide3), \(i) {
  print(i)
  print(diff_vals_name3[[i]])
  x <- data_slide3[[i]]
  bind_rows(lapply(x, \(y) {
    preprocess_fun(purrr::flatten(y), n_per_block = n_per_block_vals)
  }), .id = "rep") |>
    mutate(diff = diff_vals_name3[[i]])
})) |>
  relocate(rep, .after = n_per_block) |>
  # mutate(diff = factor(diff, levels = c("Exact Same", "0", as.character(diff_vals[3:length(diff_vals)]))))
  mutate(diff = factor(diff, levels = as.character(diff_vals_spec)))
readr::write_csv(
  data_slide3_df,
  paste0("data_slide3_df.csv")
)

# plot
p <- data_slide3_df |>
  filter(norm_type == "Frobenius") |>
  ggplot(aes(x = n_per_block, y = p_value)) +
  geom_line(aes(group = as.factor(rep)), colour = "red", alpha = 0.3) +
  geom_smooth() +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
  labs(
    x = "Number of points per block",
    y = "p-value"
  ) +
  ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
  # scale_x_continuous(breaks = n_per_block_vals) +
  scale_x_continuous(breaks = n_per_block_vals[seq(1, length(n_per_block_vals), by = 2)]) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  cecl_theme(nejm_pal = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~diff)
ggsave(plot = p, filename = "plot_thresh_blocks.png", height = 10, width = 13)

# also plot smooth lines for each, without red lines in the background
p2 <- data_slide3_df |>
  # filter(norm_type == "Frobenius") |>
  # group_by(diff, n_per_block) |>
  # slice(1:10) |>
  # ungroup() |>
  ggplot(aes(x = n_per_block, y = p_value, colour = diff, group = diff)) +
  geom_smooth() +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
  labs(
    x = "Number of points per block",
    y = "p-value",
    colour = "Correlation increase"
  ) +
  scale_x_continuous(breaks = n_per_block_vals) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
  cecl_theme() +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(plot = p2, filename = "plot_thresh_blocks2.png", height = 10, width = 12)

#### Comparison plot ####


# TODO Also compare with output from previous method (different # exceedances)
df_compare <- bind_rows(
  readr::read_csv("data_slide_df.csv") |> mutate(test = "Different lengths, different exceedances"),
  readr::read_csv("data_slide2_df.csv") |> mutate(test = "Same length, different exceedances"),
  readr::read_csv("data_slide3_df.csv") |>
    mutate(
      test = "Same length, same exceedances",
      diff = as.character(diff)
    )
) |>
  filter(norm_type == "Frobenius")

p_compare <- df_compare |>
  ggplot(aes(x = n_per_block, y = p_value, colour = test, group = test)) +
  # geom_line(aes(group = as.factor(rep)), colour = "red", alpha = 0.3) +
  geom_smooth() +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
  facet_wrap(~diff) +
  labs(
    x = "Number of points per block",
    y = "p-value"
  ) +
  # ggtitle("Permutation test p-values for control case, for different block sizes and seeds") +
  # scale_x_continuous(breaks = n_per_block_vals) +
  scale_x_continuous(breaks = n_per_block_vals[seq(1, length(n_per_block_vals), by = 2)]) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  cecl_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  plot = p_compare,
  filename = "plot_compare.png",
  height = 10,
  width = 13
)

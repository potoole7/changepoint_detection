#### Testing CE fitting stability ####

# TODO For one location and one model, arrange data by size of conditioning
# variable and sense check that when we move non-extreme points we get the
# same CE fit (for the first ~80% of points)

# TODO For models where only one point is swapped, shouldn't CE fits
# be the same for all locations 80% of the time?
# So the divergence should also be the same? Interested in looking into this!!
# Important here to swap just one point for first location, but
# to look at divergence for *all* locations, as before we onl tested for
# one

# TODO Change print for `cecl_dep` to print variables and conditioning vars
# (do the same for marginal and divergence!)

# CE fit seems to vary *a lot* during permutation tests, even for the same
# data in both blocks!
# Investigate by only moving each individual point and seeing how much
# CE fit changes, and how much the distance between CEs changes

# TODO Look into why norm values aren't the same 80% of the time, even
# when CE parameter values are
# NOTE It may be that alpha, beta, mu, sigma estimates do not necessarily have
# to be the same "runs"
# TODO Look at just one location and one model for the divergence, as
# that's what we do when looking at the difference in CE fit!
# TODO Only look at runs where alpha, beta are the same, and see what norm
# values are like there!
# TODO May be more important to do where we only move one single point, rather
# than a point from each locatinn

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(evgam) # for alternative CE fitting
library(ggplot2)
library(parallel)

source("src/00_functions.R")

#### metadata ####

seed <- 123
n_vars <- 2 # number of variables
n_locs <- 12 # number of locations
# n_locs <- 60
n_clust <- 3 # number of clusters
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per loc times n_locs per group
n_per_block <- n * n_clust # number of obs per block (two blocks here)
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group

# diff <- 0.4 # difference in correlation for local change
diff <- 0 # try for no difference, see what happens!
mc_cores <- detectCores() - 1 # number of cores for parallel processing

#### simulate data ####

# increase lowest correlation by 0.05 until it reaches middle correlation
set.seed(seed)
# start from 0 to have control example (checks difference due to seed)
# increments <- c(0, 1) # TODO Add more increments once code works
increments <- 0
data_loc <- bind_rows(lapply(seq_along(increments), \(i) {
  cor_spec <- cor_t_vec
  cor_spec[[1]] <- cor_spec[[1]] + (diff * increments[[i]])
  print(cor_spec)
  # generate data
  ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n, laplace_trans = TRUE))
  ret |>
    mutate(
      name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
      block = i # time block
    )
}))

data_loc_test <- data_loc |>
  filter(block == 1) |>
  bind_rows(
    data_loc |>
      filter(block == 1) |>
      mutate(block = 2)
  )

#### Precalculations ####

dep_val <- qlaplace(0.8) # for Laplace marginals

# Calculate cap for sampling from Laplace distribution
laplace_cap <- qlaplace(0.99)

# function to generate Laplace samples
rlaplace_trunc <- \(n, thresh_max, y_max) {
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


#### Test 1: Check if modelling is deterministic ####

# fit CE model and compute distance using CeCl 100 times
n_reps <- 100

set.seed(seed)
# dep_lst <- lapply(seq_len(n_reps), \(i) {
dep_lst <- mclapply(seq_len(n_reps), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_reps)))
  data_laplace_block <- trans_fun(data_loc_test, n_vars = 2, laplace_trans = TRUE)
  dep <- lapply(
    data_laplace_block,
    cecl_dep,
    cond_val    = dep_val,
    fit_no_keef = TRUE
  )
}, mc.cores = mc_cores)

dep_coef_df <- bind_rows(mclapply(dep_lst, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}), .id = "iteration") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

dep_coef_df |>
  group_by(model, block, name) |>
  distinct(a, b) |>
  ungroup() |>
  count(model, block, name) |>
  arrange(desc(n))
# all 1s, so model fitting is deterministic

# But what about if we reordered the samples?
set.seed(123)
dep_lst_reord <- mclapply(seq_len(n_reps), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_reps)))
  data_loc_test_reord <- data_loc_test |>
    group_by(name, block) |>
    slice_sample(prop = 1) |>
    ungroup()

  # test to see if they're the same
  # data.frame(data_loc_test$X1, data_loc_test_reord$X1) |>
  #   cor() |>
  #   print()

  data_laplace_block <- trans_fun(
    data_loc_test_reord,
    n_vars = 2, laplace_trans = TRUE
  )
  dep <- lapply(
    data_laplace_block,
    cecl_dep,
    cond_val    = dep_val,
    fit_no_keef = TRUE
  )
})

dep_coef_df_reord <- bind_rows(mclapply(dep_lst_reord, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}), .id = "iteration") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

dep_coef_df_reord |>
  group_by(model, block, name) |>
  distinct(a, b) |>
  ungroup() |>
  count(model, block, name) |>
  arrange(desc(n))
# all 1s, so model fitting is also deterministic even when reordering samples


#### Test 2: Swap one point and inspect CE ####

# Move single point (from each location) from block 1 to block 2
# Then refit the CE models

# first, just move a single point from block 1 to block 2 from each location
dep_lst_mv <- mclapply(seq_len(n_per_loc), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  data_loc_mv <- data_loc_test |>
    group_by(name, block) |>
    mutate(row = row_number()) |>
    ungroup() |>
    mutate(
      block = ifelse(block == 1 & row == i, 2, block)
    ) |>
    select(-row)

  # # check
  # data_loc_mv |>
  #   count(block, name) |>
  #   print(n = Inf)

  # proceed as before
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  lapply(
    data_laplace_block,
    cecl_dep,
    cond_val    = dep_val,
    fit_no_keef = TRUE
  )
})

# save results
# saveRDS(dep_lst_mv, "dep_lst_mv.RDS")
# dep_lst_mv <- readRDS("dep_lst_mv.RDS")

# also run with aLow, nruns to see if that makes it more stable
dep_lst_mv2 <- mclapply(seq_len(n_per_loc), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  data_loc_mv <- data_loc_test |>
    group_by(name, block) |>
    mutate(row = row_number()) |>
    ungroup() |>
    mutate(
      block = ifelse(block == 1 & row == i, 2, block)
    ) |>
    select(-row)

  # # check
  # data_loc_mv |>
  #   count(block, name) |>
  #   print(n = Inf)

  # proceed as before
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  lapply(
    data_laplace_block,
    cecl_dep,
    cond_val = dep_val,
    fit_no_keef = TRUE,
    aLow = 0,
    nruns = 2
  )
})

# join coefficients
dep_coef_df_mv <- bind_rows(mclapply(dep_lst_mv, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}), .id = "row") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

# plot histogram of different estimates (for loc_1 first)
dep_coef_df_mv_spec <- dep_coef_df_mv |>
  filter(name == "loc_1", model == "X1 | X2") |>
  # select(block, a, b) |>
  # select(block, a, b, m, s, ll) |>
  select(block, a, b, m, s, ll, row) |>
  # pivot_longer(cols = c(a, b), names_to = "parameter") |>
  pivot_longer(cols = c(a, b, m, s, ll), names_to = "parameter") |>
  # mutate(value = round(value, 5)) |>
  mutate(parameter = ifelse(parameter == "ll", "LogLik", parameter)) |>
  identity()

dep_coef_df_mv_spec |>
  ggplot(aes(x = value, fill = block)) +
  geom_histogram(position = "dodge", bins = 30) +
  # facet_wrap(parameter ~ block, scales = "free") +
  facet_wrap(parameter ~ block, scales = "free", ncol = 5) +
  labs(
    # title = "alpha, beta for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point",
    title = "CE params for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point",
    # x = "Parameter estimate (5 s.f)",
    x = "Parameter estimate"
  ) +
  cecl_theme()

# count percentage for each value of a and b
dep_coef_df_mv_spec |>
  group_by(block, parameter, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(block, parameter) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  arrange(desc(pct), parameter, block) |>
  head(n = 10)
# Around 80% of CE estimates are the same, which agrees with using
# 80th quantile of Laplace distribution as conditioning value!

# TODO find rows where alpha and beta are the same for both blocks,
# to check if norm values are the same there
# first, find modes for CE parameters
mode_fun <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

dep_coef_df_mode <- dep_coef_df_mv_spec |>
  group_by(block, parameter) |>
  filter(value == mode_fun(value)) |>
  ungroup()

# only keep those that are the same for all parameters
dep_coef_df_mode <- dep_coef_df_mode |>
  group_by(block, row) |>
  filter(n_distinct(value) == 5) |>
  ungroup()

mode_rows <- as.numeric(unique(dep_coef_df_mode$row))
# check that the length is about 800
length(mode_rows)

# also do for "better" fitting with aLow, nruns
dep_coef_df_mv2 <- bind_rows(mclapply(dep_lst_mv2, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}), .id = "row") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

# plot histogram of different estimates (for loc_1 first)
dep_coef_df_mv2_spec <- dep_coef_df_mv2 |>
  filter(name == "loc_1", model == "X1 | X2") |>
  # select(block, a, b) |>
  select(block, a, b, m, s, ll) |>
  # pivot_longer(cols = c(a, b), names_to = "parameter") |>
  pivot_longer(cols = c(a, b, m, s, ll), names_to = "parameter") |>
  # mutate(value = round(value, 5)) |>
  mutate(parameter = ifelse(parameter == "ll", "LogLik", parameter)) |>
  identity()

dep_coef_df_mv2_spec |>
  ggplot(aes(x = value, fill = block)) +
  geom_histogram(position = "dodge", bins = 30) +
  # facet_wrap(parameter ~ block, scales = "free") +
  facet_wrap(parameter ~ block, scales = "free", ncol = 5) +
  labs(
    # title = "alpha, beta for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point",
    title = "''Better'' X1 | X2, loc_1, 80th Laplace quantile, moving 1 point",
    # x = "Parameter estimate (5 s.f)",
    x = "Parameter estimate"
  ) +
  cecl_theme()

# count percentage for each value of a and b
dep_coef_df_mv2_spec |>
  group_by(block, parameter, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(block, parameter) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  arrange(desc(pct), parameter, block) |>
  head(n = 10)
# Around 80% of CE estimates are the same, which agrees with using
# 80th quantile of Laplace distribution as conditioning value!

#### Divergences ####

# Now calculate divergences for each fit
# TODO Streamline more with previous code, repeated bits
# could make first part it's own function, rest would be the same as before
div_lst_mv <- mclapply(seq_along(dep_lst_mv), \(i) {
  # div_lst_mv <- lapply(seq_along(dep_lst_mv), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  dep_spec <- dep_lst_mv[[i]]
  # also need to pull marginals
  data_loc_mv <- data_loc_test |>
    group_by(name, block) |>
    mutate(row = row_number()) |>
    ungroup() |>
    mutate(
      block = ifelse(block == 1 & row == i, 2, block)
    ) |>
    select(-row)

  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  dist <- lapply(seq_along(dep_spec), \(j) {
    cecl_dist(
      dep_spec[[j]],
      data_laplace_block[[j]],
      # n_mc = 500,
      # seed = seed # TODO See if the laplace samples are the same each time?
      laplace_sample = laplace_sample
    )
  })

  # calculate norms for each individual distance matrix, to check if they're
  # the same for each run
  frob_dist <- sapply(dist, \(x) {
    norm(
      as.matrix(x$dist_mat),
      type = "F"
    )
  })
  inf_dist <- sapply(dist, \(x) {
    norm(
      as.matrix(x$dist_mat),
      type = "I"
    )
  })

  # calculate norms
  frob <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "F"
  )
  inf <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "I"
  )

  return(list(
    # "dist" = dist, # bit memory intensive!
    "frob_dist1" = frob_dist[[1]],
    "frob_dist2" = frob_dist[[2]],
    "inf_dist1" = inf_dist[[1]],
    "inf_dist2" = inf_dist[[2]],
    "frob" = frob,
    "inf" = inf
  ))
})

# save
saveRDS(div_lst_mv, "div_lst_mv.RDS")
# load
# div_lst_mv <- readRDS("div_lst_mv.RDS")

# Now we want to compare norm values
div_df_mv <- bind_rows(mclapply(div_lst_mv, \(x){
  data.frame(
    frob_dist1 = x$frob_dist1,
    frob_dist2 = x$frob_dist2,
    inf_dist1 = x$inf_dist1,
    inf_dist2 = x$inf_dist2,
    frob = x$frob,
    inf = x$inf
  )
}), .id = "row")

# plot histogram
div_df_mv_long <- div_df_mv |>
  # pivot_longer(c(frob, inf), names_to = "norm") |>
  # pivot_longer(c(frob, frob_dist, inf, inf_dist), names_to = "norm") |>
  pivot_longer(
    c(frob, frob_dist1, frob_dist2, inf, inf_dist1, inf_dist2),
    names_to = "norm"
  ) |>
  # mutate(norm = ifelse(norm == "frob", "Frobenius", "Infinity"))
  mutate(
    norm = case_when(
      norm == "frob" ~ "Frobenius",
      norm == "frob_dist1" ~ "Frobenius (M_1)",
      norm == "frob_dist2" ~ "Frobenius (M_2)",
      norm == "inf" ~ "Infinity",
      norm == "inf_dist1" ~ "Infinity (M_1)",
      norm == "inf_dist2" ~ "Infinity (M_2)"
    )
  )

div_df_mv_long |>
  ggplot(aes(x = value, fill = norm)) +
  geom_histogram(position = "dodge", bins = 30, show.legend = FALSE) +
  labs(
    title = "Norms when moving 1 point (all sites, both models)",
    x = "Norm value"
  ) +
  facet_wrap(~norm, scales = "free") +
  cecl_theme()

# count percentage at each value
div_df_mv_long |>
  filter(norm %in% c("Frobenius", "Infinity")) |>
  mutate(value = round(value, 3)) |>
  group_by(norm, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(norm) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  arrange(desc(pct))

#### Divergence for same CE values ####

# TODO Investigate
#   block   row   name  var   cond_var   dth      a     b     m     s    ll
# 1 Block 1 165   loc_1 X1    X2       0.916 0.0739 0.551    NA    NA    NA
# 2 Block 2 165   loc_1 X1    X2       0.916 0.0739 0.551    NA    NA    NA
# 3 Block 2 578   loc_1 X1    X2       0.916 0.0739 0.551    NA    NA    NA
# 4 Block 1 681   loc_1 X1    X2       0.916 0.0739 0.551    NA    NA    NA


# pull models where CE are the same, convert to `cecl_dep` objects
dep_lst_mode <- dep_coef_df_mode |>
  mutate(
    name = "loc_1",
    var = "X1",
    cond_var = "X2",
    dth = dep_val,
    parameter = ifelse(parameter == "LogLik", "ll", parameter)
  ) |>
  pivot_wider(names_from = parameter, values_from = value) |>
  relocate(row, block, name, var, cond_var) |>
  relocate(dth, .after = everything()) |>
  identity() |>
  group_split(row, .keep = FALSE) |>
  lapply(\(x) {
    group_split(x, block, .keep = FALSE) |>
      lapply(\(y) {
        as_cecl_dep(y)
      })
  })

mode_rows_df <- dep_coef_df_mode |>
  select(row, block) |>
  mutate(block = stringr::str_replace(block, "Block ", "")) |>
  mutate(across(everything(), as.numeric)) |>
  arrange(block, row)

mode_rows1 <- mode_rows_df |>
  filter(block == 1) |>
  pull(row) |>
  unique()
n_mode <- length(mode_rows1)

# div_lst_mode <- mclapply(seq_along(dep_lst_mv), \(i) {
div_lst_mode <- lapply(seq_along(dep_lst_mode), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_mode)))
  dep_spec <- dep_lst_mode[[i]]
  # also need to pull marginals
  data_loc_mv <- data_loc_test |>
    filter(name == "loc_1") |> # only looking at location 1
    group_by(name, block) |>
    mutate(row = row_number()) |>
    ungroup() |>
    mutate(
      # block = ifelse(block == 1 & row == i, 2, block)
      block = ifelse(block == 1 & row == mode_rows1[[i]], 2, block)
    ) |>
    select(-row)

  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  # TODO Check that there is only one divergence matrix, not two
  # (as we're only conditioning on a single value)
  # TODO Check that distance matrix is created correctly
  # dist <- lapply(seq_along(dep_spec), \(j) {
  #   cecl_dist(
  #     dep_spec[[j]],
  #     data_laplace_block[[j]],
  #     n_mc = 500,
  #     seed = seed
  #   )
  # })
  dist <- tryCatch(
    {
      lapply(seq_along(dep_spec), \(j) {
        cecl_dist(
          dep_spec[[j]],
          data_laplace_block[[j]],
          # n_mc = 500,
          # seed = seed
          laplace_sample = laplace_sample
        )
      })
    },
    error = \(e) {
      print(sprintf("Error in iteration %s: %s", i, e$message))
      browser()
      debugonce(cecl_dist)
      lapply(seq_along(dep_spec), \(j) {
        cecl_dist(
          dep_spec[[j]],
          data_laplace_block[[j]],
          # n_mc = 500,
          # seed = seed
          laplace_sample = laplace_sample
        )
      })
    }
  )

  # calculate norms for each individual distance matrix, to check if they're
  # the same for each run
  frob_dist <- sapply(dist, \(x) {
    norm(
      as.matrix(x$dist_mat),
      type = "F"
    )
  })
  inf_dist <- sapply(dist, \(x) {
    norm(
      as.matrix(x$dist_mat),
      type = "I"
    )
  })

  # calculate norms
  frob <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "F"
  )
  inf <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "I"
  )

  return(list(
    # "dist" = dist, # bit memory intensive!
    "frob_dist1" = frob_dist[[1]],
    "frob_dist2" = frob_dist[[2]],
    "inf_dist1" = inf_dist[[1]],
    "inf_dist2" = inf_dist[[2]],
    "frob" = frob,
    "inf" = inf
  ))
})

# extract norm values for each run, check if they're the same
div_df_mode <- bind_rows(mclapply(div_lst_mode, \(x){
  data.frame(
    frob_dist1 = x$frob_dist1,
    frob_dist2 = x$frob_dist2,
    inf_dist1 = x$inf_dist1,
    inf_dist2 = x$inf_dist2,
    frob = x$frob,
    inf = x$inf
  )
}), .id = "row")

# check that all are 0
div_df_mode |>
  select(-row) |>
  unlist() |>
  unique()
# They are!

# So for the same random sample and CE parameters, divergence is 0
# But what about if we only move one single point from one location, rather than one point from each location?
# If we take the divergence across all locations for this setup,
# shouldn't we get about 80% 0s for divergence??

#### Test2b: Swap single point from *one* location ####

# Do for one single point and single model
# TODO Also functionalise
dep_lst_mv_single <- mclapply(seq_len(n_per_loc), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  data_loc_mv <- data_loc_test |>
    mutate(block = ifelse(
      block == 1 & name == "loc_1" & row_number() == i, 2, block
    ))
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )
  lapply(
    data_laplace_block,
    cecl_dep,
    cond_val    = dep_val,
    fit_no_keef = TRUE
  )
}, mc.cores = mc_cores)

# save
saveRDS(dep_lst_mv_single, "dep_lst_mv_single.RDS")
# load
# dep_lst_mv_single <- readRDS("dep_lst_mv_single.RDS")

dep_coef_df_mv_single <- bind_rows(mclapply(dep_lst_mv_single, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}, mc.cores = mc_cores), .id = "row") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

# plot histogram of different estimates (for loc_1 first)
dep_coef_df_mv_single_spec <- dep_coef_df_mv_single |>
  filter(name == "loc_1", model == "X1 | X2") |>
  # select(block, a, b) |>
  select(block, a, b, m, s, ll, row) |>
  # pivot_longer(cols = c(a, b), names_to = "parameter") |>
  pivot_longer(cols = c(a, b, m, s, ll), names_to = "parameter") |>
  mutate(parameter = ifelse(parameter == "ll", "LogLik", parameter)) |>
  # mutate(value = round(value, 5)) |>
  identity()

dep_coef_df_mv_single_spec |>
  ggplot(aes(x = value, fill = block)) +
  geom_histogram(position = "dodge", bins = 30) +
  facet_wrap(parameter ~ block, scales = "free", ncol = 5) +
  labs(
    # title = "alpha, beta for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point",
    # title = "alpha, beta for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point from loc_1 only",
    title = "CE params for X1 | X2, loc_1, 80th Laplace quantile, moving 1 point from loc_1 only",
    # x = "Parameter estimate (5 s.f)",
    x = "Parameter estimate"
  ) +
  cecl_theme()
# alpha estimates the same, but beta estimates less variable (less extreme)

# count percentage for each value of a and b
dep_coef_df_mv_single_spec |>
  group_by(block, parameter, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(block, parameter) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  # arrange(desc(pct)) |>
  arrange(desc(pct), parameter, block) |>
  head(n = 10)
# Again, ~80% of estimates are the same


### Test2c: Arrange by conditioning value ####

#### Divergence when swapping one point from one location ####

# TODO Just make this a function above, don't be lazy!!
div_lst_mv_single <- mclapply(seq_along(dep_lst_mv_single), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  dep_spec <- dep_lst_mv_single[[i]]
  data_loc_mv <- data_loc_test |>
    mutate(block = ifelse(
      block == 1 & name == "loc_1" & row_number() == i, 2, block
    ))
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  dist <- lapply(seq_along(dep_spec), \(j) {
    cecl_dist(
      dep_spec[[j]],
      data_laplace_block[[j]],
      # n_mc = 500,
      # seed = seed
      laplace_sample = laplace_sample
    )
  })
  frob <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "F"
  )
  inf <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "I"
  )
  return(list(
    "frob" = frob,
    "inf" = inf
  ))
})

# save
saveRDS(div_lst_mv_single, "div_lst_mv_single.RDS")
# load
# div_lst_mv_single <- readRDS("div_lst_mv_single.RDS")

div_df_mv_single <- bind_rows(mclapply(div_lst_mv_single, \(x){
  data.frame(
    frob = x$frob,
    inf = x$inf
  )
}), .id = "row")

# plot histogram
div_df_mv_single_long <- div_df_mv_single |>
  pivot_longer(c(frob, inf), names_to = "norm") |>
  mutate(norm = ifelse(norm == "frob", "Frobenius", "Infinity"))

div_df_mv_single_long |>
  ggplot(aes(x = value, fill = norm)) +
  geom_histogram(position = "dodge", bins = 30) +
  labs(
    # title = "Norms when moving 1 point (all sites, both models)",
    title = "Norms when moving 1 point (loc_1 only, both models)",
    x = "Norm value"
  ) +
  facet_wrap(~norm, scales = "free") +
  cecl_theme()

div_df_mv_single_long |>
  mutate(value = round(value, 3)) |>
  group_by(norm, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(norm) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  arrange(desc(pct))

summary(div_df_mv_single_long$value)
# ~66% are 0, and most others are very small

#### Divergence for same CE values (where single point swapped from 1 loc) ####

# But what if we only look at cases where each variable is equal to the mode?
# dep_coef_single_df_mode <- dep_coef_df_mv_single_spec |>
#   group_by(name, block, parameter) |>
#   filter(value == mode_fun(value)) |>
#   ungroup()
# # only keep those that are the same for all parameters
# dep_coef_single_df_mode <- dep_coef_single_df_mode |>
#   group_by(block, row) |>
#   filter(n_distinct(value) == 5) |>
#   ungroup()
#
# dep_lst_single_mode <- dep_coef_single_df_mode |>
#   mutate(
#     name = "loc_1",
#     var = "X1",
#     cond_var = "X2",
#     dth = dep_val,
#     parameter = ifelse(parameter == "LogLik", "ll", parameter)
#   ) |>
#   pivot_wider(names_from = parameter, values_from = value) |>
#   relocate(row, block, name, var, cond_var) |>
#   relocate(dth, .after = everything()) |>
#   identity() |>
#   group_split(row, .keep = FALSE) |>
#   lapply(\(x) {
#     group_split(x, block, .keep = FALSE) |>
#       lapply(\(y) {
#         as_cecl_dep(y)
#       })
#   })

# Only keep rows where loc 1 has the same params (other locs is the same)
mode_vals_df <- dep_coef_df_mv_single |>
  filter(name == "loc_1") |>
  group_by(block, name, var, cond_var) |>
  filter(
    a == mode_fun(a),
    b == mode_fun(b),
    m == mode_fun(m),
    s == mode_fun(s),
    ll == mode_fun(ll)
  ) |>
  ungroup()

dep_coef_single_df_mode <- dep_coef_df_mv_single |>
  filter(row %in% mode_vals_df$row)

dep_lst_single_mode <- dep_coef_single_df_mode |>
  group_split(row, .keep = FALSE) |>
  lapply(\(x) {
    group_split(x, block, .keep = FALSE) |>
      lapply(\(y) {
        as_cecl_dep(y)
      })
  })

mode_rows_single_df <- dep_coef_single_df_mode |>
  select(row, block) |>
  mutate(block = stringr::str_replace(block, "Block ", "")) |>
  mutate(across(everything(), as.numeric)) |>
  arrange(block, row)

mode_single_rows1 <- mode_rows_single_df |>
  filter(block == 1) |>
  pull(row) |>
  unique()
n_mode_single <- length(mode_single_rows1)


div_lst_single_mode <- mclapply(seq_along(dep_lst_single_mode), \(i) {
  # div_lst_single_mode <- lapply(seq_along(dep_lst_single_mode), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_mode_single)))
  dep_spec <- dep_lst_single_mode[[i]]
  data_loc_mv <- data_loc_test |>
    # mutate(block = ifelse(
    #   block == 1 & name == "loc_1" & row_number() == i, 2, block
    # ))
    mutate(block = ifelse(
      block == 1 & name == "loc_1" & row_number() == mode_single_rows1[[i]],
      2,
      block
    ))
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  dist <- lapply(seq_along(dep_spec), \(j) {
    cecl_dist(
      dep_spec[[j]],
      data_laplace_block[[j]],
      # n_mc = 500,
      # seed = seed
      laplace_sample = laplace_sample
    )
  })
  frob <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "F"
  )
  inf <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "I"
  )
  return(list(
    "frob" = frob,
    "inf" = inf
  ))
})

# extract norm values for each run, check if they're the same
div_df_single_mode <- bind_rows(mclapply(div_lst_single_mode, \(x){
  data.frame(
    # frob_dist1 = x$frob_dist1,
    # frob_dist2 = x$frob_dist2,
    # inf_dist1 = x$inf_dist1,
    # inf_dist2 = x$inf_dist2,
    frob = x$frob,
    inf = x$inf
  )
}), .id = "row")

# check that all are 0
div_df_mode |>
  select(-row) |>
  unlist() |>
  unique()
# They are also 0!
# Conclusion: In the case where we swap one row from all locations, it's much
# more unpredictable as to whether the CE parameters will be the same,
# but where we swap just one point, we know we get the same CE
# parameters 80% of the time, and also get 0 divergence as well, as
# desired

#### Divergence for ordered points ####

# For one location and one model, arrange data by size of conditioning
# variable and sense check that when we move non-extreme points we get the
# same CE fit (for the first ~80% of points)

# order points by conditioning variable
# data_loc_test_ord <- data_loc_test |>
#   filter(name == "loc_1") |>
#   arrange(block, X2) # will ensure points are swapped in order

data_loc_test_ord <- data_loc_test |>
  arrange(block, name, X2)

dep_lst_mv_single_ord <- mclapply(seq_len(n_per_loc), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  data_loc_mv <- data_loc_test_ord |>
    mutate(block = ifelse(
      # only swap points from first location
      block == 1 & name == "loc_1" & row_number() == i, 2, block
    ))
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )
  lapply(
    data_laplace_block,
    cecl_dep,
    cond_val = dep_val,
    fit_no_keef = TRUE,
    cond_var = "X2"
  )
})

# save
saveRDS(dep_lst_mv_single_ord, "dep_lst_mv_single_ord.RDS")
# load
# dep_lst_mv_single <- readRDS("dep_lst_mv_single.RDS")

dep_coef_df_mv_single_ord <- bind_rows(mclapply(dep_lst_mv_single_ord, \(x) {
  bind_rows(lapply(lapply(x, coef), as_tibble), .id = "block")
}, mc.cores = mc_cores), .id = "row") |>
  mutate(
    model = paste(var, cond_var, sep = " | "),
    block = paste("Block", block)
  )

# plot estimates vs row
dep_coef_df_mv_single_ord |>
  filter(name == "loc_1") |>
  pivot_longer(cols = c(a, b, m, s, ll), names_to = "parameter") |>
  mutate(parameter = ifelse(parameter == "ll", "LogLik", parameter)) |>
  ggplot(aes(x = as.numeric(row), y = value, color = block)) +
  geom_line() +
  geom_vline(xintercept = n_per_loc * 0.8, linetype = "dashed") +
  # facet_wrap(~parameter, scales = "free", ncol = 5) +
  facet_wrap(block ~ parameter, scales = "free", ncol = 5) +
  cecl_theme()

# Nice plot, shows how estimates are exactly the same for the first 80%
# of swaps, since they only swap non-extreme points in the conditioning
# variable, hence not effecting the CE fit.
# Interesting that for extreme swaps, the parameter estimates for each block
# look like mirror images of each other, which makes sense as the swapped point is
# extreme in opposite directions for each block, so would expect the CE fit to be
# affected in opposite directions for each block.

# also look at divergence
i <- 3
i <- 900
div_lst_mv_single_ord <- mclapply(seq_along(dep_lst_mv_single_ord), \(i) {
  system(sprintf("echo %s", paste0(i, "/", n_per_loc)))
  dep_spec <- dep_lst_mv_single_ord[[i]]
  data_loc_mv <- data_loc_test_ord |>
    mutate(block = ifelse(
      block == 1 & name == "loc_1" & row_number() == i, 2, block
    ))
  data_laplace_block <- trans_fun(
    data_loc_mv,
    n_vars = 2, laplace_trans = TRUE
  )

  # debugonce(cecl_dist)
  dist <- lapply(seq_along(dep_spec), \(j) {
    cecl_dist(
      dep_spec[[j]],
      data_laplace_block[[j]],
      # n_mc = 500,
      # seed = seed,
      laplace_sample = laplace_sample,
      cond_var = "X2"
    )
  })
  frob <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "F"
  )
  inf <- norm(
    as.matrix(dist[[2]]$dist_mat -
      dist[[1]]$dist_mat),
    type = "I"
  )
  return(list(
    "frob" = frob,
    "inf" = inf
  ))
})

# save
saveRDS(div_lst_mv_single_ord, "div_lst_mv_single_ord.RDS")
# load
# div_lst_mv_single_ord <- readRDS("div_lst_mv_single_ord.RDS")

div_df_mv_single_ord <- bind_rows(mclapply(div_lst_mv_single_ord, \(x){
  data.frame(
    frob = x$frob,
    inf = x$inf
  )
}), .id = "row")

# plot histogram
div_df_mv_single_ord_long <- div_df_mv_single_ord |>
  pivot_longer(c(frob, inf), names_to = "norm") |>
  mutate(
    row = as.numeric(row),
    norm = ifelse(norm == "frob", "Frobenius", "Infinity")
  )

div_df_mv_single_ord_long |>
  ggplot(aes(x = value, fill = norm)) +
  geom_histogram(position = "dodge", bins = 30) +
  labs(
    # title = "Norms when moving 1 point (all sites, both models)",
    title = "Norms when moving 1 point (loc_1 only, both models)",
    x = "Norm value"
  ) +
  facet_wrap(~norm, scales = "free") +
  cecl_theme()

# plot points
div_df_mv_single_ord_long |>
  # filter(row < n_per_loc * 0.8) |>
  ggplot(aes(x = row, y = value, color = norm)) +
  geom_point() +
  geom_vline(xintercept = n_per_loc * 0.8, linetype = "dashed") +
  facet_wrap(~norm, scales = "free") +
  cecl_theme()

# Any points before 80% mark that have non-zero divergence?
div_df_mv_single_ord_long |>
  group_by(norm) |>
  filter(value > mode_fun(value), row < n_per_loc * 0.8) |>
  ungroup() |>
  pivot_wider(names_from = norm, values_from = value) |>
  print(n = Inf)
# Only 4, reflecting 79.6% 0s


div_df_mv_single_ord_long |>
  mutate(value = round(value, 3)) |>
  group_by(norm, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(norm) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  arrange(desc(pct))
# ~80% are 0, and the rest are very small, as expected!

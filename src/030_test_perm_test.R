#### Testing permutation test method for determining changepoints ####

# TODO Look into whether Laplace transformation should be done at the start,
# rather than after permutating!

# Run permutation test under control case where Laplace sample and
# threshold is the same across locations and variables! (done)

# TODO Add example where data in blocks is identical (true control case)

# Want to understand sources of uncertainty in CE modelling and permutation test
# In particular, uncertainty due to:
# - Different dependence thresholds being used,
# - simulating from Laplace distribution for calculating divergence
# Want to determine if we can reduce the norm for the control case by
# controlling these uncertainties
# This may have to be done by increasing n_mc for distance calculation
# Or using the same threshold for all permutations

# For control case (same block length, just different run), frobenius and
# infinity norm differ by about 0.02, lets see if we can reduce this

# (done) Try with larger n_mc for distance calc (seems to help!)
# (done) Also try setting the seed before simulating from Laplace, or
# simulating outside the distance function and passing these values in
# (completely removes uncertainty!)
# (done) Allow single value for conditional threshold to be provided
# (done) Allow user-inputed Laplace cap, so that Laplace simulations will be
# the same when same seed and dependence threshold are used

# Currently, can have the same Laplace samples for each location (standard)
# and each variable, and also for each block (so for two different datasets).
# This means the only difference will be in the CE parameters estimated from
# data, and randomness due to threshold selection and Laplace sampling is removed.

# TODO Return & plot alpha, beta values for both setups

# (done) See if for same seed, Laplace sample is the same for all locations and
# variables, or different for all locations but the same for the control
# case

# Allow user inputted Laplace samples altogether! More intuitive than
# setting Laplace cap manually! (done)

# For each case, what do we want to be the same, and what do we want to
# differ?
# - Want threshold, Laplace simulation seed to be the same
# - Want CE parameters (esp. alpha and beta) to differ, as these are
#  estimated from data
# TODO Fill in result

# TODO Should I also have control case where data is the very same??

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(parallel)

#### metadata ####

seed <- 123
n_vars <- 2 # number of variables
n_locs <- 12 # number of locations
# n_locs <- 60
n_clust <- 3 # number of clusters
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per loc times n_locs per group
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group

# diff <- 0.4 # difference in correlation for local change
diff <- 0 # try for no difference, see what happens!
mc_cores <- detectCores() - 1 # number of cores for parallel processing

#### Functions ####

# Function to generate multivariate t data with specified correlation
# TODO Add option to transform to Laplace from GPD here?
gen_t <- \(cor_t, n, laplace_trans = FALSE) {
  # generate data
  cop_t <- tCopula(param = cor_t, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  # TODO Move GPD parameters to start of script
  ret <- data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))

  if (laplace_trans) {
    vars <- names(ret)
    gpd_vals <- data.frame(
      xi = -0.05,
      sigma = 1,
      thresh = 0
    ) |>
      tidyr::crossing(var = vars) |>
      mutate(var = factor(var, levels = vars), name = "dummy") |>
      mutate(name = "dummy") |>
      group_split(name, .keep = FALSE) |>
      lapply(\(x) group_split(x, var, .keep = FALSE))
    names(gpd_vals) <- "dummy"
    gpd_vals <- lapply(gpd_vals, `names<-`, vars)

    ret_lap <- trans_marg(
      marginal = gpd_vals,
      data_df = mutate(ret, name = "dummy"),
      mult_col = "name", vars = vars
    )
    return(as.data.frame(ret_lap[[1]]))
  }
  return(ret)
}

# Function to convert data to cecl_marg object
trans_fun <- \(data, n_vars, laplace_trans = FALSE) {
  # Calculate CE and cluster for all time blocks:

  # We know the exact GPD values, so use these to transform to Laplace
  locs <- uniqe(data$name)
  vars <- names(data)[1:n_vars]
  if (laplace_trans == FALSE) {
    gpd_vals <- data.frame(
      xi = -0.05,
      sigma = 1,
      thresh = 0,
      # change name and var to factors to ensure correct ordering when splitting
      name = factor(locs, levels = locs)
    ) |>
      tidyr::crossing(vars) |>
      rename(var = vars) |>
      mutate(var = factor(var, levels = vars)) |>
      group_split(name, .keep = FALSE) |>
      lapply(\(x) group_split(x, var, .keep = FALSE))

    gpd_vals <- lapply(gpd_vals, `names<-`, vars)
    names(gpd_vals) <- locs

    # perform transformation to standard Laplace margins (function from `CeCl`)
    # for each time block individually
    data_laplace_loc <- data |>
      mutate(block = factor(block, levels = unique(block))) |>
      group_split(block) |>
      lapply(\(x) {
        ret <- trans_marg(gpd_vals, x, vars = vars)
        names(ret) <- locs
        ret
      })
    # if data is already transformed, just split by block
  } else {
    # TODO Debug this, have the same as previous part!
    data_laplace_loc <- data |>
      group_split(block, .keep = FALSE)
  }
  names(data_laplace_loc) <- unique(data$block)

  # check you've transformed correctly
  # par(mfrow = c(1, 2))
  # plot(data_laplace_loc[[1]][[1]])
  # plot(data_laplace_loc[[2]][[1]])
  # par(mfrow = c(1, 1))

  # convert to cecl_marg objects
  marg_loc <- lapply(data_laplace_loc, as_cecl_marg)
}

# Function to fit CE, compute distances
calc_dist <- \(marg, n_mc = 500, cond_prob = 0.9, cond_val = NULL, ...) {
  if (!is.null(cond_val)) {
    cond_prob <- NULL
  }

  # fit CE models
  dep <- lapply(
    marg,
    cecl_dep,
    cond_prob   = cond_prob,
    cond_val    = cond_val,
    fit_no_keef = TRUE
  )

  # skip bad permutations
  if (any(sapply(dep, \(x) any(is.na(unlist(x$dependence)))))) {
    return(NA)
  }

  # compute distances
  dist <- lapply(seq_along(dep), \(i) {
    cecl_dist(
      dep[[i]],
      marg[[i]],
      # n_mc = 500
      n_mc = n_mc,
      ...
    )
  })
}

# Function to calculate norm or clustering diff
compare_blocks <- \(dist, type = c("norm", "clustering"), norm_type = "F") {
  if (type == "norm") {
    # norm of difference
    return(norm(
      as.matrix(dist[[2]]$dist_mat -
        dist[[1]]$dist_mat),
      type = norm_type
    ))
  } else if (type == "clustering") {
    # clustering solution
    clust <- lapply(
      dist,
      cecl_clust,
      k = n_clust,
      # cluster_mem = sort(rep(1:n_clust, each = n_locs / n_clust)),
      seed = seed
    )
    # calculate difference in clustering solutions
    return(mclust::adjustedRandIndex(
      clust[[1]]$pam$clustering,
      clust[[2]]$pam$clustering
    ))
  }
}

# Function to run single iteration of exploration for a given block position i
single_run_explore <- \(data_loc, i, laplace_trans = FALSE, ...) {
  data_block <- data_loc %>%
    group_by(name) |>
    mutate(block = ifelse(row_number() <= i, "1", "2")) |>
    ungroup()

  # check blocking is correct
  # table(data_block$block, data_block$name)

  # convert to Laplace and store marginals
  # TODO Should we do this before blocking???
  data_laplace_block <- trans_fun(data_block, n_vars, laplace_trans)

  # check you've transformed correctly
  # par(mfrow = c(1, 2))
  # plot(data_laplace_block[[1]][[1]][[1]])
  # plot(data_laplace_block[[2]][[1]][[1]])
  # par(mfrow = c(1, 1))

  dist <- suppressMessages(calc_dist(data_laplace_block, ...))

  # plot distance matrices
  # ggplot(dist[[1]], which = "image")
  # ggplot(dist[[2]], which = "image")

  norm_val_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
  norm_val_inf <- compare_blocks(dist, type = "norm", norm_type = "I")

  return(list(
    frob = norm_val_frob,
    inf = norm_val_inf
  ))
}


# plot norms
plot_norm <- \(block_vec, norm_vals_frob, norm_vals_inf) {
  data.frame(
    block_pos = block_vec,
    norm_frob = norm_vals_frob,
    norm_inf = norm_vals_inf
  ) |>
    tidyr::pivot_longer(
      cols = c(norm_frob, norm_inf),
      names_to = "norm_type",
      values_to = "norm_value"
    ) |>
    mutate(norm_type = ifelse(norm_type == "norm_frob", "Frobenius", "Infinity")) |>
    ggplot(aes(x = block_pos, y = norm_value, colour = norm_type)) +
    geom_line() +
    geom_point() +
    facet_wrap(~norm_type, scales = "free") +
    geom_vline(
      xintercept = (nrow(data_loc_orig) / n_locs) / 2,
      lty = 2, lwd = 2, colour = "red"
    ) +
    cecl_theme()
}

# function to find peaks in norm values and return corresponding block positions
find_peaks <- \(norm_vals, block_vec) {
  peaks <- c()
  for (j in 2:(length(norm_vals) - 1)) {
    if (norm_vals[j] > norm_vals[j - 1] && norm_vals[j] > norm_vals[j + 1]) {
      peaks <- c(peaks, block_vec[j])
    }
  }
  return(peaks)
}

# Function to test norm differences for various setups
test_norm_diff <- \(data_loc, block_vec, n_mc = 500, cond_val = 0, ...) {
  # norm_vals <- mclapply(block_vec, \(i) {
  norm_vals <- lapply(block_vec, \(i) {
    # single_run_explore(data_loc, i)
    # single_run_explore(data_loc, i, n_mc = 10000) # helps!
    # reduced by factor of 10!
    # single_run_explore(data_loc, i, n_mc = 10000, cond_val = dep_val)
    single_run_explore(data_loc, i, n_mc = n_mc, cond_val = cond_val, ...)
    # }, mc.cores = mc_cores)
  })

  # extract frobenius and infinity norms
  norm_vals_frob <- sapply(norm_vals, \(x) x$frob)
  norm_vals_inf <- sapply(norm_vals, \(x) x$inf)

  print(norm_vals_frob)
  print(norm_vals_inf)

  return(c("F" = abs(diff(norm_vals_frob)), "I" = abs(diff(norm_vals_inf))))
}


# # Function to perform a single permutation test
# Function to generate data for a single permutation test
perm_test_prep <- \(
  data_laplace_loc,
  cond_prob = 0.9,
  cond_val = NULL,
  # type = c("norm", "clustering"),
  # norm_type = "F"
  ...
) {
  # combine data into single data.frame with block indicator
  data_laplace_combined <- bind_rows(lapply(seq_along(data_laplace_loc), \(i) {
    bind_rows(lapply(data_laplace_loc[[i]]$transformed, as.data.frame), .id = "name") %>%
      mutate(block = i)
  })) |>
    select(X1, X2, name, block)

  # type <- match.arg(type, several.ok = FALSE)
  # permute rows within each location
  data_perm <- data_laplace_combined %>%
    group_by(name) %>%
    group_modify(\(df, key) {
      # original block sizes for this location
      n1 <- sum(df$block == levels(factor(df$block))[1])
      n2 <- sum(df$block == levels(factor(df$block))[2])

      # pool rows and randomly split
      idx <- sample.int(nrow(df))
      idx1 <- idx[seq_len(n1)]
      idx2 <- idx[-seq_len(n1)]

      rbind(
        df[idx1, ] %>% mutate(block = levels(factor(df$block))[1]),
        df[idx2, ] %>% mutate(block = levels(factor(df$block))[2])
      )
    }) %>%
    ungroup()

  # convert to cecl_marg objects (UNCHANGED)
  data_laplace_perm <- data_perm %>%
    mutate(block = factor(block, levels = unique(block))) %>%
    group_split(block) %>%
    lapply(\(x) {
      ret <- lapply(locs, \(loc) {
        subset(x, name == loc) %>%
          select(X1, X2)
      })
      names(ret) <- locs
      ret
    })
  names(data_laplace_perm) <- levels(data_perm$block)

  marg_perm <- lapply(data_laplace_perm, as_cecl_marg)

  # fit CE and compute distances
  dep_perm <- lapply(
    marg_perm,
    cecl_dep,
    cond_prob = cond_prob,
    cond_val = cond_val,
    fit_no_keef = TRUE
  )

  # skip bad permutations
  if (any(sapply(dep_perm, \(x) any(is.na(unlist(x$dependence)))))) {
    return(NA)
  }

  # compute distances
  dist_perm <- lapply(seq_along(dep_perm), \(i) {
    cecl_dist(
      dep_perm[[i]],
      marg_perm[[i]],
      n_mc = 500,
      ...
    )
  })
}
#
#   # fit CE, compute distances and calculate norm or clustering diff
#   compare_blocks(marg_perm, type, norm_type)
# }

# Function to process permutation test results
test_fun <- \(perm_test_res, norm_type = c("frob", "inf"), stat = c("mean", "median")) {
  norm_type <- match.arg(norm_type, several.ok = FALSE)
  stat <- match.arg(stat, several.ok = FALSE)

  sapply(perm_test_res, \(res) {
    if (norm_type == "frob") {
      norm_orig <- res$norm_orig_frob
      perm_norms <- res$perm_norms_frob
    } else if (norm_type == "inf") {
      norm_orig <- res$norm_orig_inf
      perm_norms <- res$perm_norms_inf
    }

    if (stat == "mean") {
      return(abs(mean(perm_norms) - norm_orig))
    } else if (stat == "median") {
      return(abs(median(perm_norms) - norm_orig))
    }
  })
}


#### Generate data ####

# increase lowest correlation by 0.05 until it reaches middle correlation
set.seed(seed)
# start from 0 to have control example (checks difference due to seed)
# increments <- c(0, 0:8)
# increments <- c(0, 0:1)
increments <- c(0, 1) # TODO Add more increments once code works
data_loc <- bind_rows(lapply(seq_along(increments), \(i) {
  cor_spec <- cor_t_vec
  cor_spec[[1]] <- cor_spec[[1]] + (diff * increments[[i]])
  print(cor_spec)
  # generate data
  # ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n))
  # debugonce(gen_t)
  ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n, laplace_trans = TRUE))
  ret |>
    mutate(
      name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
      block = i # time block
    )
}))

# check data is correct
table(data_loc$block)
table(data_loc$name)

# plot
data_loc |>
  # should have increasing cor
  filter(name == "loc_1") |>
  ggplot() +
  geom_point(aes(X1, X2)) +
  facet_wrap(~block) +
  cecl_theme()

data_loc |>
  filter(name == "loc_1") |>
  group_by(block) |>
  summarise(cor = cor(X1, X2), .groups = "drop")
# ggplot() +
# cecl_theme()

# save original data and remove block info
data_loc_orig <- data_loc
data_loc <- select(data_loc, -block)

# We know the exact GPD values, so use these to transform to Laplace
locs <- unique(data_loc$name)
vars <- names(data_loc)[1:n_vars]

#### Compute norms for control case (where blocks are exact same) ####

# size of blocks will be 1000 in both cases
block_vec <- c(1000, 1000)

# TODO Move to it's own code block
# use the same dependence threshold across all variables and locations
(dep_val <- qgpd(0.9, u = 0, xi = -0.05, sigma = 1)) # use 90th quantile of GPD
# similar to 90th quantile of data, unsurprisingly
data_loc_orig |>
  filter(name == "loc_1") |>
  group_by(block) |>
  summarise(
    q90 = quantile(X1, 0.9),
    .groups = "drop"
  )

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

laplace_sample <- rlaplace_trunc(
  n = 500,
  thresh_max = dep_val,
  y_max = laplace_cap
)

# use the same data in both blocks
data_loc_test <- data_loc_orig |>
  filter(block == 1) |>
  bind_rows(
    data_loc_orig |>
      filter(block == 1) |>
      mutate(block = 2)
  )

# check that data is the exact same
data_loc_test |>
  # should have increasing cor
  filter(name == "loc_1") |>
  select(X1, block) |>
  # convert to wide format
  group_split(block, .keep = FALSE) |>
  bind_cols() |>
  setNames(c("X11", "X12")) |>
  mutate(diff = X11 - X12) |>
  summarise(max(abs(diff)))

data_loc_test |>
  # should have increasing cor
  filter(name == "loc_1") |>
  ggplot() +
  geom_point(aes(X1, X2)) +
  facet_wrap(~block) +
  cecl_theme()

# should be slightly different for different Laplace samples and thresholds
# debugonce(single_run_explore)
debugonce(trans_fun)
test_norm_diff(data_loc_test, block_vec, cond_val = NULL, laplace_trans = TRUE) # almost 0, good
# should be less different for same threshold
test_norm_diff(data_loc_test, block_vec, cond_val = dep_val)
# should be identical for same Laplace sample and threshold
test_norm_diff(
  data_loc_test, block_vec,
  cond_val = dep_val, laplace_sample = laplace_sample
)
# successful!

# alpha, beta estimates should also be identical!
data_block_test <- data_loc_test %>%
  group_by(name) |>
  mutate(block = ifelse(row_number() <= block_vec[1], "1", "2")) |>
  ungroup()
debugonce(trans_fun)
data_laplace_block_test <- trans_fun(data_block_test, n_vars = 2)

dep <- lapply(
  data_laplace_block_test,
  cecl_dep,
  cond_prob   = NULL,
  cond_val    = dep_val,
  fit_no_keef = TRUE
)

# coefficients are the same!! as desired :)
coef(dep[[1]])[, c("a", "b", "m", "s", "dth")] -
  coef(dep[[2]])[, c("a", "b", "m", "s", "dth")]


#### Compute norms for control case (data re-simulated with same corr) ####

# vector of block positions to explore (i.e., size of first block)
# block_vec <- seq(min_size, max_size, by = 100)
# block_vec <- c(1000, 1001)
block_vec <- c(1000, 1000) # control setting


# run exploration across all block positions
# TODO Debug, see what sources of uncertainty can be controlled!
set.seed(seed)
# debugonce(single_run_explore)
# debugonce(calc_dist)


# run once for various setups
test_norm_diff(data_loc, block_vec, n_mc = 500, cond_val = NULL)
test_norm_diff(data_loc, block_vec, n_mc = 10000, cond_val = NULL)
test_norm_diff(data_loc, block_vec, n_mc = 10000, cond_val = dep_val)
test_norm_diff(data_loc, block_vec, n_mc = 10000, cond_val = dep_val, seed = 123)
test_norm_diff(data_loc, block_vec, n_mc = 500, seed = 123) # same seed removes variability
# same with using same laplace samples for all (not just initial seed)
test_norm_diff(data_loc, block_vec, cond_val = dep_val, laplace_sample = laplace_sample)

# conclusion: Using same seed/laplace samples completely removes variability in this case

# Check that estimated CE parameters are different, even for same block size
data_block <- data_loc %>%
  group_by(name) |>
  mutate(block = ifelse(row_number() <= block_vec[1], "1", "2")) |>
  ungroup()
data_laplace_block <- trans_fun(data_block, n_vars = 2)

dep <- lapply(
  data_laplace_block,
  cecl_dep,
  cond_prob   = NULL,
  cond_val    = dep_val,
  fit_no_keef = TRUE
)

# coefficients different (but thresholds the same), as desired
coef(dep[[1]])[, c("a", "b", "m", "s", "dth")] -
  coef(dep[[2]])[, c("a", "b", "m", "s", "dth")]

# TODO Seem extremely different though, why??
ggplot(dep[[1]], which = "scatter", var = "X1", cond_var = "X2")
ggplot(dep[[2]], which = "scatter", var = "X1", cond_var = "X2")

# also check that ???
dist <- suppressMessages(calc_dist(
  data_laplace_block,
  cond_val = dep_val,
  laplace_sample = laplace_sample, # norms are lower when using same samples
  seed = 123
))
ggplot(dist[[1]], which = "image")
ggplot(dist[[2]], which = "image")

dist_diff$dist_mat <- abs(dist[[2]]$dist_mat - dist[[1]]$dist_mat)
image(as.matrix(dist_diff$dist_mat))

# calculate norms
# TODO These values are very large, why for such similar data???
(norm_val_frob <- compare_blocks(dist, type = "norm", norm_type = "F"))
(norm_val_inf <- compare_blocks(dist, type = "norm", norm_type = "I"))



#### Compute norms for time 1000 and 1001 ####

# This time, we hope to see a small difference in norms, as there is a change in
# correlation structure between these two time blocks
block_vec1 <- c(1000, 1001)

test_norm_diff(block_vec1, seed = 123) # some variability, as desired!
test_norm_diff(block_vec1, seed = 123, cond_val = dep_val)

test1 <- replicate(100, test_norm_diff(block_vec1, seed = 123))
test2 <- replicate(100, test_norm_diff(block_vec1, seed = 123, cond_val = dep_val))

# check results
list(test1 = test1, test2 = test2) |>
  sapply(\(x) apply(x, 1, mean))


#### Test permutation test ####

# grid_vals <- seq(950, 1050, by = 10) # smaller grid to test
# grid_vals <- c(1000, 1001)
# grid_vals <- c(1000, 1000)
grid_vals <- 1000
n_perm <- 100 # number of permutations

set.seed(seed)
perm_test_res <- lapply(grid_vals, \(i) {
  # print how many iterations completed (as a percentage)
  print(paste0(substr(which(grid_vals == i) / length(grid_vals) * 100, 0, 5), "% complete\n"))

  # Add block based on i
  data_block <- data_loc %>%
    group_by(name) |>
    mutate(block = ifelse(row_number() <= i, "1", "2")) |>
    ungroup()

  # check blocking is correct
  # table(data_block$block, data_block$name)

  # convert to Laplace and store marginals
  # TODO Should this be done once a the start of this script???
  data_laplace_block <- trans_fun(data_block, n_vars)

  # compute permutation distances for original data
  norm_vals <- lapply(seq_len(n_perm), \(j) {
    dist <- perm_test_prep(
      data_laplace_block,
      laplace_sample = laplace_sample,
      cond_val = dep_val,
      cond_prob = NULL
    )
    list(
      "frob" = compare_blocks(dist, type = "norm", norm_type = "F"),
      "inf"  = compare_blocks(dist, type = "norm", norm_type = "I")
    )
  })

  # extract frobenius and infinity norms
  perm_norms_frob <- sapply(norm_vals, \(x) x$frob)
  perm_norms_inf <- sapply(norm_vals, \(x) x$inf)

  # compute distances for original data
  dist <- suppressMessages(calc_dist(data_laplace_block))
  norm_orig_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
  norm_orig_inf <- compare_blocks(dist, type = "norm", norm_type = "I")

  # compute p-values
  # TODO Are these p-values for one-sided or two-sided tests?
  (p_value_frob <- mean(perm_norms_frob >= norm_orig_frob))
  (p_value_inf <- mean(perm_norms_inf >= norm_orig_inf))

  return(list(
    p_value_frob    = p_value_frob,
    p_value_inf     = p_value_inf,
    norm_orig_frob  = norm_orig_frob,
    norm_orig_inf   = norm_orig_inf,
    perm_norms_frob = perm_norms_frob,
    perm_norms_inf  = perm_norms_inf
  ))
})

hist(
  # perm_test_res[[1]]$perm_norms_inf,
  perm_test_res[[1]]$perm_norms_frob,
  xlim = c(
    0,
    # max(c(perm_test_res[[1]]$perm_norms_inf, perm_test_res[[1]]$norm_orig_inf))
    max(c(perm_test_res[[1]]$perm_norms_frob, perm_test_res[[1]]$norm_orig_frob))
  )
)
abline(v = perm_test_res[[1]]$norm_orig_inf, col = "red", lwd = 2)

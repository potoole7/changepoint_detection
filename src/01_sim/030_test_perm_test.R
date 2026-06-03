#### Testing permutation test method for determining changepoints ####

# TODO dep_val should be 90th quantile of Laplace, not GPD!!
# Change and investigate effect

# TODO Make sense of all of this !!!

# TODO Should I run permutation test on data that has no difference?
# Or have I done this already?!?

# TODO CE extremes fits vary quite a bit when permutating; perhaps
# fitting an alternative (such as with evgam) would be more stable?
# TODO Try this out!! Would need to generate spatial data

# Main problems:
# - Difference in norm values for neighbouring blocks can be large
# - Norm values even under no difference can be very large
#   - Good news is that this seems to be because of different CE fits,
#     rather than uncertainty introduced by different Laplace samples etc (in our method)

# repeat for more data to hopefully reduce norm values! (done)

# TODO Do sliding test with data with literally no difference! See if
# it improves somewhat in stability

# TODO Read notes from Christian meeting!

# Run permutation test many times to see what differences there are between configurations
# - Standard,
# - Pre-calculated dependence quantile
# - Set seed w or w/o pre-calculated dependence quantile
# - Pre-calculated Laplace sample
# See why/how these make a difference!
# - Create function to do permutation test (done)
# - Allow function to take pre-calculated dependence quantile,
# seed number and pre-calculated Laplace sample (done)
# - Run! (done)
# - Investigate possible bugs in permutating? Norm values seem odd for set seeds

# TODO Are p-values for one-sided or two-sided tests??

# TODO Also run multiple times with different pre-calculated DQU values

# TODO Important: See if setting dependence quantile the same and/or using
# pre-specified Laplace sample is working correctly, as it seems to
# produce larger norm values than just using the same seed across all simulations
# but letting the quantile and Laplace sample vary
# Setting threshold: Is successful in using same dependence threshold (i.e.
# code works), but why does it increase the size of differences?

# Differences for control case with identical data are very smaller
# (norms ~ 0.005), so at least that control makes sense!
# However, norms are much larger (~0.4-0.6) when doing a simple control case
# where data is re-sampled from Laplace distribution with no difference in
# marginals or dependence correlation!
# Note that this difference is much lower when using the same seed
# (reduced to 0.~13), but interestingly not for using the same Laplace sample
# (increases to ~0.65)!
# This may be worrying, right?
# Or does the idea of a significant difference just change then?

# Look into whether Laplace transformation should be done at the start,
# rather than after permutating! (done)

# TODO Investigate why setting seed and *not* using a set dependence threshold
# results in lowest norm values
# TODO Also look back into difference in norm values for different block sizes
# which shouldn't be too large!
# TODO Write some conclusions on this!! (and moving permutation test script)
# TODO Look into differences in CE parameter values which may result in having
# quite large norm values even for very similar data
# TODO Read notes from Christian meeting

# Run permutation test under control case where Laplace sample and
# threshold is the same across locations and variables! (done)

# Add example where data in blocks is identical (true control case) (done)

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
library(tidyr)
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
n_per_block <- n * n_clust # number of obs per block (two blocks here)
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group

# diff <- 0.4 # difference in correlation for local change
diff <- 0 # try for no difference, see what happens!
mc_cores <- detectCores() - 1 # number of cores for parallel processing

#### Functions ####

# Function to generate multivariate t data with specified correlation
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
  locs <- unique(data$name)
  vars <- names(data)[1:n_vars]
  # If not previously transformed to Laplace margins, do so now by locatoin
  # TODO Confusing/counter-intuitive argument name?
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
        # transforms each location individually
        ret <- trans_marg(gpd_vals, x, vars = vars)
        names(ret) <- locs
        ret
      })
    # if data is already transformed, just split by block
  } else {
    # data already
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
  # TODO Do before blocking? Think I do with laplace_trans = TRUE!
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
# Idea is that for similar/same block sizes, differences shouldn't be large
test_norm_diff <- \(data_loc, block_vec, n_mc = 500, cond_val = 0, ret_all = TRUE, ...) {
  # norm_vals <- mclapply(block_vec, \(i) {
  norm_vals <- lapply(block_vec, \(i) {
    single_run_explore(data_loc, i, n_mc = n_mc, cond_val = cond_val, ...)
    # }, mc.cores = mc_cores)
  })

  # extract frobenius and infinity norms
  norm_vals_frob <- sapply(norm_vals, \(x) x$frob)
  norm_vals_inf <- sapply(norm_vals, \(x) x$inf)

  # print values of norms
  print(norm_vals_frob)
  print(norm_vals_inf)

  # return differences (and also norms if desired)
  diff <- c("F" = abs(diff(norm_vals_frob)), "I" = abs(diff(norm_vals_inf)))
  if (ret_all) {
    return(list("diff" = diff, "norm_vals" = norm_vals))
  } else {
    return(diff)
  }
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
  if (!is.null(cond_val)) {
    cond_prob <- NULL
  }

  # combine data into single data.frame with block indicator
  data_laplace_combined <- bind_rows(lapply(seq_along(data_laplace_loc), \(i) {
    bind_rows(lapply(data_laplace_loc[[i]]$transformed, as.data.frame), .id = "name") %>%
      mutate(block = i)
  })) |>
    select(X1, X2, name, block)

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

      df |>
        mutate(
          row = row_number(),
          block = case_when(
            row %in% idx1 ~ levels(factor(df$block))[1],
            row %in% idx2 ~ levels(factor(df$block))[2]
          )
        )
    }) %>%
    ungroup() |>
    select(-row)

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

  # fit CE
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

# Function to perform permutation test for a given block position i
perm_test_fun <- \(data_loc, grid_vals, n_perm = 100, laplace_trans = TRUE, ...) {
  lapply(grid_vals, \(i) {
    # print how many iterations completed (as a percentage)
    print(paste0(
      substr(which(grid_vals == i) / length(grid_vals) * 100, 0, 5), "% complete\n"
    ))

    # Add block based on i
    data_block <- data_loc %>%
      group_by(name) |>
      mutate(block = ifelse(row_number() <= i, "1", "2")) |>
      ungroup()

    # check blocking is correct
    # table(data_block$block, data_block$name)

    # convert to Laplace and store marginals
    data_laplace_block <- trans_fun(data_block, n_vars, laplace_trans)

    # compute permutation distances for original data
    norm_vals <- lapply(seq_len(n_perm), \(j) {
      # norm_vals <- mclapply(seq_len(n_perm), \(j) {
      dist <- perm_test_prep(
        data_laplace_block,
        # laplace_sample = laplace_sample,
        ...
      )
      list(
        "frob" = compare_blocks(dist, type = "norm", norm_type = "F"),
        "inf"  = compare_blocks(dist, type = "norm", norm_type = "I")
      )
    })
    # }, mc.cores = mc_cores)

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

#### Precalculations ####

# TODO Move to it's own code block
# use the same dependence threshold across all variables and locations
(dep_val <- qgpd(0.8, u = 0, xi = -0.05, sigma = 1)) #
# (dep_val <- qgpd(0.9, u = 0, xi = -0.05, sigma = 1)) # use 90th quantile of GPD

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

set.seed(seed)
laplace_sample <- rlaplace_trunc(
  n = 500, # TODO Increase later
  thresh_max = dep_val,
  y_max = laplace_cap
)

# Also make a longer Laplace sample
laplace_sample2 <- rlaplace_trunc(
  n = 2000,
  thresh_max = dep_val,
  y_max = laplace_cap
)

#### Compute norms for control case (where blocks are exact same) ####

dep_orig <- dep # store original dependence objects

# size of blocks will be 1000 in both cases
block_vec <- c(n_per_loc, n_per_loc)

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
set.seed(seed)
test_norm_diff(data_loc_test, block_vec, cond_val = NULL, laplace_trans = TRUE)
# [1] 0.005217737 0.024858710 # F norm values
# [1] 0.005333177 0.028638695 # I norm values
# F          I
# 0.01964097 0.02330552
# almost 0, good
# TODO should be less different (and smaller) for same threshold, why not?!
# May be because dep_val is larger than the threshold selected by letting
# it vary?
test_norm_diff(
  data_loc_test, block_vec,
  cond_val = dep_val,
  # cond_val = qgpd(0.8, u = 0, xi = -0.05, sigma = 1),
  laplace_trans = TRUE
)
# 0.04724554 0.05732726
# [1] 0.06477626 0.00734077
# [1] 0.074324763 0.008784031
# F          I
# 0.05743549 0.06554073
# using same seed gives 0 norms (i.e. succesfully found to be identical)
# TODO Fill in comment below
# Therefore, ...
# debugonce(cecl_dist)
test_norm_diff(data_loc_test, block_vec, laplace_trans = TRUE, seed = seed)
# [1] 0 0
# [1] 0 0
# F I
# 0 0
# Also identical for same Laplace sample and threshold
test_norm_diff(
  data_loc_test,
  block_vec,
  cond_val       = dep_val,
  laplace_sample = laplace_sample,
  laplace_trans  = TRUE
)
# [1] 0 0
# [1] 0 0
# F I
# 0 0

# Run these tests many times to better identify patterns
test1 <- replicate(100, test_norm_diff(data_loc_test, block_vec, cond_val = NULL, laplace_trans = TRUE, ret_all = TRUE))
test2 <- replicate(100, test_norm_diff(
  data_loc_test,
  block_vec,
  cond_val = dep_val,
  laplace_trans = TRUE,
  ret_all = TRUE
))
test3 <- replicate(100, test_norm_diff(data_loc_test, block_vec, laplace_trans = TRUE, seed = seed, ret_all = TRUE))
test4 <- replicate(100, test_norm_diff(
  data_loc_test,
  block_vec,
  cond_val       = dep_val,
  laplace_sample = laplace_sample,
  laplace_trans  = TRUE,
  ret_all        = TRUE
))

# extract norm values and differences, compare for different setups
norm_df <- bind_rows(lapply(list(test1, test2, test3, test4), \(x) {
  idx <- seq(2, length(test1), by = 2)
  vals <- unlist(x[idx])
  data.frame(
    "frob" = vals[names(vals) == "frob"],
    "inf"  = vals[names(vals) == "inf"]
  )
}), .id = "indicator")

norm_df |>
  group_by(indicator) |>
  mutate(row = row_number()) |>
  ungroup() |>
  ggplot(aes(x = row)) +
  geom_line(aes(y = frob), colour = "red") +
  geom_line(aes(y = inf), colour = "blue") +
  facet_wrap(~indicator)

norm_df |>
  group_by(indicator) |>
  summarise(across(c("frob", "inf"), \(x) median(x)), .groups = "drop")

# check differences

# # test changing fixed dependence threshold value
# # TODO Would have to run many times to check how this effects things!
# debugonce(cecl_dist)
dep_vals <- lapply(c(0.8, 0.85, 0.9), \(x) {
  qgpd(x, u = 0, xi = -0.05, sigma = 1)
})

set.seed(seed)
test_dep <- replicate(100, lapply(dep_vals, \(x) {
  test_norm_diff(data_loc_test, block_vec, cond_val = x, laplace_trans = TRUE)
}))

# alpha, beta estimates should also be identical!
data_block_test <- data_loc_test %>%
  group_by(name) |>
  mutate(block = ifelse(row_number() <= block_vec[1], "1", "2")) |>
  ungroup()
data_laplace_block_test <- trans_fun(data_block_test, n_vars = 2, laplace_trans = TRUE)

dep <- lapply(
  data_laplace_block_test,
  cecl_dep,
  cond_prob   = NULL,
  cond_val    = dep_val,
  fit_no_keef = TRUE
)

# coefficients are the same!! as desired :)
all.equal(coef(dep[[1]]), coef(dep[[2]]))
# Therefore, non-zero norms for identical data is due to variability in
# Laplace sample being used


#### Permutation tests for various setups ####

# grid_vals <- 1000
grid_vals <- n_per_loc #
n_perm <- 100 # number of permutations
set.seed(seed)

# TODO ALso add n_mc argument? 500 vs (how many?) more? Would be NA if
# using laplace_sample! (also use two different laplace samples)

run_df <- tidyr::crossing(
  grid_vals      = grid_vals,
  n_perm         = n_perm,
  dep_val        = c(NA, dep_val), # can't use NULL here, can sub in later
  seed           = c(NA, seed),
  laplace_sample = c(NA, paste(laplace_sample, collapse = "---")) # to allow crossing
) |>
  # can't have both seed and laplace_sample
  filter(
    !(!is.na(seed) & !is.na(laplace_sample))
  )

# # run for single setup to test
# set.seed(123)
# perm_test_fun(
#   data_loc       = data_loc_test,
#   grid_vals      = grid_vals,
#   n_perm         = n_perm,
#   # cond_val       = NULL,
#   cond_val       = dep_val,
#   # seed           = NULL,
#   seed           = 123,
#   laplace_sample = NULL
#   # laplace_sample = laplace_sample
# )

set.seed(seed)
# break run_df into list, run each row separately
run_df_fun <- \(data_loc_test, run_df) {
  perm_test_res_lst <- run_df |>
    group_split(row_number(), .keep = FALSE) |>
    lapply(function(x) {
      grid_vals <- x[["grid_vals"]]
      n_perm <- x[["n_perm"]]
      n_mc <- 500 # default
      if (!is.null(x[["n_mc"]])) {
        n_mc <- x[["n_mc"]]
      }

      cond_val <- if (is.na(x[["dep_val"]])) NULL else as.numeric(x[["dep_val"]])
      seed <- if (is.na(x[["seed"]])) NULL else as.numeric(x[["seed"]])

      laplace_sample <- x[["laplace_sample"]]
      laplace_sample <- if (is.na(laplace_sample)) {
        NULL
      } else {
        as.numeric(strsplit(laplace_sample, "---", fixed = TRUE)[[1]])
      }

      print(paste("grid_vals:", grid_vals))
      print(paste("n_perm:", n_perm))
      print(paste("n_mc:", n_mc))
      print(paste("cond_val:", ifelse(is.null(cond_val), "NULL", cond_val)))
      print(paste("seed:", ifelse(is.null(seed), "NULL", seed)))
      print(paste("laplace_sample:", ifelse(is.null(laplace_sample), "NULL", "Provided")))

      res <- perm_test_fun(
        data_loc       = data_loc_test,
        grid_vals      = grid_vals,
        n_perm         = n_perm,
        cond_val       = cond_val,
        seed           = seed,
        laplace_sample = laplace_sample,
        n_mc           = n_mc
      )
      print("Ran!")
      return(res)
    })
}

perm_test_res_lst <- run_df_fun(data_loc_test, run_df)

# extract norm values
norm_df <- bind_rows(lapply(seq_len(nrow(run_df)), \(i) {
  res <- perm_test_res_lst[[i]][[1]]
  data.frame(
    row = i,
    grid_vals = run_df$grid_vals[i],
    n_perm = run_df$n_perm[i],
    cond_val = ifelse(is.na(run_df$dep_val[i]), NA, run_df$dep_val[i]),
    seed = ifelse(is.na(run_df$seed[i]), NA, run_df$seed[i]),
    laplace_sample = ifelse(is.na(run_df$laplace_sample[i]), "NULL", "Provided"),
    frob = res$perm_norms_frob,
    inf = res$perm_norms_inf
  )
}))

# plot
norm_df_long <- norm_df |>
  tidyr::pivot_longer(
    cols = c(frob, inf),
    names_to = "norm_type",
    values_to = "norm_value"
  ) |>
  mutate(norm_type = ifelse(norm_type == "frob", "Frobenius", "Infinity"))

norm_df_long |>
  ggplot(aes(x = norm_value)) +
  geom_histogram(bins = 30, fill = "lightblue", colour = "black") +
  facet_grid(
    norm_type ~ row,
    scales = "fixed"
  ) +
  scale_x_continuous(limits = c(0, 1.5)) +
  cecl_theme()

# Find lowest mean norm values by setting
norm_df_long |>
  group_by(row, norm_type) |>
  summarise(
    mean_norm = mean(norm_value),
    min_norm = min(norm_value),
    max_norm = max(norm_value),
    .groups = "drop"
  ) |>
  arrange(mean_norm) |>
  # tidyr::pivot_wider(names_from = norm_type, values_from = mean_norm) |>
  identity()

# A tibble: 12 × 5
#     row norm_type mean_norm min_norm max_norm
#       <int> <chr>         <dbl>    <dbl>    <dbl>
# 1     3 Frobenius     0.556    0.267    0.957
# 2     1 Frobenius     0.567    0.261    0.950
# 3     2 Frobenius     0.569    0.324    1.12
# 4     5 Frobenius     0.574    0.251    1.07
# 5     6 Frobenius     0.606    0.251    1.09
# 6     4 Frobenius     0.663    0.271    1.20
# 7     3 Infinity      0.679    0.280    1.44
# 8     1 Infinity      0.685    0.280    1.38
# 9     5 Infinity      0.705    0.290    1.60
# 10    2 Infinity      0.717    0.351    1.95
# 11    6 Infinity      0.740    0.279    1.46
# 12    4 Infinity      0.874    0.399    1.86
# means around 0.5/0.6 for Frob and 0.6/0.7 for Inf, with min values around 0.2/0.3

#### Repeat for larger sample ####

# TODO See if increasing sample size changes results
# TODO See if increasing n_perm changes results
# TODO See if changing n_mc changes results
# n_per_loc2 <- 2000 # increased from 1000
# n_per_loc2 <- 5000
n_per_loc2 <- 10000
n2 <- n_per_loc2 * (n_locs / n_clust)
# n_perm2 <- 100
n_perm2 <- 500 # increase number of permutations too
# n_mc <- 500
n_mc2 <- 2000 # increase number of MC samples to use!

set.seed(seed)
increments <- c(0, 1)
data_loc2 <- bind_rows(lapply(seq_along(increments), \(i) {
  cor_spec <- cor_t_vec
  cor_spec[[1]] <- cor_spec[[1]] + (diff * increments[[i]])
  print(cor_spec)
  # generate data
  ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n2, laplace_trans = TRUE))
  ret |>
    mutate(
      name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
      block = i # time block
    )
}))

block_vec2 <- c(n_per_loc2, n_per_loc2)

# use the same data in both blocks
data_loc_test2 <- data_loc2 |>
  filter(block == 1) |>
  bind_rows(
    data_loc_orig |>
      filter(block == 1) |>
      mutate(block = 2)
  )

# change setup from initial permutation tests
run_df2 <- run_df |>
  mutate(
    grid_vals = n_per_loc2, # middle grid value now has to change
    n_perm = n_perm2, # number of permutations increased
    n_mc = n_mc2, # Monte Carlo samples increases
    # Use pre-calculated Laplace sample of size 2000
    laplace_sample = case_when(
      !is.na(.data[["laplace_sample"]]) ~ paste(laplace_sample2, collapse = "---"),
      TRUE ~ laplace_sample
    )
  ) |>
  relocate(n_mc, .after = n_perm)

perm_test_res_lst2 <- run_df_fun(data_loc_test2, run_df2)

# extract norm values
norm_df2 <- bind_rows(lapply(seq_len(nrow(run_df2)), \(i) {
  res <- perm_test_res_lst2[[i]][[1]]
  data.frame(
    row = i,
    grid_vals = run_df2$grid_vals[i],
    n_perm = run_df2$n_perm[i],
    n_mc = run_df2$n_mc[i],
    cond_val = ifelse(is.na(run_df2$dep_val[i]), NA, run_df2$dep_val[i]),
    seed = ifelse(is.na(run_df2$seed[i]), NA, run_df2$seed[i]),
    laplace_sample = ifelse(is.na(run_df2$laplace_sample[i]), "NULL", "Provided"),
    frob = res$perm_norms_frob,
    inf = res$perm_norms_inf
  )
}))

# plot
norm_df_long2 <- norm_df2 |>
  tidyr::pivot_longer(
    cols = c(frob, inf),
    names_to = "norm_type",
    values_to = "norm_value"
  ) |>
  mutate(norm_type = ifelse(norm_type == "frob", "Frobenius", "Infinity"))

norm_df_long2 |>
  ggplot(aes(x = norm_value)) +
  geom_histogram(bins = 30, fill = "lightblue", colour = "black") +
  facet_grid(
    norm_type ~ row,
    scales = "fixed"
  ) +
  scale_x_continuous(limits = c(0, 1.5)) +
  cecl_theme()

# Find lowest mean norm values by setting
norm_df_long2 |>
  group_by(row, norm_type) |>
  summarise(
    mean_norm = mean(norm_value),
    min_norm = min(norm_value),
    max_norm = max(norm_value),
    .groups = "drop"
  ) |>
  arrange(mean_norm) |>
  identity()

# for 2000 samples:
#     row norm_type mean_norm min_norm max_norm
#   <int> <chr>         <dbl>    <dbl>    <dbl>
# 1     2 Frobenius     0.483    0.219    0.832
# 2     1 Frobenius     0.484    0.240    0.854
# 3     3 Frobenius     0.485    0.203    0.734
# 4     5 Frobenius     0.499    0.214    0.902
# 5     4 Frobenius     0.566    0.276    1.05
# 6     6 Frobenius     0.572    0.215    1.10
# 7     1 Infinity      0.581    0.230    1.37
# 8     2 Infinity      0.586    0.239    1.27
# 9     3 Infinity      0.594    0.236    1.06
# 10    5 Infinity      0.618    0.255    1.33
# 11    6 Infinity      0.690    0.270    1.43
# 12    4 Infinity      0.701    0.308    1.61
# for 5000 samples:
#      row norm_type mean_norm min_norm max_norm
#     <int>  <chr>         <dbl>    <dbl>    <dbl>
# 1      1 Frobenius     0.433    0.222    0.644
# 2      3 Frobenius     0.441    0.229    0.756
# 3      2 Frobenius     0.459    0.259    0.841
# 4      5 Frobenius     0.499    0.279    1.05
# 5      4 Frobenius     0.532    0.226    1.08
# 6      6 Frobenius     0.535    0.291    0.874
# 7      3 Infinity      0.550    0.288    1.25
# 8      1 Infinity      0.551    0.257    1.00
# 9      2 Infinity      0.566    0.260    1.09
# 10     5 Infinity      0.645    0.278    1.91
# 11     4 Infinity      0.675    0.228    1.33
# 12     6 Infinity      0.687    0.308    1.37
# for 10000 samples:
#      row norm_type mean_norm min_norm max_norm
#    <int> <chr>         <dbl>    <dbl>    <dbl>
# 1      5 Frobenius     0.410    0.199    0.700
# 2      2 Frobenius     0.427    0.215    0.803
# 3      1 Frobenius     0.431    0.216    0.769
# 4      3 Frobenius     0.448    0.231    0.866
# 5      5 Infinity      0.508    0.224    0.967
# 6      4 Frobenius     0.510    0.240    0.996
# 7      6 Frobenius     0.511    0.211    1.16
# 8      1 Infinity      0.544    0.230    1.33
# 9      2 Infinity      0.547    0.200    1.45
# 10     3 Infinity      0.575    0.244    1.20
# 11     4 Infinity      0.650    0.247    1.62
# 12     6 Infinity      0.658    0.235    1.60
# For 10000 samples and 500 permutations for the permutation test
#     row norm_type mean_norm min_norm max_norm
#   <int> <chr>         <dbl>    <dbl>    <dbl>
# 1     1 Frobenius     0.419    0.180    0.928
# 2     2 Frobenius     0.426    0.146    0.938
# 3     3 Frobenius     0.428    0.130    0.912
# 4     5 Frobenius     0.436    0.173    0.927
# 5     6 Frobenius     0.522    0.188    1.14
# 6     4 Frobenius     0.522    0.175    1.03
# 7     1 Infinity      0.528    0.222    1.62
# 8     2 Infinity      0.540    0.185    1.32
# 9     3 Infinity      0.541    0.162    1.56
# 10    5 Infinity      0.552    0.165    1.44
# 11    4 Infinity      0.670    0.205    1.64
# 12    6 Infinity      0.676    0.241    2.03
# For 10000 samples, 500 permutations and 2000 MC samples
# A tibble: 12 × 5
#     row norm_type mean_norm min_norm max_norm
#    <int> <chr>         <dbl>    <dbl>    <dbl>
# 1     1 Frobenius     0.421    0.182    0.931
# 2     2 Frobenius     0.427    0.146    0.939
# 3     3 Frobenius     0.431    0.163    0.819
# 4     5 Frobenius     0.434    0.162    0.802
# 5     4 Frobenius     0.513    0.151    1.19
# 6     6 Frobenius     0.519    0.236    1.06
# 7     1 Infinity      0.531    0.221    1.63
# 8     2 Infinity      0.541    0.186    1.32
# 9     3 Infinity      0.547    0.172    1.35
# 10    5 Infinity      0.548    0.207    1.36
# 11    4 Infinity      0.650    0.192    1.41
# 12    6 Infinity      0.660    0.239    1.84

# Conclusions:
# ????

#### Compare CE fits for minimally different permutations ####

# Move 0 -> 100 observations from one block to the other, see
# how CE parameter estimates change (if at all)
# The uncertainty in CE estimates has to be the only source of uncertainty
# left here, since we've kept literally everything else equal!

dep_orig <- dep # store original dependence objects

# number of values to swap between blocks
# i_vec <- 1:100
# i_vec <- 0:100
i_vec <- 0:n_per_loc

i <- 1

# res <- lapply(i_vec, \(i) {
res <- mclapply(i_vec, \(i) {
  # print(paste("i =", i))
  system(sprintf("echo %s", paste(i, "done")))
  swp_loc1 <- sample(n_per_loc, i, replace = FALSE) # swap for block 1 to 2
  swp_loc2 <- sample((n_per_loc + 1):(n_per_loc * 2), i, replace = FALSE)

  # use the same data in both blocks, but swp swap_loc obs
  data_block_test_swp <- data_loc_test %>%
    group_by(name) |>
    mutate(row = row_number()) %>%
    mutate(
      block = case_when(
        row %in% swp_loc1 ~ "2",
        row %in% swp_loc2 ~ "1",
        TRUE ~ as.character(block)
      )
    ) |>
    ungroup()

  # # check that's done correctly
  # data_block_test_swp |>
  #   filter(name == "loc_1") |>
  #   filter(row %in% c(swp_loc1, swp_loc2)) |>
  #   pull(block) # should be 2s then 1s, rather than 1s then 2s

  data_block_test_swp <- select(data_block_test_swp, -row)

  data_laplace_block_test_swp <- trans_fun(
    data_block_test_swp,
    n_vars = 2, laplace_trans = TRUE
  )

  dep_swp <- lapply(
    data_laplace_block_test_swp,
    cecl_dep,
    cond_prob   = NULL,
    cond_val    = dep_val,
    fit_no_keef = TRUE
  )

  coef_diff <- coef(dep_swp[[1]])[, c("a", "b", "m", "s", "ll")] -
    coef(dep_swp[[2]])[, c("a", "b", "m", "s", "ll")]
  frob_norm_coef <- norm(as.matrix(coef_diff), type = "F")
  inf_norm_coef <- norm(as.matrix(coef_diff), type = "I")

  # Also calculate distances and norms for these!
  dist_swp <- lapply(seq_along(dep_swp), \(j) {
    cecl_dist(
      dep_obj  = dep_swp[[j]],
      marg_obj = data_laplace_block_test_swp[[j]],
      seed     = seed
    )
  })

  dist_diff_swp <- as.matrix(dist_swp[[2]]$dist_mat - dist_swp[[1]]$dist_mat)
  frob_norm_dist <- norm(dist_diff_swp, type = "F")
  inf_norm_dist <- norm(dist_diff_swp, type = "I")

  return(list(
    coef_diff = coef_diff,
    frob_coef = frob_norm_coef,
    inf_coef  = inf_norm_coef,
    frob_dist = frob_norm_dist,
    inf_dist  = inf_norm_dist
  ))
  # })
}, mc.cores = mc_cores)

# extract norm results
norm_res_df <- bind_rows(lapply(res, \(x) {
  data.frame(
    frob_coef = x$frob_coef,
    inf_coef  = x$inf_coef,
    frob_dist = x$frob_dist,
    inf_dist  = x$inf_dist
  )
}), .id = "n_swapped")

# plot results
norm_res_df

norm_res_df |>
  pivot_longer(
    cols = c(frob_coef, inf_coef, frob_dist, inf_dist),
    names_to = "norm_type",
    values_to = "norm_value"
  ) |>
  # slice(-1) |>
  mutate(
    across(c(n_swapped, norm_value), as.numeric),
    norm_type = case_when(
      norm_type == "frob_coef" ~ "Frobenius (coef)",
      norm_type == "inf_coef" ~ "Infinity (coef)",
      norm_type == "frob_dist" ~ "Frobenius (dist)",
      norm_type == "inf_dist" ~ "Infinity (dist)"
    )
  ) |>
  ggplot(aes(x = n_swapped, y = norm_value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~norm_type, scales = "free_y") +
  cecl_theme()


#### Compute norms for control case (data re-simulated with same corr) ####

# vector of block positions to explore (i.e., size of first block)
# block_vec <- seq(min_size, max_size, by = 100)
# block_vec <- c(1000, 1001)
# block_vec <- c(1000, 1000) # control setting, block size stays the same
block_vec <- c(n_per_loc, n_per_loc)

# run exploration across all block positions
set.seed(seed)
# debugonce(single_run_explore)
# debugonce(calc_dist)
# run once for various setups
test_norm_diff(data_loc, block_vec, n_mc = 500, laplace_trans = TRUE, cond_val = NULL)
# [1] 0.4363335 0.4380059 # F norm values
# [1] 0.4808696 0.4778487 # I norm values
# F           I
# 0.001672405 0.003020898
# norm values fairly similar (~0.001 diff), but quite large (~0.4) under no real difference!!
# increase MC samples
test_norm_diff(data_loc, block_vec, n_mc = 10000, laplace_trans = TRUE, cond_val = NULL)
# [1] 0.4363614 0.4357258
# [1] 0.4803298 0.4796136
# F            I
# 0.0006355648 0.0007161336
# difference in norms is significantly lower, but values are similar (~0.4)
# Using the same dependence quantile value across all locations, vars and blocks
test_norm_diff(data_loc, block_vec, n_mc = 10000, laplace_trans = TRUE, cond_val = dep_val)
# [1] 0.5958024 0.5971617
# [1] 0.6880635 0.6983510
# F           I
# 0.001359326 0.010287562
# norm values get much larger! From 0.4 to ~0.6, why is that?
# TODO Investigate!
# use same seed for sampling from Laplace distribution
test_norm_diff(data_loc, block_vec, n_mc = 10000, laplace_trans = TRUE, cond_val = dep_val, seed = seed)
# [1] 0.5968418 0.5968418
# [1] 0.6970967 0.6970967
# F I
# 0 0
# TODO Why does using the same seed but different dep thresh result in lowest
test_norm_diff(data_loc, block_vec, n_mc = 500, laplace_trans = TRUE, seed = seed) # same seed removes variability
# [1] 0.136286 0.136286
# [1] 0.1339932 0.1339932
# F I
# 0 0
# same with using same laplace samples for all (not just initial seed)
test_norm_diff(data_loc, block_vec, cond_val = dep_val, laplace_sample = laplace_sample)
# [1] 0.6590588 0.6590588
# [1] 0.7143129 0.7143129
# F I
# 0 0

# conclusion: Using same seed/laplace samples completely removes variability in this case
# TODO Look back into this conclusion!

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

test_norm_diff(block_vec1, seed = seed) # some variability, as desired!
test_norm_diff(block_vec1, seed = seed, cond_val = dep_val)

test1 <- replicate(100, test_norm_diff(block_vec1, seed = seed))
test2 <- replicate(100, test_norm_diff(block_vec1, seed = seed, cond_val = dep_val))

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
perm_test_res <- perm_test_fun(grid_vals)

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

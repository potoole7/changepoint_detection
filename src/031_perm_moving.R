#### Permutation test to identify a *single* changepoint in a time series ####

# For no change, get p-values of nearly 1, which is good! Works for
# testing no change case (i.e. control example)
# for diff = 0.4, all p-values are 0, so change is very significant
# for diff = 0.2, p-values for Frobenius are all 0, for infinity norm they are
# significant but not all 0, and the smallest p-value is at the true changepoint
# for diff = 0.1, p-values are not significant, and worst for true changepoint!
# TODO Try with smaller difference in correlation

# Finding the changepoint with our "exploration" shows numerous peaks,
# but the true changepoint is often among them, but not the largest!
# This may be due to the effect of blocking being larger than the effect of
# the changepoint itself.
# TODO How to improve this? Or estimate the effect of blocking!

# Speed up permutation test by separating testing and data
# splitting (done)
# TODO Test with different correlation changes
# TODO Investigate whether p-values are calculated correctly
# TODO Investigate whether I need to do Laplace transformation after
# Test permutation test for fine grid (done)
# TODO How to decide on CE threshold?

# Currently, block size appears to have a larger effect than the
# changepoint!
# Try increasing difference

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
# diff <- 0.8
# diff <- 0 # try for no difference, see what happens!
diff <- 0.2 # small difference
# diff <- 0.1 # smaller difference
mc_cores <- detectCores() - 1 # number of cores for parallel processing

#### Functions ####

# Function to generate multivariate t data with specified correlation
gen_t <- \(cor_t, n) {
  # generate data
  cop_t <- tCopula(param = cor_t, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}

# Function to convert data to cecl_marg object
trans_fun <- \(data, n_vars) {
  # Calculate CE and cluster for all time blocks:

  # We know the exact GPD values, so use these to transform to Laplace
  locs <- unique(data$name)
  vars <- names(data)[1:n_vars]
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
calc_dist <- \(marg) {
  # fit CE models
  dep <- lapply(
    marg,
    cecl_dep,
    cond_prob = 0.9,
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
      n_mc = 500
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
single_run_explore <- \(data_loc, i) {
  data_block <- data_loc %>%
    group_by(name) |>
    mutate(block = ifelse(row_number() <= i, "1", "2")) |>
    ungroup()

  # check blocking is correct
  # table(data_block$block, data_block$name)

  # convert to Laplace and store marginals
  # TODO Should we do this before blocking???
  data_laplace_block <- trans_fun(data_block, n_vars)

  # check you've transformed correctly
  # par(mfrow = c(1, 2))
  # plot(data_laplace_block[[1]][[1]][[1]])
  # plot(data_laplace_block[[2]][[1]][[1]])
  # par(mfrow = c(1, 1))

  dist <- suppressMessages(calc_dist(data_laplace_block))

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

# # Function to perform a single permutation test
# Function to generate data for a single permutation test
perm_test_prep <- \(
  data_laplace_loc # ,
  # type = c("norm", "clustering"),
  # norm_type = "F"
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
    cond_prob = 0.9,
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
      n_mc = 500
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
  ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n))
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


#### Exploratory ####

# want to move block delimiter by 1 and compute norms at each, as exploration
# of where changepoint may be

min_size <- 200 # minimum size of each block
max_size <- data_loc |>
  filter(name == "loc_1") |>
  nrow() - min_size

# vector of block positions to explore (i.e., size of first block)
block_vec <- seq(min_size, max_size, by = 100)

# run exploration across all block positions
set.seed(seed)
norm_vals <- mclapply(block_vec, \(i) {
  single_run_explore(data_loc, i)
}, mc.cores = mc_cores)

# extract frobenius and infinity norms
norm_vals_frob <- sapply(norm_vals, \(x) x$frob)
norm_vals_inf <- sapply(norm_vals, \(x) x$inf)

plot_norm(block_vec, norm_vals_frob, norm_vals_inf)
# Frobenius norm finds peak at 1000 (and at 700),
# infinity has a smaller peak for true changepoint than at 700

# find peaks in norm values
find_peaks(norm_vals_frob, block_vec)
find_peaks(norm_vals_inf, block_vec)

# get block positions of maximum norms (excluding first value)
(max_frob <- block_vec[which.max(norm_vals_frob[-1]) + 1])
(max_inf <- block_vec[which.max(norm_vals_inf[-1]) + 1])

# run again for blocks from 600 to 1200, with increments of 10
block_vec_full <- seq(
  min(c(max_frob, max_inf)) - 100,
  max(c(max_frob, max_inf)) + 100,
  by = 10
)

norm_vals_full <- mclapply(block_vec_full, \(i) {
  single_run_explore(data_loc, i)
}, mc.cores = mc_cores)

plot_norm(
  block_vec_full,
  sapply(norm_vals_full, \(x) x$frob),
  sapply(norm_vals_full, \(x) x$inf)
)
# doesn't seem to find the changepoint very well...
# Maybe transform to Laplace using empirical method in each loop?

find_peaks(
  sapply(norm_vals_full, \(x) x$frob),
  block_vec_full
)
find_peaks(
  sapply(norm_vals_full, \(x) x$inf),
  block_vec_full
)
# lots of peaks! But both include the treu chanegepoint at 1000 at least

# finally, try from 900 to 1100 by 1
block_vec_fine <- seq(900, 1100, by = 1)
norm_vals_fine <- mclapply(block_vec_fine, \(i) {
  single_run_explore(data_loc, i)
}, mc.cores = mc_cores)
plot_norm(
  block_vec_fine,
  sapply(norm_vals_fine, \(x) x$frob),
  sapply(norm_vals_fine, \(x) x$inf)
)

find_peaks(
  sapply(norm_vals_fine, \(x) x$frob),
  block_vec_fine
)
find_peaks(
  sapply(norm_vals_fine, \(x) x$inf),
  block_vec_fine
)
# lots of peaks determined! Nothing very conclusive though


#### Permutation test (test) ####

# for points in "full" grid, we will run permutation test 100 times & determine
# p-values at each
# grid_vals <- block_vec_full
# grid_vals <- c(990, 1000, 1010) # smaller grid to test
grid_vals <- seq(950, 1050, by = 10) # smaller grid to test
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
    dist <- perm_test_prep(data_laplace_block)
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

# saveRDS(perm_test_res, file = "perm_test_res_moving.rds")
saveRDS(
  perm_test_res,
  # file = "perm_test_res_moving.rds"
  file = paste0("perm_test_res_moving_diff_", diff, ".rds")
)
# load
perm_test_res <- readRDS(
  file = paste0("perm_test_res_moving_diff_", diff, ".rds")
)

# plotting, we can see that change is very significant!
hist(perm_test_res[[1]]$perm_norms_inf, xlim = c(0, perm_test_res[[1]]$norm_orig_inf))
abline(v = perm_test_res[[1]]$norm_orig_inf, col = "red", lwd = 2)

# all p-values are 0!
# TODO Choose smaller value of difference in correlation, too large here!
cbind(
  grid_vals,
  sapply(perm_test_res, \(x) x$p_value_inf)
) |>
  as.data.frame()

cbind(
  grid_vals,
  sapply(perm_test_res, \(x) x$p_value_frob)
) |>
  as.data.frame()

# p-values aren't telling us much, try calculate calculate difference from
# mean/median permutation value

test_fun(perm_test_res, norm_type = "inf", stat = "mean") |>
  cbind(grid_vals) |>
  as.data.frame() |>
  setNames(c("diff_mean", "block_pos")) |>
  arrange(desc(diff_mean)) # smallest difference at true changepoint!

test_fun(perm_test_res, norm_type = "frob", stat = "mean") |>
  cbind(grid_vals) |>
  as.data.frame() |>
  setNames(c("diff_mean", "block_pos")) |>
  arrange(desc(diff_mean)) # 4th smallest difference at true changepoint!

# try with median
test_fun(perm_test_res, norm_type = "inf", stat = "median") |>
  cbind(grid_vals) |>
  as.data.frame() |>
  setNames(c("diff_median", "block_pos")) |>
  arrange(desc(diff_median)) # also worst!

test_fun(perm_test_res, norm_type = "frob", stat = "median") |>
  cbind(grid_vals) |>
  as.data.frame() |>
  setNames(c("diff_median", "block_pos")) |>
  arrange(desc(diff_median)) # also 8th worst

# Method is poor for detecting changepoint here, likely due to effect of
# blocking being larger than effect of changepoint
# Or perhaps effect of changepoint is far too high (for diff = 0.4), and so
# permutation test always finds significant difference

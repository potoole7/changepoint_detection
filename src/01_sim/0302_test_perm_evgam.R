#### Test permutation test for spatial data ####

# TODO also just pull alpha and beta values, not just norms??

# Also look at differences for just alpha and beta (running now)
# Perform permutation test
# - Allow fitting CE with evgam in `cacl_dist` (done)
# - Run on server (done)
# - Analyse results (done)

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(parallel)
library(evgam)
library(ggridges)
source("src/00_functions.R")

#### metadata ####

seed <- 123
n_vars <- 2 # number of variables
# n_locs <- 12 # number of locations
n_locs <- 60
# n_locs <- 60
n_clust <- 3 # number of clusters
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per loc times n_locs per group
n_per_block <- n * n_clust # number of obs per block (two blocks here)
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group

diff <- 0 # try for no difference in correlation (but stochastically different data!
# diff <- 0.4
# use_test_case <- TRUE # try for absolutely no difference (i.e. same exact data in both blocks)
use_test_case <- FALSE
mc_cores <- detectCores() - 1 # number of cores for parallel processing
# mc_cores <- 14


#### Generate data ####

increments <- c(0, 1)
data_loc <- bind_rows(lapply(seq_along(increments), \(i) {
  cor_spec <- cor_t_vec
  cor_spec[[1]] <- cor_spec[[1]] + (diff * increments[[i]])

  # generate data with spatially varying correlation
  set.seed(seed)
  clust_centres <- data.frame(
    # cor = cor_t_vec,
    cor = cor_spec,
    x = runif(n_clust, 0, 5),
    y = runif(n_clust, 0, 5)
  )

  # how many locations per cluster
  locs_per_clust <- n_locs / n_clust
  stopifnot(locs_per_clust == floor(locs_per_clust))

  # generate locations around cluster centres, if not provided
  locs <- clust_centres %>%
    mutate(cluster = row_number()) %>%
    tidyr::uncount(locs_per_clust) %>% # replicate centres locs_per_clust times
    group_by(cor, x, y) %>%
    # generate locations around centre with some noise
    mutate(
      loc_id = row_number(),
      x_loc  = x + rnorm(n(), 0, 0.8),
      y_loc  = y + rnorm(n(), 0, 0.8)
    ) %>%
    ungroup() %>%
    mutate(loc_global_id = row_number()) %>%
    select(cluster, loc_global_id, x_loc, y_loc)


  print(cor_spec)
  # generate data
  # ret <- do.call(rbind, lapply(cor_spec, gen_t, n = n, laplace_trans = TRUE))
  ret <- gen_spatial_data(
    cor_t_vec = cor_spec,
    n_clust = n_clust,
    n_locs = n_locs,
    n_per_loc = n_per_loc,
    sd_loc = 0.8,
    clust_centres = clust_centres,
    locs = locs,
    laplace_trans = TRUE
  )
  ret |>
    mutate(
      name = rep(paste0("loc_", 1:n_locs), each = nrow(ret) / n_locs),
      block = i # time block
    )
}))


locs <- unique(data_loc$name)

# check data is correct
table(data_loc$block)
table(data_loc$name)

# should only be of length n_locs (and only have one per location)
sort(unique(paste(data_loc$name, data_loc$x_loc, sep = " --- "))) |>
  head()

# plot
data_loc |>
  ggplot(aes(x = x_loc, y = y_loc, colour = rho)) +
  geom_point(size = 4) +
  facet_wrap(~block) +
  cecl_theme(nejm_pal = FALSE)
# looks good!

# save original data and remove block info
data_loc_orig <- data_loc
data_loc <- select(data_loc, -block)


#### Precalculations ####

# TODO Move to it's own code block
# use the same dependence threshold across all variables and locations
# (dep_val <- qgpd(0.8, u = 0, xi = -0.05, sigma = 1)) #
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


#### Compare CE fits for minimally different permutations ####

# use the exact same data in both blocks, as a test
if (use_test_case && diff == 0) {
  data_loc_test <- data_loc_orig |>
    filter(block == 1) |>
    bind_rows(
      data_loc_orig |>
        filter(block == 1) |>
        mutate(block = 2)
    )
} else {
  data_loc_test <- data_loc_orig
}

# Move 0 -> 100 observations from one block to the other, see
# how CE parameter estimates change (if at all)
# The uncertainty in CE estimates has to be the only source of uncertainty
# left here, since we've kept literally everything else equal!

# number of values to swap between blocks
# i_vec <- 0:n_per_loc
i_vec <- 0:100

# testing function
# i <- 1
# data <- data_loc_test

# Function to calculate norms for coefficient differences
norm_fun <- \(coef, which = c("a", "b", "m", "s")) {
  # coef_diff <- coef[[1]][, c("a", "b", "m", "s")] -
  #   coef[[2]][, c("a", "b", "m", "s")]
  coef_diff <- coef[[1]][, which] - coef[[2]][, which]
  frob_norm_coef <- norm(as.matrix(coef_diff), type = "F")
  inf_norm_coef <- norm(as.matrix(coef_diff), type = "M")

  return(list(
    coef_diff = coef_diff,
    frob = frob_norm_coef,
    inf = inf_norm_coef
  ))
}


# TODO Speed up without breaking code
incremental_perm <- \(data, i_vec, fit_evgam = FALSE, nruns = 1, use_start = FALSE, aLow = -1, mc_cores = getOption("mc.cores", 2L)) {
  # Parallel setup
  # Can't use if having start values!
  # apply_fun <- if (use_start) lapply else parallel::mclapply

  # initial start values (overridden in loop)
  # TODO (or are they because of mapply?)
  start0 <- list(
    c("a" = 0.01, "b" = 0.01),
    c("a" = 0.01, "b" = 0.01)
  )

  # function within apply
  res_fun <- \(i) {
    # browser()
    start <- start0
    swp_loc1 <- sample(n_per_loc, i, replace = FALSE) # swap for block 1 to 2
    swp_loc2 <- sample((n_per_loc + 1):(n_per_loc * 2), i, replace = FALSE)

    # use the same data in both blocks, but swp swap_loc obs
    data_block_test_swp <- data %>%
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

    data_block_test_swp_fit <- select(
      data_block_test_swp,
      -any_of(c("row", "cluster", "x_loc", "y_loc", "rho"))
    )

    data_laplace_block_test_swp <- trans_fun(
      data_block_test_swp_fit,
      n_vars = 2, laplace_trans = TRUE
    )

    # fit with "vanilla" CE model
    # print(start)
    if (fit_evgam == FALSE) {
      dep_swp <- lapply(seq_along(data_laplace_block_test_swp), \(k) {
        cecl_dep(
          # data_block_test_swp[[k]],
          data_laplace_block_test_swp[[k]],
          cond_prob = NULL,
          cond_val = dep_val,
          fit_no_keef = TRUE,
          nruns = nruns,
          aLow = aLow,
          start = start[[k]]
        )
      })

      # update start values with those for previous iteration
      # if (use_start == TRUE) {
      #   # start0 <<- coef(dep_swp)
      #   start0 <<- lapply(dep_swp, coef)
      # }
    } else if (fit_evgam == TRUE) {
      # first, put marginal data into form for evgam fitting
      marg_join_lst <- lapply(seq_along(data_laplace_block_test_swp), \(k) {
        x <- data_laplace_block_test_swp[[k]]
        bind_rows(lapply(seq_along(x$transformed), \(j) {
          data.frame(x$transformed[[j]]) |>
            mutate(name = paste0("loc_", j))
        })) |>
          left_join(
            data |>
              filter(block == k) |>
              distinct(name, x_loc, y_loc),
            by = "name"
          ) |>
          rename(x = x_loc, y = y_loc) |>
          select(-name)
      })

      dep_swp <- lapply(marg_join_lst, \(x) {
        x1_cond <- fit_evgam(
          df = x,
          dep_val = qlaplace(p = 0.8),
          var = "X2",
          cond_var = "X1"
        )

        x2_cond <- fit_evgam(
          df = x,
          dep_val = qlaplace(p = 0.8),
          var = "X1",
          cond_var = "X2"
        )

        # join and change to `cecl_dep` format
        coef_evgam <- bind_rows(x1_cond$predictions, x2_cond$predictions) |>
          arrange(name, var, cond_var)
        rownames(coef_evgam) <- NULL
        as_cecl_dep(coef_evgam)
      })
    }

    # can't compare ll as it's not estimated in evgam method
    coef <- lapply(dep_swp, coef)
    # coef_diff <- coef[[1]][, c("a", "b", "m", "s")] -
    #   coef[[2]][, c("a", "b", "m", "s")]
    # frob_norm_coef <- norm(as.matrix(coef_diff), type = "F")
    # inf_norm_coef <- norm(as.matrix(coef_diff), type = "I")
    diff <- norm_fun(coef, which = c("a", "b", "m", "s"))
    coef_diff <- diff$coef_diff
    frob_norm_coef <- diff$frob
    inf_norm_coef <- diff$inf

    # Now look at just alpha and beta (i.e. ignore changes in nuisance params)
    diff_a <- norm_fun(coef, which = c("a"))
    coef_diff_a <- diff$coef_diff
    frob_norm_coef_a <- diff$frob
    inf_norm_coef_a <- diff$inf

    diff_b <- norm_fun(coef, which = c("b"))
    coef_diff_b <- diff$coef_diff
    frob_norm_coef_b <- diff$frob
    inf_norm_coef_b <- diff$inf

    diff_m <- norm_fun(coef, which = c("m"))
    coef_diff_m <- diff$coef_diff
    frob_norm_coef_m <- diff$frob
    inf_norm_coef_m <- diff$inf

    diff_s <- norm_fun(coef, which = c("s"))
    coef_diff_s <- diff$coef_diff
    frob_norm_coef_s <- diff$frob
    inf_norm_coef_s <- diff$inf

    diff_ab <- norm_fun(coef, which = c("a", "b"))
    coef_diff_ab <- diff$coef_diff
    frob_norm_coef_ab <- diff$frob
    inf_norm_coef_ab <- diff$inf

    # Also calculate distances and norms for these!
    dist_swp <- lapply(seq_along(dep_swp), \(j) {
      cecl_dist(
        dep_obj = dep_swp[[j]],
        marg_obj = data_laplace_block_test_swp[[j]],
        # seed     = seed
        laplace_sample = laplace_sample
      )
    })

    dist_diff_swp <- as.matrix(dist_swp[[2]]$dist_mat - dist_swp[[1]]$dist_mat)
    frob_norm_dist <- norm(dist_diff_swp, type = "F")
    inf_norm_dist <- norm(dist_diff_swp, type = "M")

    system(sprintf("echo %s", paste(i, "done")))

    return(list(
      coef = coef,
      coef_diff = coef_diff,
      frob_coef = frob_norm_coef,
      inf_coef = inf_norm_coef,
      frob_norm_coef_a = frob_norm_coef_a,
      inf_norm_coef_a = inf_norm_coef_a,
      frob_norm_coef_b = frob_norm_coef_b,
      inf_norm_coef_b = inf_norm_coef_b,
      frob_norm_coef_m = frob_norm_coef_m,
      inf_norm_coef_m = inf_norm_coef_m,
      frob_norm_coef_s = frob_norm_coef_s,
      inf_norm_coef_s = inf_norm_coef_s,
      frob_norm_coef_ab = frob_norm_coef_ab,
      inf_norm_coef_ab = inf_norm_coef_ab,
      frob_dist = frob_norm_dist,
      inf_dist = inf_norm_dist
    ))
    return(res)
  }

  # run once and use start values from initial fit for subsequent runs
  if (use_start) {
    res_first <- res_fun(i_vec[[1]])
    start0 <- res_first$coef
    i_vec <- i_vec[-1]
  }

  ret <- mclapply(i_vec, res_fun)
  # ret <- lapply(i_vec, res_fun)
  if (use_start) {
    ret <- c(list(res_first), ret)
  }
  return(ret)
}

# run for vanilla and evgam CE fits
res_ce <- incremental_perm(data_loc_test, i_vec, fit_evgam = FALSE) # standard
# In optimisation fit model a second time
res_ce2 <- incremental_perm(data_loc_test, i_vec, fit_evgam = FALSE, nruns = 2)
# Also add that alpha minimum value must be 0, since we know corr is positive
res_ce2_alow <- incremental_perm(data_loc_test, i_vec, fit_evgam = FALSE, nruns = 2, aLow = 0)
# use start values form previous iteration (slow, difficult to parallelise)
res_ce2_alow_st <- incremental_perm(data_loc_test, i_vec, fit_evgam = FALSE, nruns = 2, use_start = TRUE, aLow = 0)
# res_evgam <- incremental_perm(data_loc_test, i_vec, fit_evgam = TRUE)
# res_lst <- list("cecl" = res_ce, "evgam" = res_evgam)
res_lst <- list(
  "ce"          = res_ce,
  "ce2"         = res_ce2,
  "ce2_alow"    = res_ce2_alow,
  "ce2_alow_st" = res_ce2_alow_st,
  "evgam"       = res_evgam
)

# save results; quite slow to run for evgam case!
saveRDS(res_lst, paste0("0302_incremental_perm_results_diff_", diff, ".rds"))

# diff <- 0
# diff <- 0.4
res_lst <- readRDS(paste0("0302_incremental_perm_results_diff_", diff, ".rds"))
# res_ce <- res_lst$cecl
# res_evgam <- res_lst$evgam

# extract norm results
# norm_res_df <- bind_rows(lapply(res_ce, \(x) {
#   data.frame(
#     frob_coef = x$frob_coef,
#     inf_coef  = x$inf_coef,
#     frob_dist = x$frob_dist,
#     inf_dist  = x$inf_dist
#   )
# }), .id = "n_swapped")

# lst_names <- c("cecl", "evgam")
lst_names <- c("ce", "ce2", "ce2_alow", "ce2_alow_st", "evgam")
# lst_names <- c("ce", "ce2", "ce2_alow", "ce2_alow_st")
# lst_names <- c("cecl")
# norm_res_df <- bind_rows(lapply(seq_along(res_lst), \(i) {
norm_res_df <- bind_rows(mclapply(seq_along(res_lst), \(i) {
  x <- res_lst[[i]]
  bind_rows(lapply(x, \(y) {
    data.frame(
      frob_coef    = y$frob_coef,
      inf_coef     = y$inf_coef,
      frob_coef_a  = y$frob_norm_coef_a,
      inf_coef_a   = y$inf_norm_coef_a,
      frob_coef_b  = y$frob_norm_coef_b,
      inf_coef_b   = y$inf_norm_coef_b,
      frob_coef_m  = y$frob_norm_coef_m,
      inf_coef_m   = y$inf_norm_coef_m,
      frob_coef_s  = y$frob_norm_coef_s,
      inf_coef_s   = y$inf_norm_coef_s,
      frob_coef_ab = y$frob_norm_coef_ab,
      inf_coef_ab  = y$inf_norm_coef_ab,
      frob_dist    = y$frob_dist,
      inf_dist     = y$inf_dist,
      method       = lst_names[[i]]
    )
  }), .id = "n_swapped")
}))

# plot results
norm_res_df

norm_res_df_long <- norm_res_df |>
  pivot_longer(
    # cols = c(frob_coef, inf_coef, frob_dist, inf_dist),
    # cols = c(frob_coef, inf_coef, frob_coef_ab, inf_coef_ab, frob_dist, inf_dist),
    # cols = c(frob_coef, inf_coef, frob_coef_a, inf_coef_a, frob_coef_b, inf_coef_b, frob_coef_ab, inf_coef_ab, frob_dist, inf_dist),
    # cols = c(frob_coef, inf_coef, frob_coef_a, inf_coef_a, frob_coef_b, inf_coef_b, frob_coef_m, inf_coef_m, frob_coef_s, inf_coef_s, frob_coef_ab, inf_coef_ab, frob_dist, inf_dist),
    cols = c(frob_coef, inf_coef, frob_coef_a, inf_coef_a, frob_coef_b, inf_coef_b, frob_coef_m, inf_coef_m, frob_coef_s, inf_coef_s, frob_dist, inf_dist),
    names_to = "norm_type",
    values_to = "norm_value"
  ) |>
  # slice(-1) |>
  mutate(
    across(c(n_swapped, norm_value), as.numeric),
    method = case_when(
      method == "ce" ~ "CE",
      method == "ce2" ~ "CE (nruns=2)",
      method == "ce2_alow" ~ "CE (nruns=2, aLow=0)",
      method == "ce2_alow_st" ~ "CE (nruns=2, aLow=0, start)",
      method == "evgam" ~ "evgam CE"
    ),
    norm_type = case_when(
      norm_type == "frob_coef" ~ "Frobenius (coef)",
      norm_type == "inf_coef" ~ "Infinity (coef)",
      norm_type == "frob_coef_a" ~ "Frobenius (coef) (alpha)",
      norm_type == "inf_coef_a" ~ "Infinity (coef) (alpha)",
      norm_type == "frob_coef_b" ~ "Frobenius (coef) (beta)",
      norm_type == "inf_coef_b" ~ "Infinity (coef) (beta)",
      norm_type == "frob_coef_m" ~ "Frobenius (coef) (mu)",
      norm_type == "inf_coef_m" ~ "Infinity (coef) (mu)",
      norm_type == "frob_coef_s" ~ "Frobenius (coef) (sigma)",
      norm_type == "inf_coef_s" ~ "Infinity (coef) (sigma)",
      norm_type == "frob_coef_ab" ~ "Frobenius (coef) (alpha,beta)",
      norm_type == "inf_coef_ab" ~ "Infinity (coef) (alpha,beta)",
      norm_type == "frob_dist" ~ "Frobenius (dist)",
      norm_type == "inf_dist" ~ "Infinity (dist)"
    )
  ) |>
  # ensure method is a factor with correct ordering for plotting
  mutate(
    method = factor(method, levels = unique(method))
  )

# fitted lines through dots
norm_res_df_long |>
  ggplot(aes(x = n_swapped, y = norm_value, colour = method)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  facet_wrap(~norm_type, scales = "free") +
  labs(
    x = "# observations swapped between blocks",
    y = "Norm value"
  ) +
  cecl_theme()

# also plot histograms
norm_res_df_long |>
  ggplot(aes(x = norm_value, fill = method)) +
  geom_density(alpha = 0.5) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.5, bins = 30) +
  facet_wrap(~norm_type, scales = "free") +
  labs(x = "Norm value", y = "Density") +
  cecl_theme()

# also plot ridges
norm_res_df_long |>
  ggplot(aes(x = norm_value, y = method, fill = method)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~norm_type, scales = "free") +
  labs(x = "Norm value", y = "Method") +
  cecl_theme() +
  # remove long method names from y-axis
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  )


#### Run permutation test with vanilla and evgam ####

# Note: differences for the same seed are both small
# TODO Improve summary of these results, what is `diff` here?

# half in each block, run twice to see stochasitc differences
block_vec <- c(n_per_loc, n_per_loc)

set.seed(seed)
# run once (vanilla CE)
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_evgam = FALSE
)
# [1] 0.03208471 0.12292635
# [1] 0.06448799 0.15908541
# $diff
# F          I
# 0.09084164 0.09459742
# CE with starting values etc
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_start = TRUE, nruns = 2, aLow = 0,
  use_evgam = FALSE
)
# [1] 0.30150128 0.04040002
# [1] 0.36832971 0.04641947
# $diff
# F         I
# 0.2611013 0.3219102
# evgam
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_evgam = TRUE, data_loc = data_loc_test
)
# [1] 0.1559021 0.1352027
# [1] 0.1566529 0.1339805
# $diff
# F          I
# 0.02069936 0.02267248

# result:
# Norm values are definitely smaller when using evgam to fit the CE model
# Differences between the two for the same block size though are roughly
# the same
# Somehow, "improving" CE with nruns etc makes it worse? Could just be seed ...

#### Compute norms for times 1000 and 1001 ####

# This time, we hope to see a small difference in norms, as there is a change in
# correlation structure between these two time blocks
block_vec1 <- c(1000, 1001)
set.seed(seed)
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec1,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_evgam = FALSE
)
# [1] 0.03208471 0.15811811
# [1] 0.06448799 0.31915638
# $diff
# F         I
# 0.1260334 0.2546684
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_start = TRUE, nruns = 2, aLow = 0,
  use_evgam = FALSE
)
# [1] 0.30150128 0.04040002
# [1] 0.36832971 0.04641947
# $diff
# F         I
# 0.2611013 0.3219102
test_norm_diff(
  data = select(data_loc_test, X1, X2, name, block),
  block_vec1,
  n_mc = 500, laplace_trans = TRUE,
  cond_val = dep_val, # seed = seed,
  use_evgam = TRUE, data_loc = data_loc_test
)
# [1] 0.1559021 0.2156945
# [1] 0.1566529 0.3671678
# $diff
# F          I
# 0.05979246 0.21051486
# results: Again, norms are smaller, and differences are slightly smaller

# TODO Run many times!
# test1 <- replicate(100, test_norm_diff(block_vec1, seed = seed))
# test2 <- replicate(100, test_norm_diff(block_vec1, seed = seed, cond_val = dep_val))
#
# # check results
# list(test1 = test1, test2 = test2) |>
#   sapply(\(x) apply(x, 1, mean))

#### Test permutation test ####

grid_vals <- seq(950, 1050, by = 10)
# grid_vals <- c(1000, 1001)
# grid_vals <- c(1000, 1000)
# grid_vals <- c(999, 1000, 1001)
# TODO Increase!
n_perm <- 100 # number of permutations

set.seed(seed)
# perm_test_res <- perm_test_fun(grid_vals)
# debugonce(perm_test_fun)
perm_test_res_ce <- perm_test_fun(
  data = select(data_loc_test, X1, X2, name, block),
  grid_vals = grid_vals,
  n_perm = n_perm,
  laplace_trans = TRUE,
  cond_val = dep_val, laplace_sample = laplace_sample # seed = seed
)

# Optimise CeCl with aLow, nruns, use_start
perm_test_res_ce2 <- perm_test_fun(
  data = select(data_loc_test, X1, X2, name, block),
  grid_vals = grid_vals,
  n_perm = n_perm,
  laplace_trans = TRUE,
  aLow = 0, nruns = 2, use_start = TRUE,
  cond_val = dep_val, laplace_sample = laplace_sample # seed = seed
)

# perform using evgam
perm_test_res_evgam <- perm_test_fun(
  data = select(data_loc_test, X1, X2, name, block),
  grid_vals = grid_vals,
  n_perm = n_perm,
  laplace_trans = TRUE,
  cond_val = dep_val, laplace_sample = laplace_sample, # seed = seed,
  use_evgam = TRUE, data_loc = data_loc_test
)

perm_test_res <- list(
  cecl  = perm_test_res_ce,
  cecl2 = perm_test_res_ce2,
  evgam = perm_test_res_evgam
)
# save results
saveRDS(
  perm_test_res,
  paste0("0302_perm_test_results_diff_", diff, ".rds")
)

# diff <- 0
# diff <- 0.4
perm_test_res <- readRDS(paste0("0302_perm_test_results_diff_", diff, ".rds"))
perm_test_res_ce <- perm_test_res[[1]]

# plot
hist(
  # perm_test_res_ce[[1]]$perm_norms_frob,
  perm_test_res_ce[[1]]$perm_norms_inf,
  xlim = c(
    0,
    # max(c(perm_test_res_ce[[1]]$perm_norms_frob, perm_test_res_ce[[1]]$norm_orig_frob))
    max(c(perm_test_res_ce[[1]]$perm_norms_inf, perm_test_res_ce[[1]]$norm_orig_inf))
  )
)
# abline(v = perm_test_res_ce[[1]]$norm_orig_frob, col = "red", lwd = 2)
abline(v = perm_test_res_ce[[1]]$norm_orig_inf, col = "red", lwd = 2)
perm_test_res_ce[[1]]$p_value_inf

sapply(perm_test_res_ce, `[[`, "p_value_frob")
# If difference = 0, then p-values all 1; if diff = 0.4, then p-values all 0!


# extract results
# i <- j <- 1
perm_test_res_df <- bind_rows(lapply(seq_along(perm_test_res), \(i) {
  x <- perm_test_res[[i]]
  bind_rows(lapply(seq_along(x), \(j) {
    data.frame(
      norm_value = c(x[[j]]$perm_norms_inf, x[[j]]$perm_norms_frob),
      norm_orig = c(
        rep(x[[j]]$norm_orig_inf, length(x[[j]]$perm_norms_inf)),
        rep(x[[j]]$norm_orig_frob, length(x[[j]]$perm_norms_frob))
      ),
      p_value = c(
        rep(x[[j]]$p_value_inf, length(x[[j]]$perm_norms_inf)),
        rep(x[[j]]$p_value_frob, length(x[[j]]$perm_norms_frob))
      ),
      norm_type = rep(c("Infinity", "Frobenius"), each = length(x[[j]]$perm_norms_inf)),
      method = names(perm_test_res)[i],
      grid_val = grid_vals[[j]]
    )
  }))
}))

# plot
perm_test_res_df |>
  ggplot(aes(x = norm_value, fill = method)) +
  geom_density(alpha = 0.5) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.5,
    bins = 30
  ) +
  facet_wrap(~norm_type, scales = "free") +
  labs(x = "Norm value", y = "Density") +
  cecl_theme()
# values are reliably lower for evgam method, as desired, which is good!

# plot with ridges for each grid value
perm_test_res_df |>
  ggplot(aes(x = norm_value, y = as.factor(grid_val), fill = method)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  facet_wrap(~norm_type, scales = "free") +
  labs(x = "Norm value", y = "Grid value") +
  cecl_theme()

# which is the highest peak?
peak_df <- perm_test_res_df %>%
  group_by(norm_type, grid_val, method) %>%
  summarise(
    peak_norm_value = {
      d <- density(norm_value, na.rm = TRUE)
      d$x[which.max(d$y)]
    },
    peak_density = {
      d <- density(norm_value, na.rm = TRUE)
      max(d$y)
    },
    .groups = "drop"
  )
peak_df |>
  arrange(norm_type, method, desc(peak_density))

# plot
peak_df |>
  ggplot(aes(x = grid_val, y = peak_norm_value, colour = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(method ~ norm_type, scales = "free") +
  geom_vline(
    xintercept = (nrow(data_loc_orig) / n_locs) / 2, # TODO Need to add this as an argument!
    lty = 2, lwd = 2, colour = "red"
  ) +
  labs(
    x = "Grid value",
    y = "Peak norm value"
  ) +
  cecl_theme()
# not at 1000 for any of them, unfortunately ..

perm_test_res_df |>
  # pivot_longer(c(norm_orig, norm_value, p_value)) |>
  ggplot(aes(x = grid_val, y = norm_orig, colour = method)) +
  # ggplot(aes(x = grid_val, y = value, colour = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~norm_type, scales = "free") +
  # facet_wrap(~ name + norm_type, scales = "free") +
  geom_vline(
    xintercept = (nrow(data_loc_orig) / n_locs) / 2, # TODO Need to add this as an argument!
    lty = 2, lwd = 2, colour = "red"
  ) +
  scale_y_continuous(limits = c(0, NA)) + # set y-axis to start at 0
  labs(
    x = "Grid value",
    y = "Original norm value"
  ) +
  cecl_theme()
# Still not quite right, norm value should peak at the point where the blocks are equal size
# as that's the "correct" changepoint

# for no difference, norm value is very low, as desired, since data is the same!
# TODO Investigate, shouldn't the seed be the same, meaning the distance should be 0?

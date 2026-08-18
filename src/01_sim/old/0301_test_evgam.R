#### ####

#### Testing permutation test method for determining changepoints ####

# TODO CE extremes fits vary quite a bit when permutating; perhaps
# fitting an alternative (such as with evgam) would be more stable?
# TODO Try this out!! Would need to generate spatial data

# TODO Generate spatial data (as in other script)
# TODO Fit "vanilla" CE with CeCl
# TODO Fit alternative CE with evgam (also use CeCl::as_cecl_dep to convert)
# TODO Test that method still clusters correctly (for "vanilla" and/or evgam)
# TODO Compare results (evgam fits should be more stable??)

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(evgam) # for alternative CE fitting
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

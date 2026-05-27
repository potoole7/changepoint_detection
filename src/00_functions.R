#### Functions ####

# Function to generate multivariate t data with specified correlation
gen_t <- \(cor_t, n, n_vars = 2, df_t = 3, laplace_trans = FALSE) {
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
calc_dist <- \(
  marg,
  n_mc = 500,
  cond_prob = 0.9,
  cond_val = NULL,
  use_dth = FALSE,
  cond_var = NULL,
  use_evgam = FALSE,
  data_loc = NULL,
  aLow = -1,
  nruns = 1,
  ret_dep = FALSE,
  start = list(c(a = 0.01, b = 0.01), c(a = 0.01, b = 0.01)),
  ...
) {
  if (!is.null(cond_val)) {
    cond_prob <- NULL
  }

  # fit CE models
  if (use_evgam == FALSE) {
    dep <- lapply(seq_along(marg), \(i) {
      cecl_dep(
        obj         = marg[[i]],
        cond_prob   = cond_prob,
        cond_val    = cond_val,
        cond_var    = cond_var,
        aLow        = aLow,
        nruns       = nruns,
        fit_no_keef = TRUE,
        start       = start[[i]]
      )
    })
    # skip bad permutations
    if (any(sapply(dep, \(x) any(is.na(unlist(x$dependence)))))) {
      return(NA)
    }
  } else if (use_evgam) {
    if (is.null(data_loc)) {
      stop("data_loc must be provided when use_evgam = TRUE")
    }

    # first, put marginal data into form for evgam fitting
    marg_join_lst <- lapply(seq_along(marg), \(k) {
      x <- marg[[k]]
      bind_rows(lapply(seq_along(x$transformed), \(j) {
        data.frame(x$transformed[[j]]) |>
          mutate(name = paste0("loc_", j))
      })) |>
        left_join(
          # data |>
          data_loc |>
            filter(block == k) |>
            distinct(name, x_loc, y_loc),
          by = "name"
        ) |>
        rename(x = x_loc, y = y_loc) |>
        select(-name)
    })

    dep <- lapply(marg_join_lst, \(x) {
      x1_cond <- NULL
      if (!is.null(cond_var) && cond_var == "X1") {
        x1_cond <- fit_evgam(
          df = x,
          dep_val = cond_val,
          var = "X2",
          cond_var = "X1"
        )
      }

      x2_cond <- NULL
      if (!is.null(cond_var) && cond_var == "X2") {
        x2_cond <- fit_evgam(
          df = x,
          dep_val = cond_val,
          var = "X1",
          cond_var = "X2"
        )
      }

      # join and change to `cecl_dep` format
      coef_evgam <- bind_rows(x1_cond$predictions, x2_cond$predictions) |>
        arrange(name, var, cond_var)
      rownames(coef_evgam) <- NULL
      as_cecl_dep(coef_evgam)
    })
  } else {
    stop("Invalid value for use_evgam")
  }

  # TODO Implement (and add message)
  thresh_max <- NULL
  if (use_dth == TRUE) {
    thresh <- lapply(dep, \(x) {
      lapply(x$dependence, CeCl:::pull_thresh_trans)
    })
    # take max
    thresh_max <- max(unlist(thresh))
  }

  # compute distances
  dist <- lapply(seq_along(dep), \(i) {
    cecl_dist(
      dep[[i]],
      marg[[i]],
      n_mc = n_mc,
      dth = thresh_max,
      ...
    )
  })
  if (ret_dep) {
    return(list("dist" = dist, "dep" = dep))
  } else {
    return(dist)
  }
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
single_run_explore <- \(data, i, laplace_trans = FALSE, ...) {
  data_block <- data %>%
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

  dist_obj <- suppressMessages(calc_dist(data_laplace_block, ...))
  dist <- dist_obj
  # cond <- length(dist_obj) == 2
  cond <- !is.null(names(dist_obj)) && all(names(dist_obj) == c("dist", "dep"))
  if (cond) {
    dist <- dist$dist
  }

  # plot distance matrices
  # ggplot(dist[[1]], which = "image")
  # ggplot(dist[[2]], which = "image")

  norm_val_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
  norm_val_inf <- compare_blocks(dist, type = "norm", norm_type = "M")
  norm_val_spec <- compare_blocks(dist, type = "norm", norm_type = "2")

  if (cond) {
    return(list(
      frob = norm_val_frob,
      inf = norm_val_inf,
      spec = norm_val_spec,
      dep = dist_obj$dep
    ))
  } else {
    return(list(
      frob = norm_val_frob,
      inf = norm_val_inf,
      spec = norm_val_spec
    ))
  }
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
      xintercept = (nrow(data_loc_orig) / n_locs) / 2, # TODO Need to add this as an argument!
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
test_norm_diff <- \(data, block_vec, n_mc = 500, cond_val = 0, ret_all = TRUE, use_start = FALSE, ...) {
  # initial parameter values for CE
  start <- list(
    c(a = 0.01, b = 0.01),
    c(a = 0.01, b = 0.01)
  )
  if (use_start == TRUE) {
    norm1 <- single_run_explore(
      data, block_vec[1],
      n_mc = n_mc, cond_val = cond_val, ret_dep = TRUE, ...
    )
    dep <- norm1$dep
    norm1 <- norm1[!names(norm1) == "dep"]
    start <- lapply(dep, coef)
    block_vec <- block_vec[-1]
  }

  # norm_vals <- mclapply(block_vec, \(i) {
  if (length(block_vec) >= 1) {
    norm_vals <- lapply(block_vec, \(i) {
      single_run_explore(
        data, i,
        n_mc = n_mc, cond_val = cond_val, start = start, ...
      )
      # }, mc.cores = mc_cores)
    })
    # add norm values from first run to the rest
    if (use_start == TRUE) {
      norm_vals <- c(list(norm1), norm_vals)
    }
  } else {
    norm_vals <- list(norm1)
  }

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

# Function to generate data for a single permutation test
perm_test_prep <- \(
  data_laplace_loc,
  cond_prob = 0.9,
  cond_val = NULL,
  cond_var = NULL,
  use_dth = FALSE,
  use_evgam = FALSE,
  data_loc = NULL,
  aLow = -1,
  nruns = 1,
  start = list(c(a = 0.01, b = 0.01), c(a = 0.01, b = 0.01)),
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

  locs <- unique(data_perm$name)

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
  # TODO Change!!
  # dep_perm <- lapply(
  #   marg_perm,
  #   cecl_dep,
  #   cond_prob = cond_prob,
  #   cond_val = cond_val,
  #   fit_no_keef = TRUE
  # )
  # fit CE models
  if (use_evgam == FALSE) {
    # dep_perm <- lapply(
    #   marg_perm,
    #   cecl_dep,
    #   cond_prob   = cond_prob,
    #   cond_val    = cond_val,
    #   fit_no_keef = TRUE,
    #   aLow        = aLow,
    #   nruns       = nruns
    # )
    dep_perm <- lapply(seq_along(marg_perm), \(i) {
      cecl_dep(
        obj         = marg_perm[[i]],
        cond_prob   = cond_prob,
        cond_val    = cond_val,
        cond_var    = cond_var,
        aLow        = aLow,
        nruns       = nruns,
        fit_no_keef = TRUE,
        start       = start[[i]]
      )
    })
    # skip bad permutations
    if (any(sapply(dep_perm, \(x) any(is.na(unlist(x$dependence)))))) {
      return(NA)
    }
  } else if (use_evgam) {
    if (is.null(data_loc)) {
      stop("data_loc must be provided when use_evgam = TRUE")
    }

    # first, put marginal data into form for evgam fitting
    marg_join_lst <- lapply(seq_along(marg_perm), \(k) {
      x <- marg_perm[[k]]
      bind_rows(lapply(seq_along(x$transformed), \(j) {
        data.frame(x$transformed[[j]]) |>
          mutate(name = paste0("loc_", j))
      })) |>
        left_join(
          data_loc |>
            filter(block == k) |>
            distinct(name, x_loc, y_loc),
          by = "name"
        ) |>
        rename(x = x_loc, y = y_loc) |>
        select(-name)
    })

    dep_perm <- lapply(marg_join_lst, \(x) {
      x1_cond <- fit_evgam(
        df = x,
        dep_val = cond_val,
        var = "X2",
        cond_var = "X1"
      )

      x2_cond <- fit_evgam(
        df = x,
        dep_val = cond_val,
        var = "X1",
        cond_var = "X2"
      )

      # join and change to `cecl_dep` format
      coef_evgam <- bind_rows(x1_cond$predictions, x2_cond$predictions) |>
        arrange(name, var, cond_var)
      rownames(coef_evgam) <- NULL
      as_cecl_dep(coef_evgam)
    })
  } else {
    stop("Invalid value for use_evgam")
  }


  thresh_max <- NULL
  if (use_dth == TRUE) {
    thresh <- lapply(dep_perm, \(x) {
      lapply(x$dependence, CeCl:::pull_thresh_trans)
    })
    # take max
    thresh_max <- max(unlist(thresh))
  }

  # compute distances

  dist_perm <- lapply(seq_along(dep_perm), \(i) {
    cecl_dist(
      dep_perm[[i]],
      marg_perm[[i]],
      dth = thresh_max,
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
# TODO Optimise better to work in parallel, not great at the moment
perm_test_fun <- \(
  data,
  grid_vals,
  n_perm = 100,
  n_per_block = NULL,
  laplace_trans = TRUE,
  use_start = FALSE,
  ret_dep = FALSE,
  use_dth = FALSE,
  ...
) {
  # initial parameter values for CE
  start <- list(
    c(a = 0.01, b = 0.01),
    c(a = 0.01, b = 0.01)
  )
  # if use_start is TRUE, get starting values from fitting CE to original data
  if (use_start) {
    data_start <- data
    if (!"block" %in% names(data)) {
      data_start <- data_start |>
        mutate(block = ifelse(row_number() <= n() / 2, "1", "2"))
    }
    data_laplace_start <- trans_fun(data_start, n_vars, laplace_trans)

    dots <- list(...)
    # if (!"cond_val" %in% names(dots)) {
    #   cond_val_start <- qlaplace(0.8)
    # } else {
    #   cond_val_start <- dots$cond_val
    # }
    if (!"cond_val" %in% names(dots)) {
      cond_val_start <- NULL
    } else {
      cond_val_start <- dots$cond_val
    }
    if (!"cond_prob" %in% names(dots)) {
      cond_prob_start <- NULL
    } else {
      cond_prob_start <- dots$cond_prob
    }

    if (!"cond_var" %in% names(dots)) {
      cond_var_start <- NULL
    } else {
      cond_var_start <- dots$cond_var
    }
    if (!"nruns" %in% names(dots)) {
      nruns_start <- 1
    } else {
      nruns_start <- dots$nruns
    }
    if (!"aLow" %in% names(dots)) {
      aLow_start <- -1
    } else {
      aLow_start <- dots$aLow
    }

    dep_start <- lapply(seq_along(data_laplace_start), \(k) {
      cecl_dep(
        data_laplace_start[[k]],
        cond_val = cond_val_start,
        cond_prob = cond_prob_start,
        fit_no_keef = TRUE,
        nruns = nruns_start,
        cond_var = cond_var_start,
        aLow = aLow_start,
        start = start[[k]]
      )
    })

    # start values
    start <- lapply(dep_start, coef)
  }

  # fun inside apply
  ret_fun <- \(i) {
    # browser()
    # Add block based on i
    data_block <- data %>%
      group_by(name) |>
      mutate(block = ifelse(row_number() <= i, "1", "2")) |>
      ungroup()

    # if desired, take n_per_block obs for each location and block
    if (!is.null(n_per_block)) {
      # data_block <- data_block %>%
      #   group_by(name, block) %>%
      #   slice_tail(n = n_per_block) %>%
      #   ungroup()
      data_block <- data_block |>
        group_split(block, .keep = TRUE) |>
        lapply(\(x) {
          y <- x %>%
            group_by(name)
          # for first block, take last n_per_block obs
          if (as.numeric(y$block[1]) == 1) {
            y <- y |>
              slice_tail(n = n_per_block)
            # for second block, take first n_per_block obs
          } else if (as.numeric(y$block[1]) == 2) {
            y <- y |>
              slice_head(n = n_per_block)
          }
          return(ungroup(y))
        }) |>
        bind_rows()
    }

    # check blocking is correct
    # table(data_block$block, data_block$name)

    # convert to Laplace and store marginals
    data_laplace_block <- trans_fun(data_block, n_vars, laplace_trans)

    # compute distances for original data (before permutation)
    # debugonce(cecl_dep)
    dist <- suppressMessages(calc_dist(
      data_laplace_block,
      start = start, ret_dep = ret_dep, use_dth = use_dth, ...
    ))
    if (ret_dep == TRUE) {
      dep <- dist$dep
      dist <- dist$dist
    }

    norm_orig_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
    norm_orig_inf <- compare_blocks(dist, type = "norm", norm_type = "M")

    # compute permutation distances for original data
    norm_vals <- lapply(seq_len(n_perm), \(j) {
      # norm_vals <- mclapply(seq_len(n_perm), \(j) {
      # print progress of permutations
      # sprintf(system("echo %s", j), "Permutation")
      dist <- perm_test_prep(
        data_laplace_block,
        start = start,
        use_dth = use_dth,
        # laplace_sample = laplace_sample,
        ...
      )
      list(
        "frob" = compare_blocks(dist, type = "norm", norm_type = "F"),
        "inf"  = compare_blocks(dist, type = "norm", norm_type = "M")
      )
    })

    # extract frobenius and infinity norms
    perm_norms_frob <- sapply(norm_vals, \(x) x$frob)
    perm_norms_inf <- sapply(norm_vals, \(x) x$inf)

    # compute p-values
    # TODO Are these p-values for one-sided or two-sided tests?
    (p_value_frob <- mean(perm_norms_frob >= norm_orig_frob))
    (p_value_inf <- mean(perm_norms_inf >= norm_orig_inf))

    # print how many iterations completed (as a percentage)
    # print(paste0(
    #   substr(which(grid_vals == i) / length(grid_vals) * 100, 0, 5), "% complete\n"
    # ))
    # system(sprintf(
    #   "%s%% complete\n",
    #   round(which(grid_vals == i) / length(grid_vals) * 100, 2)
    # ))
    system(sprintf(
      'echo "\n%s\n"',
      paste0(round(which(grid_vals == i) / length(grid_vals), 3) * 100, "% completed", collapse = "")
    ))

    ret_list <- list(
      p_value_frob    = p_value_frob,
      p_value_inf     = p_value_inf,
      norm_orig_frob  = norm_orig_frob,
      norm_orig_inf   = norm_orig_inf,
      perm_norms_frob = perm_norms_frob,
      perm_norms_inf  = perm_norms_inf
    )
    if (ret_dep == TRUE) {
      dep_vals_df <- lapply(dep, coef) |>
        bind_rows(.id = "block")
      ret_list <- c(ret_list, list(dep_vals = dep_vals_df))
    }
    return(ret_list)
  }
  return(mclapply(grid_vals, ret_fun))
  # return(lapply(grid_vals, ret_fun))
}

seed_fun <- \(seed) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }

    set.seed(seed)

    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
  }
}

# Function to generate spatial data with varying correlation
gen_spatial_data <- \(
  cor_t_vec,
  n_clust = 3,
  n_locs = 60,
  n_per_loc = 1000,
  sd_loc = 0.8, # spatial spread around cluster centres
  ell = 0.5, # length scale for spatial smoothing (bigger == smoother)
  eps = 1e-6, # to keep rho inside (-1, 1)
  seed = NULL,
  clust_centres = NULL,
  locs = NULL,
  laplace_trans = FALSE
) {
  # generate points for clsuter centres, if not provided
  if (is.null(clust_centres)) {
    # set seed, if desired
    seed_fun(seed)
    clust_centres <- data.frame(
      cor = cor_t_vec,
      x = runif(n_clust, 0, 5),
      y = runif(n_clust, 0, 5)
    )
  }

  # plot
  # ggplot(clust_centres, aes(x = x, y = y, color = as.factor(cor))) +
  #   geom_point(size = 4) +
  #   labs(color = "Correlation") +
  #   theme_minimal()

  # generate locations around cluster centres

  # smooth correlation function for these points, while keeping centres fixed

  # how many locations per cluster
  locs_per_clust <- n_locs / n_clust
  stopifnot(locs_per_clust == floor(locs_per_clust))
  if (is.null(locs)) {
    # generate locations around cluster centres, if not provided
    seed_fun(seed)
    locs <- clust_centres %>%
      mutate(cluster = row_number()) %>%
      tidyr::uncount(locs_per_clust) %>% # replicate centres locs_per_clust times
      group_by(cor, x, y) %>%
      # generate locations around centre with some noise
      mutate(
        loc_id = row_number(),
        x_loc  = x + rnorm(n(), 0, sd_loc),
        y_loc  = y + rnorm(n(), 0, sd_loc)
      ) %>%
      ungroup() %>%
      mutate(loc_global_id = row_number()) %>%
      select(cluster, loc_global_id, x_loc, y_loc)
  }

  # head(locs)
  #
  # ggplot() +
  #   geom_point(data = locs, aes(x = x_loc, y = y_loc), alpha = 0.6) +
  #   geom_point(data = clust_centres, aes(x = x, y = y, colour = as.factor(cor)), size = 4) +
  #   theme_minimal() +
  #   labs(colour = "Centre corr")

  # length scale (bigger => smoother, more overlap among centres)
  # ell <- 2.0

  # distance matrix: rows = locations, cols = centres
  D <- sp::spDists(
    x = as.matrix(locs[, c("x_loc", "y_loc")]),
    y = as.matrix(clust_centres[, c("x", "y")]),
    longlat = FALSE
  )

  # kernel weights + normalize to sum to 1 per location
  W <- exp(-(D^2) / (2 * ell^2))
  W <- W / rowSums(W)

  # blended rho for each location
  rho_loc <- as.numeric(W %*% clust_centres$cor)

  # calculate for each cluster and return
  rho_df <- data.frame(
    "rho" = rho_loc,
    # TODO Change this line, as locs may be provided externally
    "cluster" = rep(1:n_clust, each = locs_per_clust)
  )

  rho_df |>
    group_by(cluster) |>
    summarise(mean(rho)) |>
    print()

  # keep rho inside (-1, 1) for copula validity (esp if you later allow negatives)
  # eps <- 1e-6
  rho_loc <- pmin(1 - eps, pmax(-1 + eps, rho_loc))

  locs$rho <- rho_loc

  # ggplot(locs, aes(x = x_loc, y = y_loc, colour = rho)) +
  #   geom_point(size = 2) +
  #   scale_colour_viridis_c() +
  #   theme_minimal() +
  #   labs(colour = "rho(x,y)")

  # join centres to locations for plotting
  cent_as_locs <- clust_centres %>%
    transmute(
      loc_global_id = 0L + row_number(),
      x_loc = x, y_loc = y, rho = cor, is_centre = TRUE
    ) |>
    mutate(cluster = row_number())

  locs$is_centre <- FALSE
  locs2 <- bind_rows(cent_as_locs, locs) %>%
    arrange(loc_global_id)

  # generate data at each location
  data <- bind_rows(lapply(1:n_locs, function(ell) {
    x <- gen_t(rho_loc[ell], n_per_loc, n_vars = n_vars, df_t = df_t, laplace_trans)
    x$name <- paste0("loc_", ell)
    x
  }))

  # join in spatial data
  data_spat <- data |>
    left_join(
      locs %>%
        mutate(
          name = paste0("loc_", loc_global_id)
        ) |>
        select(name, cluster, x_loc, y_loc, rho)
    )

  return(data_spat)
}

# Function to fit CE model using evgam for a given conditioning variable
# TODO Allow cond_prob to be specified here as well!
fit_evgam <- \(df, dep_val, var, cond_var, k = 10) {
  df_ev <- subset(df, get(cond_var) > dep_val)

  # check how many unique locations left
  n_orig <- length(unique(paste(df$x, df$y)))
  n_left <- length(unique(paste(df_ev$x, df_ev$y)))
  if (n_left < n_orig) {
    message(paste0(
      "After thresholding on ", cond_var, " only ", n_left,
      " unique locations left (originally ", n_orig, ")"
    ))
  }

  # fit evgam model
  fmla_condex <- list(
    as.formula(paste0(var, " ~ s(x, y, k = ", k, ")")),
    # ~1,
    as.formula(paste0(var, " ~ s(x, y, k = ", k, ")")),
    as.formula(paste0("~ s(x, y, k = ", k, ")")),
    as.formula(paste0("~ s(x, y, k = ", k, ")"))
  )

  message(paste("Fitting evgam for var", var, "conditioning on", cond_var))
  m1 <- evgam(
    fmla_condex,
    data   = df_ev,
    family = "condex",
    args   = list(x = matrix(df_ev[[cond_var]], ncol = 1))
  )
  message("Model fitted")

  # predict for locations in data
  grid_df <- data.frame(x = df_ev$x, y = df_ev$y) |>
    distinct()
  pred_df <- cbind(grid_df, predict(m1, newdata = grid_df, type = "response"))
  # match to coef.cecl_dep format
  pred_df$var <- var
  pred_df$cond_var <- cond_var
  return(list(
    "model" = m1,
    "predictions" = pred_df |>
      mutate(
        name     = paste0("loc_", 1:nrow(pred_df)),
        dth      = dep_val,
        ll       = NA
      ) |>
      rename(
        a = alpha,
        b = beta,
        m = location,
        s = scale
      ) |>
      select(name, var, cond_var, a, b, m, s, ll, dth)
  ))
}

#### Functions ####

# function for plotting alpha, beta estimates
plot_ab_map <- \(
  ab_df, data, areas, n_breaks = 8, range = c(1, 6), elev_df = NULL
) {
  ab_df_long <- ab_df
  if (!"parameter" %in% names(ab_df)) {
    ab_df_long <- ab_df |>
      tidyr::pivot_longer(a:dth, names_to = "parameter")
  }

  ab_sf <- ab_df_long |>
    left_join(distinct(data, name, lon, lat)) |>
    # st_to_sf()
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas))

  # check that b values aren't exceptionally negative (<-1)
  # also check that all values are equal to fixed value for b, if applying
  # sort(ab_sf[ab_sf$parameter == "b", ]$value)[1:10]

  # plot parameter values for each parameter and variable
  names <- c("a", "b")
  p_lst <- lapply(seq_along(names), \(i) {
    # browser()
    p <- ggplot()
    # add elevation raster
    if (!is.null(elev_df)) {
      # plot bins if available, if not plot directly (less distinct colours)
      if ("elev_bin" %in% names(elev_df)) {
        p <- p +
          geom_tile(
            data = elev_df,
            aes(x = x, y = y, fill = elev_bin),
            width = diff(range(elev_df$x)) / length(unique(elev_df$x)),
            height = diff(range(elev_df$y)) / length(unique(elev_df$y)),
            show.legend = FALSE
          ) +
          labs(x = "", y = "") +
          # scale_fill_viridis_d(
          #   name = "Elevation",
          #   option = "D",
          #   direction = 1
          # )
          scico::scale_fill_scico_d(
            palette = "bilbao",
            name = "Elevation\n(m)",
            direction = -1
          )
      } else {
        p <- p +
          geom_tile(
            data = elev_df,
            aes(x = x, y = y, fill = elevation),
            show.legend = FALSE
          ) +
          labs(x = "", y = "") +
          scale_fill_viridis(
            name = "Elevation\n(m)",
            option = "D",
            limits = c(0, 1500)
          )
      }
      p <- p + geom_sf(data = areas, fill = NA, colour = "black")
    } else {
      p <- p + geom_sf(data = areas, fill = NA, colour = "black")
    }
    p +
      # geom_sf(data = areas, fill = NA, colour = "white") +
      coord_sf(expand = FALSE) + # remove padding around plot
      ggnewscale::new_scale_fill() +
      geom_sf(
        data = ab_sf %>%
          filter(parameter == names[i]) %>%
          # mutate(vars = paste0(parameter, " - ", vars)),
          mutate(vars = paste0(parameter, " - ", cond_var)),
        aes(fill = value, size = value),
        colour = "black",
        stroke = 1,
        pch = 21
      ) +
      scale_size_continuous(
        breaks = scales::extended_breaks(n = n_breaks),
        range = range,
        guide = "legend"
      ) +
      scale_fill_gradient2(
        low = "blue3",
        high = "red3",
        na.value = "grey",
        breaks = scales::extended_breaks(n = n_breaks),
        guide = "legend"
      ) +
      guides(fill = guide_legend(), size = guide_legend()) +
      labs(fill = "", size = "") +
      facet_wrap(
        ~vars,
        ncol = 2,
        # labeller = as_labeller(c(
        #   "a - rain"       = "alpha ~ ' - ' ~ 'Precipitation | Wind Speed'",
        #   "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Precipitation'",
        #   "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Precipitation'",
        #   "b - rain"       = "beta ~ ' - ' ~ 'Precipitation | Wind Speed'",
        #   "b - wind_speed" = "beta ~ ' - ' ~ 'Wind Speed | Precipitation'"
        # ), default = label_parsed)
      ) +
      cecl_theme(nejm_pal = FALSE)
  })
  return(p_lst)
}


#' @title Plot clustering solution on map
#' @description Plot clustering solution on map
#' @param pts Spatial points object
#' @param areas Spatial polygons object
#' @param clust_obj Clustering object
#' @return ggplot object
#' @rdname plt_clust_map
#' @export
# TODO: Could make this plot/ggplot method for object
plt_clust_map <- \(
  pts,
  areas,
  clust_obj,
  plot_medoids = TRUE,
  elev_df = NULL,
  rm_elev_leg = TRUE,
  pt_size = 4
) {
  name <- clust <- medoid <- NULL

  if (inherits(clust_obj, "kmeans")) {
    clust_element <- "cluster"
    medoid_locs <- NA
  } else if (inherits(clust_obj, "pam")) {
    clust_element <- "clustering"
    medoids <- clust_obj$medoids
    medoid_locs <- NA
    if (inherits(clust_obj$medoids, "character")) {
      medoid_locs <- medoids
    } else if (!is.null(rownames(medoids))) {
      medoid_locs <- rownames(medoids)
    }
  } else {
    stop("Clustering class not currently supported")
  }

  # reorder alphabetically
  # TODO: Look into this, required??
  # clust_names <- names(clust_obj[[clust_element]])
  # if (!is.null(clust_names)) {
  #   clust_obj[[clust_element]] <- clust_obj[[clust_element]][order(clust_names)]
  # }

  # TODO: Medoids doesn#t work as rows are named, fix!
  # pts_plt <- cbind(pts, data.frame("clust" = clust_obj[[clust_elementement]])) |>
  #   dplyr::mutate(
  #     medoid = ifelse(name %in% medoid_locs, TRUE, FALSE),
  #     medoid = factor(medoid, levels = c(FALSE, TRUE))
  #   )

  # extract cluster membership for each site
  clust_df <- dplyr::tibble(
    "name" = names(clust_obj[[clust_element]]),
    "clust" = clust_obj[[clust_element]]
  )
  # join to location data
  pts_plt <- pts |>
    dplyr::left_join(clust_df, by = "name") |>
    dplyr::mutate(
      medoid = ifelse(name %in% medoid_locs, TRUE, FALSE),
      medoid = factor(medoid, levels = c(FALSE, TRUE))
    )

  point_cols <- ggsci::pal_nejm()(sum(!is.na(unique(pts_plt$clust))))
  any_na <- FALSE
  if (any(is.na(pts_plt$clust))) {
    any_na <- TRUE
    #   pts_plt <- pts_plt |>
    #     dplyr::mutate(
    #       clust = ifelse(is.na(clust), 4, clust)
    #     )
    #
    #   point_cols <- c(point_cols, "white")
    pts_plt_na <- pts_plt |>
      dplyr::filter(is.na(clust))
    pts_plt <- pts_plt |>
      dplyr::filter(!is.na(clust))
  }

  # plot locations on map, colouring by cluster
  p <- ggplot2::ggplot()

  if (!is.null(elev_df)) {
    # plot bins if available, if not plot directly (less distinct colours)
    if ("elev_bin" %in% names(elev_df)) {
      p <- p + geom_tile(
        data = elev_df,
        aes(x = x, y = y, fill = elev_bin),
        width = diff(range(elev_df$x)) / length(unique(elev_df$x)),
        height = diff(range(elev_df$y)) / length(unique(elev_df$y))
      ) +
        # scale_fill_viridis_d(
        #   name = "Elevation",
        #   option = "D",
        #   direction = 1
        # )
        scico::scale_fill_scico_d(
          name = "Elevation",
          # palette = "lajolla",
          palette = "bilbao",
          direction = -1
        )
      # scale_fill_manual(
      #   values = c(
      #     "#FFFFFF", "#C5C2B2", "#B19E68",
      #     "#A6785B", "#9B5352", "#6D1F23"
      #     # "#FFFFFF", "#CBC9C0", "#BBB287",
      #     # "#AC8F60", "#A4745A", "#914249"
      #   ),
      #   name = "Elevation"
      # )
    } else {
      p <- p +
        geom_tile(
          data = elev_df,
          aes(x = x, y = y, fill = elevation)
        ) +
        scale_fill_viridis(
          name = "Elevation\n(m)",
          option = "D",
          limits = c(0, 1500)
        )
    }

    # remove legend for elevation; handy if joining plots with patchwork!
    if (rm_elev_leg) {
      p <- p +
        guides(fill = "none")
    }

    # p <- p + geom_sf(data = areas, fill = NA, colour = "white")
    # } else {
    #   p <- p + geom_sf(data = areas, fill = NA, colour = "black")
  }
  p <- p +
    geom_sf(data = areas, fill = NA, colour = "black") +
    coord_sf(expand = FALSE) + # remove padding around plot
    ggnewscale::new_scale_fill() +
    # ggnewscale::new_scale_colour() +
    NULL
  # p <- ggplot2::ggplot(areas) +
  # p <- p +
  # ggplot2::geom_sf(colour = "black", fill = NA) +
  # ggplot2::geom_sf(areas, colour = "black", fill = NA)
  # ggplot2::geom_sf(
  #   data = pts_plt,
  #   ggplot2::aes(
  #     colour = factor(clust), shape = medoid, size = as.numeric(medoid)
  #   ),
  #   alpha = 0.8
  # ) +

  if (plot_medoids == TRUE) {
    p <- p +
      ggplot2::geom_sf(
        data = pts_plt,
        ggplot2::aes(
          # colour = factor(clust),
          fill = factor(clust),
          shape = medoid,
          size = as.numeric(medoid)
        ),
        alpha = 0.8,
        # alpha = 0.9,
        stroke = 1,
        pch = 21
      ) +
      ggplot2::scale_shape_discrete(breaks = c(1, 15)) +
      # ggplot2::scale_size_continuous(range = c(4, 8)) +
      ggplot2::scale_size_continuous(range = c(pt_size, pt_size * 2)) +
      ggplot2::guides(shape = "none", size = "none")
  } else {
    p <- p +
      ggplot2::geom_sf(
        data = pts_plt,
        # ggplot2::aes(colour = factor(clust)),
        ggplot2::aes(fill = factor(clust)),
        alpha = 0.8,
        # alpha = 0.9,
        # size = 4,
        size = pt_size,
        stroke = 1,
        pch = 21
      )
  }

  # Add white (i.e. missing) pts on top, to make sure they're visible
  if (any_na) {
    p <- p +
      ggplot2::geom_sf(
        data = pts_plt_na |>
          dplyr::mutate(clust = 1),
        # ggplot2::aes(colour = factor(clust)),
        # ggplot2::aes(fill = factor(clust)),
        fill = "white",
        alpha = 0.8,
        # alpha = 0.9,
        # size = 4,
        size = pt_size,
        stroke = 1,
        pch = 21
      )
  }

  p <- p +
    # ggplot2::labs(colour = "Cluster") +
    # ggplot2::labs(fill = "Cluster") +
    guides(fill = "none") +
    coord_sf(expand = FALSE) + # remove padding around plot
    # evc_theme(nejm_pal = FALSE) +
    cecl_theme(nejm_pal = FALSE) +
    # ggsci::scale_colour_nejm()
    # ggsci::scale_fill_nejm() +
    scale_fill_manual(
      values = point_cols
    ) +
    labs(x = "", y = "")

  return(p)
}

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
      mutate(block = factor(block, levels = unique(block))) |>
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

# function to calculate log-Euclidean SPD distance between two distance matrices
# log_euclid <- \(M1, M2, sigma = 1) {
# log_euclid <- \(M1, M2, sigma = NULL, eps = 1e-6) {
#   if (is.null(sigma)) {
#     # estimate sigma as median of upper triangular values of M1 and M2
#     vals <- c(M1[upper.tri(M1)], M2[upper.tri(M2)])
#     sigma <- median(vals[vals > 0])
#   }
#
#   # apply Gaussian kernel to convert to similarity matrices
#   K1 <- exp(-M1 / sigma)
#   K2 <- exp(-M2 / sigma)
#
#   # add small value to diagonal to ensure positive definiteness
#   K1 <- K1 + diag(eps, nrow(K1))
#   K2 <- K2 + diag(eps, nrow(K2))
#
#   # compute matrix logarithm and then Frobenius norm of difference
#   L1 <- expm::logm(K1)
#   L2 <- expm::logm(K2)
#
#   norm(L1 - L2, type = "F")
# }
log_euclid <- \(M1, M2, sigma, tol = 1e-10) {
  M1 <- (M1 + t(M1)) / 2
  M2 <- (M2 + t(M2)) / 2

  K1 <- exp(-(M1 / sigma)^2 / 2)
  K2 <- exp(-(M2 / sigma)^2 / 2)

  K1 <- (K1 + t(K1)) / 2
  K2 <- (K2 + t(K2)) / 2

  eig1 <- eigen(K1, symmetric = TRUE)
  eig2 <- eigen(K2, symmetric = TRUE)

  if (
    min(eig1$values) <= tol ||
      min(eig2$values) <= tol
  ) {
    stop("Kernel matrix is not strictly positive definite.")
  }

  L1 <- eig1$vectors %*%
    diag(log(eig1$values)) %*%
    t(eig1$vectors)

  L2 <- eig2$vectors %*%
    diag(log(eig2$values)) %*%
    t(eig2$vectors)

  norm(L1 - L2, type = "F")
}

# Function to calculate norm or clustering diff
compare_blocks <- \(
  dist,
  type = c("norm", "clustering"),
  norm_type = "F"
) {
  type <- match.arg(type)

  if (type == "norm") {
    M1 <- as.matrix(dist[[1]]$dist_mat)
    M2 <- as.matrix(dist[[2]]$dist_mat)

    if (
      !is.null(rownames(M1)) &&
        !is.null(rownames(M2))
    ) {
      if (
        !setequal(rownames(M1), rownames(M2)) ||
          !setequal(colnames(M1), colnames(M2))
      ) {
        stop(
          "The two distance matrices contain different locations."
        )
      }

      M2 <- M2[
        rownames(M1),
        colnames(M1),
        drop = FALSE
      ]
    }

    if (!identical(dim(M1), dim(M2))) {
      stop("Distance matrices have different dimensions.")
    }

    if (norm_type == "log SPD") {
      return(log_euclid(M1, M2))
    }

    A <- M2 - M1
    value <- norm(A, type = norm_type)

    if (norm_type == "2") {
      frobenius <- norm(A, type = "F")

      if (
        value >
          frobenius +
            100 * .Machine$double.eps *
              max(1, frobenius)
      ) {
        stop(
          "Spectral norm exceeds Frobenius norm: ",
          "spectral = ", value,
          ", Frobenius = ", frobenius
        )
      }
    }

    return(value)
  }

  if (type == "clustering") {
    clust <- lapply(
      dist,
      cecl_clust,
      k = n_clust,
      seed = seed
    )

    return(
      mclust::adjustedRandIndex(
        clust[[1]]$pam$clustering,
        clust[[2]]$pam$clustering
      )
    )
  }
}

# compare_blocks <- \(dist, type = c("norm", "clustering"), norm_type = "F") {
#   if (type == "norm") {
#     # norm of difference
#     if (norm_type == "log SPD") {
#       return(
#         log_euclid(as.matrix(dist[[1]]$dist_mat), as.matrix(dist[[2]]$dist_mat))
#       )
#     } else {
#       mat <- as.matrix(dist[[2]]$dist_mat - dist[[1]]$dist_mat)
#       return(norm(
#         mat,
#         type = norm_type
#       ))
#     }
#   } else if (type == "clustering") {
#     # clustering solution
#     clust <- lapply(
#       dist,
#       cecl_clust,
#       k = n_clust,
#       # cluster_mem = sort(rep(1:n_clust, each = n_locs / n_clust)),
#       seed = seed
#     )
#     # calculate difference in clustering solutions
#     return(mclust::adjustedRandIndex(
#       clust[[1]]$pam$clustering,
#       clust[[2]]$pam$clustering
#     ))
#   }
# }

# # Function to run single iteration of exploration for a given block position i
# single_run_explore <- \(data, i, laplace_trans = FALSE, ...) {
#   data_block <- data %>%
#     group_by(name) |>
#     mutate(block = ifelse(row_number() <= i, "1", "2")) |>
#     ungroup()
#
#   # check blocking is correct
#   # table(data_block$block, data_block$name)
#
#   # convert to Laplace and store marginals
#   # TODO Do before blocking? Think I do with laplace_trans = TRUE!
#   data_laplace_block <- trans_fun(data_block, n_vars, laplace_trans)
#
#   # check you've transformed correctly
#   # par(mfrow = c(1, 2))
#   # plot(data_laplace_block[[1]][[1]][[1]])
#   # plot(data_laplace_block[[2]][[1]][[1]])
#   # par(mfrow = c(1, 1))
#
#   dist_obj <- suppressMessages(calc_dist(data_laplace_block, ...))
#   dist <- dist_obj
#   # cond <- length(dist_obj) == 2
#   cond <- !is.null(names(dist_obj)) && all(names(dist_obj) == c("dist", "dep"))
#   if (cond) {
#     dist <- dist$dist
#   }
#
#   # plot distance matrices
#   # ggplot(dist[[1]], which = "image")
#   # ggplot(dist[[2]], which = "image")
#
#   mat1 <- dist[[1]]$dist_mat
#
#   norm_val_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
#   # for 2 x 2 matrix, norms are all the same, so just check values
#   dim_cond <- dim(mat1)[[1]] > 2
#   if (dim_cond) {
#     norm_val_inf <- compare_blocks(dist, type = "norm", norm_type = "M")
#     norm_val_spec <- compare_blocks(dist, type = "norm", norm_type = "2")
#     norm_val_ln_spd <- compare_blocks(dist, type = "norm", norm_type = "log SPD")
#   }
#
#   if (cond) {
#     ret <- list(
#       frob = norm_val_frob,
#       # inf = norm_val_inf,
#       # spec = norm_val_spec,
#       # ln_spd = norm_val_ln_spd,
#       dep = dist_obj$dep
#     )
#
#     if (dim_cond) {
#       ret <- list(
#         frob = norm_val_frob,
#         inf = norm_val_inf,
#         spec = norm_val_spec,
#         ln_spd = norm_val_ln_spd,
#         dep = dist_obj$dep
#       )
#     }
#     return(ret)
#   } else {
#     ret <- list(frob = norm_val_frob)
#     if (dim_cond) {
#       ret <- list(
#         frob = norm_val_frob,
#         inf = norm_val_inf,
#         spec = norm_val_spec,
#         ln_spd = norm_val_ln_spd
#       )
#     }
#     return(ret)
#   }
# }

# Calculate discrepancy statistics at one candidate boundary without
# performing permutation tests.
#
# Modes:
#   1. n_per_block supplied:
#      i is a row position and fixed observation windows are used.
#
#   2. n_years_per_block supplied:
#      i is the season_year immediately before the candidate change.
#
#   3. Both NULL:
#      i is a row position and all observations on either side are used.
single_run_explore <- \(
  data,
  i,
  n_per_block = NULL,
  n_years_per_block = NULL,
  n_vars = 2,
  laplace_trans = FALSE,
  skip_failed = TRUE,
  ...
) {
  #### Validate arguments ####

  if (!is.null(n_per_block) && !is.null(n_years_per_block)) {
    stop(
      "Supply only one of `n_per_block` and ",
      "`n_years_per_block`."
    )
  }

  row_window_mode <- !is.null(n_per_block)
  year_window_mode <- !is.null(n_years_per_block)

  if (row_window_mode) {
    if (
      length(n_per_block) != 1L ||
        is.na(n_per_block) ||
        n_per_block < 1L
    ) {
      stop("`n_per_block` must be a positive integer.")
    }

    n_per_block <- as.integer(n_per_block)
  }

  if (year_window_mode) {
    required_columns <- c(
      "name",
      "date",
      "season_year"
    )

    missing_columns <- setdiff(
      required_columns,
      names(data)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Year-based windows require the following columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    if (
      length(n_years_per_block) != 1L ||
        is.na(n_years_per_block) ||
        n_years_per_block < 1L
    ) {
      stop(
        "`n_years_per_block` must be a positive integer."
      )
    }

    n_years_per_block <- as.integer(
      n_years_per_block
    )
  }

  #### Construct blocks ####

  if (year_window_mode) {
    all_years <- sort(unique(data$season_year))
    year_position <- match(i, all_years)

    if (is.na(year_position)) {
      stop(
        "Candidate year ", i,
        " is not present in `season_year`."
      )
    }

    if (year_position < n_years_per_block) {
      stop(
        "Candidate year ", i,
        " does not have ", n_years_per_block,
        " years available on the left."
      )
    }

    if (
      length(all_years) - year_position <
        n_years_per_block
    ) {
      stop(
        "Candidate year ", i,
        " does not have ", n_years_per_block,
        " years available on the right."
      )
    }

    left_years <- all_years[
      (year_position - n_years_per_block + 1L):
      year_position
    ]

    right_years <- all_years[
      (year_position + 1L):
      (year_position + n_years_per_block)
    ]

    data_block <- data |>
      filter(
        season_year %in% c(
          left_years,
          right_years
        )
      ) |>
      mutate(
        block = if_else(
          season_year %in% left_years,
          "1",
          "2"
        )
      ) |>
      arrange(name, date)

    # Confirm that no seasonal year has been split.
    year_check <- data_block |>
      distinct(season_year, block) |>
      count(
        season_year,
        name = "n_blocks"
      )

    if (any(year_check$n_blocks != 1L)) {
      stop(
        "At least one seasonal year was assigned to ",
        "more than one block."
      )
    }

    # Confirm equal date coverage across stations within each block.
    station_block_sizes <- data_block |>
      distinct(name, block, date) |>
      count(
        name,
        block,
        name = "n_dates"
      )

    coverage_check <- station_block_sizes |>
      group_by(block) |>
      summarise(
        n_distinct_sizes = n_distinct(n_dates),
        .groups = "drop"
      )

    if (any(coverage_check$n_distinct_sizes != 1L)) {
      stop(
        "Stations have different numbers of dates within ",
        "at least one block."
      )
    }

    n_left <- data_block |>
      filter(block == "1") |>
      summarise(n = n_distinct(date)) |>
      pull(n)

    n_right <- data_block |>
      filter(block == "2") |>
      summarise(n = n_distinct(date)) |>
      pull(n)

    block_info <- list(
      mode = "season_year",
      split_after = i,
      split_before = all_years[
        year_position + 1L
      ],
      left_years = left_years,
      right_years = right_years,
      n_left = n_left,
      n_right = n_right
    )
  } else {
    #### Original row-based modes ####

    data_ordered <- data

    if ("date" %in% names(data_ordered)) {
      data_ordered <- data_ordered |>
        arrange(name, date)
    }

    n_by_station <- data_ordered |>
      count(
        name,
        name = "n_observations"
      )

    if (i < 1L) {
      stop(
        "The candidate row position must be positive."
      )
    }

    if (any(i >= n_by_station$n_observations)) {
      stop(
        "Candidate row ", i,
        " does not leave observations on both sides ",
        "for every station."
      )
    }

    data_block <- data_ordered |>
      group_by(name) |>
      mutate(
        block = if_else(
          row_number() <= i,
          "1",
          "2"
        )
      ) |>
      ungroup()

    if (row_window_mode) {
      if (i < n_per_block) {
        stop(
          "Candidate row ", i,
          " leaves fewer than ", n_per_block,
          " observations on the left."
        )
      }

      if (
        any(
          n_by_station$n_observations - i <
            n_per_block
        )
      ) {
        stop(
          "Candidate row ", i,
          " leaves fewer than ", n_per_block,
          " observations on the right."
        )
      }

      left_block <- data_block |>
        filter(block == "1") |>
        group_by(name) |>
        slice_tail(n = n_per_block) |>
        ungroup()

      right_block <- data_block |>
        filter(block == "2") |>
        group_by(name) |>
        slice_head(n = n_per_block) |>
        ungroup()

      data_block <- bind_rows(
        left_block,
        right_block
      )

      block_info <- list(
        mode = "fixed_rows",
        split_after = i,
        split_before = i + 1L,
        left_years = NULL,
        right_years = NULL,
        n_left = n_per_block,
        n_right = n_per_block
      )
    } else {
      block_info <- list(
        mode = "cumulative_rows",
        split_after = i,
        split_before = i + 1L,
        left_years = NULL,
        right_years = NULL,
        n_left = i,
        n_right = min(
          n_by_station$n_observations - i
        )
      )
    }
  }

  #### Helper for failed fits ####

  failed_result <- \(
    message,
    stage = "calc_dist"
  ) {
    list(
      success = FALSE,
      block_info = block_info,
      frob = NA_real_,
      inf = NA_real_,
      spec = NA_real_,
      # ln_spd = NA_real_,
      dep = NULL,
      error_stage = stage,
      error_message = message
    )
  }

  #### Convert to required CE structure ####

  data_laplace_block <- trans_fun(
    data_block |>
      select(
        X1,
        X2,
        name,
        block
      ),
    n_vars,
    laplace_trans
  )

  #### Calculate divergence matrices ####

  dist_attempt <- tryCatch(
    {
      suppressMessages(
        calc_dist(
          data_laplace_block,
          ...
        )
      )
    },
    error = \(e) {
      structure(
        list(
          message = conditionMessage(e)
        ),
        class = "calc_dist_error"
      )
    }
  )

  # calc_dist threw an error.
  if (inherits(dist_attempt, "calc_dist_error")) {
    error_message <- paste0(
      "`calc_dist()` failed at candidate ",
      i,
      ": ",
      dist_attempt$message
    )

    if (!skip_failed) {
      stop(error_message, call. = FALSE)
    }

    warning(error_message, call. = FALSE)

    return(
      failed_result(
        message = dist_attempt$message,
        stage = "calc_dist"
      )
    )
  }

  dist_obj <- dist_attempt

  # calc_dist returned a bare NA.
  is_na_result <- (
    is.atomic(dist_obj) &&
      length(dist_obj) == 1L &&
      is.na(dist_obj)
  )

  if (is_na_result) {
    error_message <- paste0(
      "`calc_dist()` returned NA at candidate ",
      i,
      "."
    )

    if (!skip_failed) {
      stop(error_message, call. = FALSE)
    }

    warning(error_message, call. = FALSE)

    return(
      failed_result(
        message = error_message,
        stage = "calc_dist"
      )
    )
  }

  #### Extract distance matrices ####

  has_dep <- (
    !is.null(names(dist_obj)) &&
      all(
        c("dist", "dep") %in%
          names(dist_obj)
      )
  )

  if (has_dep) {
    dist <- dist_obj$dist
  } else {
    dist <- dist_obj
  }

  #### Validate distance-matrix structure ####

  valid_dist_structure <- (
    is.list(dist) &&
      length(dist) == 2L &&
      all(
        vapply(
          dist,
          \(x) {
            is.list(x) &&
              !is.null(x$dist_mat) &&
              # is.matrix(x$dist_mat)
              inherits(x$dist_mat, "dist")
          },
          logical(1)
        )
      )
  )

  if (!valid_dist_structure) {
    error_message <- paste0(
      "Invalid distance-matrix result at candidate ",
      i,
      "."
    )

    if (!skip_failed) {
      stop(error_message, call. = FALSE)
    }

    warning(error_message, call. = FALSE)

    return(
      failed_result(
        message = error_message,
        stage = "validate_dist"
      )
    )
  }

  #### Ensure no location was silently dropped ####

  expected_n_locations <- n_distinct(
    data_block$name
  )

  correct_dimensions <- all(
    vapply(
      dist,
      \(x) {
        identical(
          as.integer(dim(x$dist_mat)),
          c(
            expected_n_locations,
            expected_n_locations
          )
        )
      },
      logical(1)
    )
  )

  if (!correct_dimensions) {
    actual_dimensions <- paste(
      vapply(
        dist,
        \(x) {
          paste(
            dim(x$dist_mat),
            collapse = " x "
          )
        },
        character(1)
      ),
      collapse = "; "
    )

    error_message <- paste0(
      "At least one location was dropped at candidate ",
      i,
      ". Expected two ",
      expected_n_locations,
      " x ",
      expected_n_locations,
      " matrices; obtained ",
      actual_dimensions,
      "."
    )

    if (!skip_failed) {
      stop(error_message, call. = FALSE)
    }

    warning(error_message, call. = FALSE)

    return(
      failed_result(
        message = error_message,
        stage = "validate_dist"
      )
    )
  }

  #### Calculate discrepancy statistics ####

  norm_attempt <- tryCatch(
    {
      mat1 <- dist[[1]]$dist_mat

      norm_results <- list(
        frob = compare_blocks(
          dist,
          type = "norm",
          norm_type = "F"
        ),
        inf = NA_real_,
        spec = NA_real_# ,
        # ln_spd = NA_real_
      )

      if (nrow(mat1) > 2L) {
        norm_results$inf <- compare_blocks(
          dist,
          type = "norm",
          norm_type = "M"
        )

        norm_results$spec <- compare_blocks(
          dist,
          type = "norm",
          norm_type = "2"
        )

        # norm_results$ln_spd <- compare_blocks(
        #   dist,
        #   type = "norm",
        #   norm_type = "log SPD"
        # )
      }

      norm_results
    },
    error = \(e) {
      structure(
        list(
          message = conditionMessage(e)
        ),
        class = "norm_error"
      )
    }
  )

  if (inherits(norm_attempt, "norm_error")) {
    error_message <- paste0(
      "Norm calculation failed at candidate ",
      i,
      ": ",
      norm_attempt$message
    )

    if (!skip_failed) {
      stop(error_message, call. = FALSE)
    }

    warning(error_message, call. = FALSE)

    return(
      failed_result(
        message = norm_attempt$message,
        stage = "compare_blocks"
      )
    )
  }

  #### Return successful result ####

  list(
    success = TRUE,
    block_info = block_info,
    frob = norm_attempt$frob,
    inf = norm_attempt$inf,
    spec = norm_attempt$spec,
    # ln_spd = norm_attempt$ln_spd,
    dep = if (has_dep) {
      dist_obj$dep
    } else {
      NULL
    },
    error_stage = NA_character_,
    error_message = NA_character_
  )
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

# Generate divergence matrices for one permutation.
#
# `permutation_groups` should be either:
#
#   NULL:
#     Each observation/date is permuted separately, but the same
#     permutation is applied at every location.
#
#   A list of length two:
#     One vector for each original block. Each vector must have one
#     value per observation/date in that block. Repeated values identify
#     observations that must remain together, e.g. season_year.
#
# When season-year vectors are supplied, whole seasonal years are
# assigned to pseudo-blocks.

perm_test_prep <- \(
  data_laplace_loc,
  permutation_groups = NULL,
  cond_prob = 0.9,
  cond_val = NULL,
  cond_var = NULL,
  use_dth = FALSE,
  use_evgam = FALSE,
  data_loc = NULL,
  aLow = -1,
  nruns = 1,
  n_mc = 500,
  start = list(
    c(a = 0.01, b = 0.01),
    c(a = 0.01, b = 0.01)
  ),
  ...
) {
  #### Validate input ####
  
  if (length(data_laplace_loc) != 2L) {
    stop(
      "`data_laplace_loc` must contain exactly two blocks."
    )
  }
  
  if (use_evgam) {
    stop(
      "The rewritten permutation function does not currently ",
      "support `use_evgam = TRUE`, because the spatial metadata ",
      "must also be reassigned to the permuted blocks."
    )
  }
  
  if (!is.null(cond_val)) {
    cond_prob <- NULL
  }
  
  transformed_blocks <- lapply(
    data_laplace_loc,
    \(x) x$transformed
  )
  
  location_names <- names(
    transformed_blocks[[1]]
  )
  
  if (is.null(location_names)) {
    stop(
      "The transformed matrices in block 1 are not named by location."
    )
  }
  
  if (
    !setequal(
      location_names,
      names(transformed_blocks[[2]])
    )
  ) {
    stop(
      "The two blocks contain different locations."
    )
  }
  
  # Put locations into identical order in both blocks.
  transformed_blocks <- lapply(
    transformed_blocks,
    \(x) x[location_names]
  )
  
  #### Check row counts within each block ####
  
  block_sizes <- vapply(
    transformed_blocks,
    \(block) {
      sizes <- vapply(
        block,
        nrow,
        integer(1)
      )
      
      if (length(unique(sizes)) != 1L) {
        stop(
          "Locations have different numbers of observations ",
          "within the same block."
        )
      }
      
      sizes[[1]]
    },
    integer(1)
  )
  
  n_blocks <- length(transformed_blocks)
  
  #### Construct permutation groups ####
  
  if (is.null(permutation_groups)) {
    # Each observation is its own group. Prefixing by the original
    # block makes the identifiers unique after pooling.
    permutation_groups <- lapply(
      seq_len(n_blocks),
      \(block_index) {
        paste0(
          "block_",
          block_index,
          "_observation_",
          seq_len(block_sizes[[block_index]])
        )
      }
    )
    
    permutation_unit <- "observation"
  } else {
    if (
      !is.list(permutation_groups) ||
      length(permutation_groups) != n_blocks
    ) {
      stop(
        "`permutation_groups` must be NULL or a list of length two."
      )
    }
    
    group_lengths <- lengths(
      permutation_groups
    )
    
    if (!identical(
      as.integer(group_lengths),
      as.integer(block_sizes)
    )) {
      stop(
        "Each `permutation_groups` vector must have one value ",
        "per observation in its corresponding block. Expected ",
        paste(block_sizes, collapse = " and "),
        " values, but received ",
        paste(group_lengths, collapse = " and "),
        "."
      )
    }
    
    if (
      any(
        vapply(
          permutation_groups,
          anyNA,
          logical(1)
        )
      )
    ) {
      stop(
        "`permutation_groups` cannot contain missing values."
      )
    }
    
    # Prefix groups by the original block. This prevents accidental
    # collisions if the same label appears in both blocks.
    permutation_groups <- Map(
      \(groups, block_index) {
        paste0(
          "block_",
          block_index,
          "_group_",
          as.character(groups)
        )
      },
      permutation_groups,
      seq_len(n_blocks)
    )
    
    permutation_unit <- "group"
  }
  
  #### Convert both blocks to one long data frame ####
  
  data_laplace_combined <- bind_rows(
    lapply(
      seq_len(n_blocks),
      \(block_index) {
        block <- transformed_blocks[
          [block_index]
        ]
        
        bind_rows(
          Map(
            \(matrix_data, location) {
              matrix_data <- as.data.frame(
                matrix_data
              )
              
              if (
                !all(
                  c("X1", "X2") %in%
                  names(matrix_data)
                )
              ) {
                stop(
                  "Each transformed matrix must contain ",
                  "columns `X1` and `X2`."
                )
              }
              
              matrix_data |>
                transmute(
                  X1,
                  X2,
                  name = location,
                  original_block = block_index,
                  observation_index = row_number(),
                  permutation_group =
                    permutation_groups[
                      [block_index]
                    ]
                )
            },
            block,
            names(block)
          )
        )
      }
    )
  )
  
  #### Generate one common permutation ####
  
  group_table <- data_laplace_combined |>
    distinct(
      permutation_group,
      original_block
    )
  
  n_groups_block_1 <- group_table |>
    filter(original_block == 1L) |>
    nrow()
  
  n_groups_block_2 <- group_table |>
    filter(original_block == 2L) |>
    nrow()
  
  pooled_groups <- sample(
    group_table$permutation_group
  )
  
  block_assignment <- tibble(
    permutation_group = pooled_groups,
    block = c(
      rep(
        "1",
        n_groups_block_1
      ),
      rep(
        "2",
        n_groups_block_2
      )
    )
  )
  
  # Apply the same group assignment to every location.
  data_perm <- data_laplace_combined |>
    select(
      X1,
      X2,
      name,
      original_block,
      observation_index,
      permutation_group
    ) |>
    left_join(
      block_assignment,
      by = "permutation_group",
      relationship = "many-to-one"
    )
  
  if (anyNA(data_perm$block)) {
    stop(
      "At least one permutation group was not assigned to a block."
    )
  }
  
  #### Validate common assignment across locations ####
  
  assignment_check <- data_perm |>
    distinct(
      permutation_group,
      name,
      block
    ) |>
    count(
      permutation_group,
      name = "n_assignments"
    )
  
  if (any(
    assignment_check$n_assignments !=
    length(location_names)
  )) {
    stop(
      "A permutation group was not represented at every location."
    )
  }
  
  block_check <- data_perm |>
    distinct(
      permutation_group,
      block
    ) |>
    count(
      permutation_group,
      name = "n_blocks"
    )
  
  if (any(block_check$n_blocks != 1L)) {
    stop(
      "A permutation group was assigned to multiple blocks."
    )
  }
  
  #### Reconstruct cecl_marg objects ####
  
  permuted_blocks <- lapply(
    c("1", "2"),
    \(new_block) {
      block_data <- data_perm |>
        filter(block == new_block)
      
      location_data <- lapply(
        location_names,
        \(location) {
          block_data |>
            filter(name == location) |>
            # Ordering is not important for CE likelihood fitting, but
            # deterministic ordering makes debugging easier.
            arrange(
              original_block,
              observation_index
            ) |>
            select(
              X1,
              X2
            ) |>
            as.data.frame()
        }
      )
      
      names(location_data) <- location_names
      
      # Every location must have the same number of rows.
      location_sizes <- vapply(
        location_data,
        nrow,
        integer(1)
      )
      
      if (length(unique(location_sizes)) != 1L) {
        stop(
          "The common permutation produced unequal location ",
          "sample sizes within pseudo-block ", new_block, "."
        )
      }
      
      as_cecl_marg(location_data)
    }
  )
  
  names(permuted_blocks) <- c("1", "2")
  
  #### Fit CE models and calculate divergence matrices ####
  
  dist_perm <- calc_dist(
    marg = permuted_blocks,
    n_mc = n_mc,
    cond_prob = cond_prob,
    cond_val = cond_val,
    use_dth = use_dth,
    cond_var = cond_var,
    use_evgam = FALSE,
    aLow = aLow,
    nruns = nruns,
    ret_dep = FALSE,
    start = start,
    ...
  )
  
  # Preserve the existing convention that a failed CE fit returns NA.
  if (
    is.atomic(dist_perm) &&
    length(dist_perm) == 1L &&
    is.na(dist_perm)
  ) {
    return(NA)
  }
  
  attr(
    dist_perm,
    "permutation_unit"
  ) <- permutation_unit
  
  attr(
    dist_perm,
    "n_groups"
  ) <- c(
    block_1 = n_groups_block_1,
    block_2 = n_groups_block_2
  )
  
  dist_perm
}

# # Function to generate data for a single permutation test
# perm_test_prep <- \(
#   data_laplace_loc,
#   cond_prob = 0.9,
#   cond_val = NULL,
#   cond_var = NULL,
#   use_dth = FALSE,
#   use_evgam = FALSE,
#   data_loc = NULL,
#   aLow = -1,
#   nruns = 1,
#   start = list(c(a = 0.01, b = 0.01), c(a = 0.01, b = 0.01)),
#   # type = c("norm", "clustering"),
#   # norm_type = "F"
#   ...
# ) {
#   if (!is.null(cond_val)) {
#     cond_prob <- NULL
#   }
# 
#   # combine data into single data.frame with block indicator
#   data_laplace_combined <- bind_rows(lapply(seq_along(data_laplace_loc), \(i) {
#     bind_rows(lapply(data_laplace_loc[[i]]$transformed, as.data.frame), .id = "name") %>%
#       mutate(block = i)
#   })) |>
#     select(X1, X2, name, block)
# 
#   # permute rows within each location
#   data_perm <- data_laplace_combined %>%
#     group_by(name) %>%
#     group_modify(\(df, key) {
#       # original block sizes for this location
#       n1 <- sum(df$block == levels(factor(df$block))[1])
#       n2 <- sum(df$block == levels(factor(df$block))[2])
# 
#       # pool rows and randomly split
#       idx <- sample.int(nrow(df))
#       idx1 <- idx[seq_len(n1)]
#       idx2 <- idx[-seq_len(n1)]
# 
#       df |>
#         mutate(
#           row = row_number(),
#           block = case_when(
#             row %in% idx1 ~ levels(factor(df$block))[1],
#             row %in% idx2 ~ levels(factor(df$block))[2]
#           )
#         )
#     }) %>%
#     ungroup() |>
#     select(-row)
# 
#   locs <- unique(data_perm$name)
# 
#   # convert to cecl_marg objects (UNCHANGED)
#   data_laplace_perm <- data_perm %>%
#     mutate(block = factor(block, levels = unique(block))) %>%
#     group_split(block) %>%
#     lapply(\(x) {
#       ret <- lapply(locs, \(loc) {
#         subset(x, name == loc) %>%
#           select(X1, X2)
#       })
#       names(ret) <- locs
#       ret
#     })
#   names(data_laplace_perm) <- levels(data_perm$block)
# 
#   marg_perm <- lapply(data_laplace_perm, as_cecl_marg)
# 
#   # fit CE
#   # TODO Change!!
#   # dep_perm <- lapply(
#   #   marg_perm,
#   #   cecl_dep,
#   #   cond_prob = cond_prob,
#   #   cond_val = cond_val,
#   #   fit_no_keef = TRUE
#   # )
#   # fit CE models
#   if (use_evgam == FALSE) {
#     # dep_perm <- lapply(
#     #   marg_perm,
#     #   cecl_dep,
#     #   cond_prob   = cond_prob,
#     #   cond_val    = cond_val,
#     #   fit_no_keef = TRUE,
#     #   aLow        = aLow,
#     #   nruns       = nruns
#     # )
#     dep_perm <- lapply(seq_along(marg_perm), \(i) {
#       cecl_dep(
#         obj         = marg_perm[[i]],
#         cond_prob   = cond_prob,
#         cond_val    = cond_val,
#         cond_var    = cond_var,
#         aLow        = aLow,
#         nruns       = nruns,
#         fit_no_keef = TRUE,
#         start       = start[[i]]
#       )
#     })
#     # skip bad permutations
#     if (any(sapply(dep_perm, \(x) any(is.na(unlist(x$dependence)))))) {
#       return(NA)
#     }
#   } else if (use_evgam) {
#     if (is.null(data_loc)) {
#       stop("data_loc must be provided when use_evgam = TRUE")
#     }
# 
#     # first, put marginal data into form for evgam fitting
#     marg_join_lst <- lapply(seq_along(marg_perm), \(k) {
#       x <- marg_perm[[k]]
#       bind_rows(lapply(seq_along(x$transformed), \(j) {
#         data.frame(x$transformed[[j]]) |>
#           mutate(name = paste0("loc_", j))
#       })) |>
#         left_join(
#           data_loc |>
#             filter(block == k) |>
#             distinct(name, x_loc, y_loc),
#           by = "name"
#         ) |>
#         rename(x = x_loc, y = y_loc) |>
#         select(-name)
#     })
# 
#     dep_perm <- lapply(marg_join_lst, \(x) {
#       x1_cond <- fit_evgam(
#         df = x,
#         dep_val = cond_val,
#         var = "X2",
#         cond_var = "X1"
#       )
# 
#       x2_cond <- fit_evgam(
#         df = x,
#         dep_val = cond_val,
#         var = "X1",
#         cond_var = "X2"
#       )
# 
#       # join and change to `cecl_dep` format
#       coef_evgam <- bind_rows(x1_cond$predictions, x2_cond$predictions) |>
#         arrange(name, var, cond_var)
#       rownames(coef_evgam) <- NULL
#       as_cecl_dep(coef_evgam)
#     })
#   } else {
#     stop("Invalid value for use_evgam")
#   }
# 
# 
#   thresh_max <- NULL
#   if (use_dth == TRUE) {
#     thresh <- lapply(dep_perm, \(x) {
#       lapply(x$dependence, CeCl:::pull_thresh_trans)
#     })
#     # take max
#     thresh_max <- max(unlist(thresh))
#   }
# 
#   # compute distances
# 
#   dist_perm <- lapply(seq_along(dep_perm), \(i) {
#     cecl_dist(
#       dep_perm[[i]],
#       marg_perm[[i]],
#       dth = thresh_max,
#       ...
#     )
#   })
# }
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

# # Function to perform permutation test for a given block position i
# # TODO Optimise better to work in parallel, not great at the moment
# perm_test_fun <- \(
#   data,
#   grid_vals,
#   n_perm = 100,
#   n_per_block = NULL,
#   n_years_per_block = NULL,
#   n_vars = 2,
#   laplace_trans = TRUE,
#   use_start = FALSE,
#   ret_dep = FALSE,
#   use_dth = FALSE,
#   verbose = FALSE,
#   ...
# ) {
#   # initial parameter values for CE
#   start <- list(
#     c(a = 0.01, b = 0.01),
#     c(a = 0.01, b = 0.01)
#   )
#   # if use_start is TRUE, get starting values from fitting CE to original data
#   if (use_start) {
#     data_start <- data
#     if (!"block" %in% names(data)) {
#       data_start <- data_start |>
#         # fit CE with half of observations in each block
#         mutate(block = ifelse(row_number() <= n() / 2, "1", "2"))
#     }
#     data_laplace_start <- trans_fun(data_start, n_vars, laplace_trans)
#
#     dots <- list(...)
#     # if (!"cond_val" %in% names(dots)) {
#     #   cond_val_start <- qlaplace(0.8)
#     # } else {
#     #   cond_val_start <- dots$cond_val
#     # }
#     if (!"cond_val" %in% names(dots)) {
#       cond_val_start <- NULL
#     } else {
#       cond_val_start <- dots$cond_val
#     }
#     if (!"cond_prob" %in% names(dots)) {
#       cond_prob_start <- NULL
#     } else {
#       cond_prob_start <- dots$cond_prob
#     }
#
#     if (!"cond_var" %in% names(dots)) {
#       cond_var_start <- NULL
#     } else {
#       cond_var_start <- dots$cond_var
#     }
#     if (!"nruns" %in% names(dots)) {
#       nruns_start <- 1
#     } else {
#       nruns_start <- dots$nruns
#     }
#     if (!"aLow" %in% names(dots)) {
#       aLow_start <- -1
#     } else {
#       aLow_start <- dots$aLow
#     }
#
#     dep_start <- lapply(seq_along(data_laplace_start), \(k) {
#       cecl_dep(
#         data_laplace_start[[k]],
#         cond_val = cond_val_start,
#         cond_prob = cond_prob_start,
#         fit_no_keef = TRUE,
#         nruns = nruns_start,
#         cond_var = cond_var_start,
#         aLow = aLow_start,
#         start = start[[k]]
#       )
#     })
#
#     # start values
#     start <- lapply(dep_start, coef)
#   }
#
#   # keep mclapply for where it's needed most
#   i_fun <- j_fun <- lapply
#   if (length(grid_vals) > n_perm) {
#     i_fun <- mclapply
#   } else {
#     j_fun <- mclapply
#   }
#
#   # fun inside apply
#   ret_fun <- \(i) {
#      if (verbose) {
#       system(sprintf(
#         'echo "\n%s\n"',
#         paste("repitition", which(grid_vals == i))
#       ))
#     }
#
#     # Add block based on i
#     data_block <- data %>%
#       group_by(name) |>
#       mutate(block = ifelse(row_number() <= i, "1", "2")) |>
#       ungroup()
#
#     # if desired, take n_per_block obs for each location and block
#     if (!is.null(n_per_block)) {
#       # give error if we want more observations than we can have at point i
#       stopifnot(i < n_per_block)
#
#       data_block <- data_block |>
#         group_split(block, .keep = TRUE) |>
#         lapply(\(x) {
#           y <- x %>%
#             group_by(name)
#           # for first block, take last n_per_block obs
#           if (as.numeric(y$block[1]) == 1) {
#             y <- y |>
#               slice_tail(n = n_per_block)
#             # for second block, take first n_per_block obs
#           } else if (as.numeric(y$block[1]) == 2) {
#             y <- y |>
#               slice_head(n = n_per_block)
#           }
#           return(ungroup(y))
#         }) |>
#         bind_rows()
#     }
#
#     # convert to Laplace and store marginals
#     data_laplace_block <- trans_fun(data_block, n_vars, laplace_trans)
#
#     # compute distances for original data (before permutation)
#     dist <- suppressMessages(calc_dist(
#       data_laplace_block,
#       start = start, ret_dep = ret_dep, use_dth = use_dth, ...
#     ))
#     if (ret_dep == TRUE) {
#       dep <- dist$dep
#       dist <- dist$dist
#     }
#
#     # compute different norms
#     norm_orig_frob <- compare_blocks(dist, type = "norm", norm_type = "F")
#     norm_orig_inf <- compare_blocks(dist, type = "norm", norm_type = "M")
#     norm_orig_spec <- compare_blocks(dist, type = "norm", norm_type = "2")
#     norm_orig_ln_spd <- compare_blocks(dist, type = "norm", norm_type = "log SPD")
#
#     # compute permutation distances for original data
#     norm_vals <- j_fun(seq_len(n_perm), \(j) {
#       # print progress of permutations
#       if (verbose) {
#         system(sprintf("echo Permutation %s", j))
#       }
#       dist <- perm_test_prep(
#         data_laplace_block,
#         start = start,
#         use_dth = use_dth,
#         ...
#       )
#       list(
#         "frob"   = compare_blocks(dist, type = "norm", norm_type = "F"),
#         "inf"    = compare_blocks(dist, type = "norm", norm_type = "M"),
#         "spec"   = compare_blocks(dist, type = "norm", norm_type = "2"),
#         "ln_spd" = compare_blocks(dist, type = "norm", norm_type = "log SPD")
#       )
#     })
#
#     # extract frobenius and infinity norms
#     perm_norms_frob <- sapply(norm_vals, \(x) x$frob)
#     perm_norms_inf <- sapply(norm_vals, \(x) x$inf)
#     perm_norms_spec <- sapply(norm_vals, \(x) x$spec)
#     perm_norms_ln_spd <- sapply(norm_vals, \(x) x$ln_spd)
#
#     # compute p-values
#     # TODO Are these p-values for one-sided or two-sided tests?
#     (p_value_frob <- mean(perm_norms_frob >= norm_orig_frob))
#     (p_value_inf <- mean(perm_norms_inf >= norm_orig_inf))
#     (p_value_spec <- mean(perm_norms_spec >= norm_orig_spec))
#     (p_value_ln_spd <- mean(perm_norms_ln_spd >= norm_orig_ln_spd))
#
#     # print how many iterations completed (as a percentage)
#     if (verbose) {
#       system(sprintf(
#         'echo "\n%s\n"',
#         paste0(round(which(grid_vals == i) / length(grid_vals), 3) * 100, "% completed", collapse = "")
#       ))
#     }
#
#     ret_list <- list(
#       p_value_frob = p_value_frob,
#       p_value_inf = p_value_inf,
#       p_value_spec = p_value_spec,
#       p_value_ln_spd = p_value_ln_spd,
#       norm_orig_frob = norm_orig_frob,
#       norm_orig_inf = norm_orig_inf,
#       norm_orig_spec = norm_orig_spec,
#       norm_orig_ln_spd = norm_orig_ln_spd,
#       perm_norms_frob = perm_norms_frob,
#       perm_norms_inf = perm_norms_inf,
#       perm_norms_spec = perm_norms_spec,
#       perm_norms_ln_spd = perm_norms_ln_spd
#     )
#     if (ret_dep == TRUE) {
#       dep_vals_df <- lapply(dep, coef) |>
#         bind_rows(.id = "block")
#       ret_list <- c(ret_list, list(dep_vals = dep_vals_df))
#     }
#     return(ret_list)
#   }
#   return(i_fun(grid_vals, ret_fun))
# }

# Function to perform permutation tests at candidate row or seasonal-year
# boundaries.
#
# Modes:
#   1. n_per_block supplied:
#      grid_vals are row positions and fixed-length observation windows are used.
#
#   2. n_years_per_block supplied:
#      grid_vals are seasonal years and fixed-length year windows are used.
#
#   3. Both are NULL:
#      grid_vals are row positions and all observations on either side are used.
#
# Do not supply both n_per_block and n_years_per_block.
perm_test_fun <- \(
  data,
  grid_vals,
  n_perm = 100,
  n_per_block = NULL,
  n_years_per_block = NULL,
  n_vars = 2,
  laplace_trans = TRUE,
  use_start = FALSE,
  ret_dep = FALSE,
  use_dth = FALSE,
  verbose = FALSE,
  ...
) {
  #### Validate arguments ####

  if (!is.null(n_per_block) && !is.null(n_years_per_block)) {
    stop(
      "Supply only one of `n_per_block` and ",
      "`n_years_per_block`."
    )
  }

  if (length(grid_vals) == 0L) {
    stop("`grid_vals` must contain at least one candidate boundary.")
  }

  if (
    length(n_perm) != 1L ||
      is.na(n_perm) ||
      n_perm < 1L
  ) {
    stop("`n_perm` must be a positive integer.")
  }

  row_window_mode <- !is.null(n_per_block)
  year_window_mode <- !is.null(n_years_per_block)
  cumulative_mode <- is.null(n_per_block) &&
    is.null(n_years_per_block)

  if (row_window_mode) {
    if (
      length(n_per_block) != 1L ||
        is.na(n_per_block) ||
        n_per_block < 1L
    ) {
      stop("`n_per_block` must be a positive integer.")
    }

    n_per_block <- as.integer(n_per_block)
  }

  if (year_window_mode) {
    required_columns <- c("name", "date", "season_year")
    missing_columns <- setdiff(required_columns, names(data))

    if (length(missing_columns) > 0L) {
      stop(
        "Year-based windows require the following missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    if (
      length(n_years_per_block) != 1L ||
        is.na(n_years_per_block) ||
        n_years_per_block < 1L
    ) {
      stop("`n_years_per_block` must be a positive integer.")
    }

    n_years_per_block <- as.integer(n_years_per_block)

    all_years <- sort(unique(data$season_year))

    if (length(all_years) < 2L * n_years_per_block) {
      stop(
        "There are only ", length(all_years),
        " seasonal years, but at least ",
        2L * n_years_per_block,
        " are needed."
      )
    }

    invalid_grid_vals <- setdiff(grid_vals, all_years)

    if (length(invalid_grid_vals) > 0L) {
      stop(
        "These year-based `grid_vals` are not present in ",
        "`season_year`: ",
        paste(invalid_grid_vals, collapse = ", ")
      )
    }
  }

  #### Initial parameter values ####

  start <- replicate(
    n_vars,
    c(a = 0.01, b = 0.01),
    simplify = FALSE
  )

  if (use_start) {
    data_start <- data

    if (!"block" %in% names(data_start)) {
      if (year_window_mode) {
        # Use the middle candidate boundary to construct representative
        # starting-value blocks containing whole seasonal years.
        start_grid_val <- grid_vals[
          ceiling(length(grid_vals) / 2)
        ]

        start_year_position <- match(
          start_grid_val,
          all_years
        )

        start_left_years <- all_years[
          (start_year_position - n_years_per_block + 1L):
          start_year_position
        ]

        start_right_years <- all_years[
          (start_year_position + 1L):
          (start_year_position + n_years_per_block)
        ]

        data_start <- data_start |>
          filter(
            season_year %in% c(
              start_left_years,
              start_right_years
            )
          ) |>
          mutate(
            block = if_else(
              season_year %in% start_left_years,
              "1",
              "2"
            )
          ) |>
          arrange(name, date)
      } else {
        # Original observation-based behaviour
        if ("date" %in% names(data_start)) {
          data_start <- data_start |>
            arrange(name, date)
        }

        data_start <- data_start |>
          group_by(name) |>
          mutate(
            block = if_else(
              row_number() <= floor(n() / 2),
              "1",
              "2"
            )
          ) |>
          ungroup()
      }
    }

    data_laplace_start <- trans_fun(
      data_start,
      n_vars,
      laplace_trans
    )

    dots <- list(...)

    cond_val_start <- if ("cond_val" %in% names(dots)) {
      dots$cond_val
    } else {
      NULL
    }

    cond_prob_start <- if ("cond_prob" %in% names(dots)) {
      dots$cond_prob
    } else {
      NULL
    }

    cond_var_start <- if ("cond_var" %in% names(dots)) {
      dots$cond_var
    } else {
      NULL
    }

    nruns_start <- if ("nruns" %in% names(dots)) {
      dots$nruns
    } else {
      1
    }

    aLow_start <- if ("aLow" %in% names(dots)) {
      dots$aLow
    } else {
      -1
    }

    dep_start <- lapply(
      seq_along(data_laplace_start),
      \(k) {
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
      }
    )

    start <- lapply(dep_start, coef)
  }

  #### Choose level of parallelisation ####

  i_fun <- lapply
  j_fun <- lapply

  if (length(grid_vals) > n_perm) {
    i_fun <- mclapply
  } else {
    j_fun <- mclapply
  }

  #### Test one candidate boundary ####

  ret_fun <- \(i) {
    # browser()
    grid_position <- match(i, grid_vals)

    if (verbose) {
      message(
        "\nCandidate ",
        grid_position,
        " of ",
        length(grid_vals),
        ": ",
        i
      )
    }

    #### Construct observed blocks ####

    if (year_window_mode) {
      # In year mode, i means that the candidate change occurs
      # between season_year i and the following available season_year.

      year_position <- match(i, all_years)

      if (year_position < n_years_per_block) {
        stop(
          "Grid value ", i,
          " does not have ", n_years_per_block,
          " seasonal years available on the left."
        )
      }

      if (
        length(all_years) - year_position <
          n_years_per_block
      ) {
        stop(
          "Grid value ", i,
          " does not have ", n_years_per_block,
          " seasonal years available on the right."
        )
      }

      left_years <- all_years[
        (year_position - n_years_per_block + 1L):
        year_position
      ]

      right_years <- all_years[
        (year_position + 1L):
        (year_position + n_years_per_block)
      ]

      data_block <- data |>
        filter(
          season_year %in% c(left_years, right_years)
        ) |>
        mutate(
          block = if_else(
            season_year %in% left_years,
            "1",
            "2"
          )
        ) |>
        arrange(name, date)

      # Confirm that no seasonal year has been split between blocks.
      year_check <- data_block |>
        distinct(season_year, block) |>
        count(season_year, name = "n_blocks")

      if (any(year_check$n_blocks != 1L)) {
        stop(
          "At least one seasonal year was assigned to ",
          "more than one block."
        )
      }

      # Confirm that every station has the same date coverage within
      # each block.
      station_block_sizes <- data_block |>
        distinct(name, block, date) |>
        count(name, block, name = "n_dates")

      coverage_check <- station_block_sizes |>
        group_by(block) |>
        summarise(
          n_distinct_sizes = n_distinct(n_dates),
          .groups = "drop"
        )

      if (any(coverage_check$n_distinct_sizes != 1L)) {
        stop(
          "Stations do not have equal date coverage within ",
          "at least one year-based block."
        )
      }

      n_left_dates <- data_block |>
        filter(block == "1") |>
        summarise(n = n_distinct(date)) |>
        pull(n)

      n_right_dates <- data_block |>
        filter(block == "2") |>
        summarise(n = n_distinct(date)) |>
        pull(n)

      block_info <- list(
        mode = "season_year",
        split_after = i,
        split_before = all_years[year_position + 1L],
        left_years = left_years,
        right_years = right_years,
        n_left = n_left_dates,
        n_right = n_right_dates
      )
    } else {
      # In row mode, i is a row position within each station.

      data_ordered <- data

      if ("date" %in% names(data_ordered)) {
        data_ordered <- data_ordered |>
          arrange(name, date)
      }

      n_by_station <- data_ordered |>
        count(name, name = "n_observations")

      if (i < 1L) {
        stop("All row-based `grid_vals` must be positive.")
      }

      if (any(i >= n_by_station$n_observations)) {
        stop(
          "Grid value ", i,
          " does not leave observations on both sides ",
          "for every station."
        )
      }

      data_block <- data_ordered |>
        group_by(name) |>
        mutate(
          block = if_else(
            row_number() <= i,
            "1",
            "2"
          )
        ) |>
        ungroup()

      if (row_window_mode) {
        if (i < n_per_block) {
          stop(
            "Grid value ", i,
            " leaves fewer than ", n_per_block,
            " observations in the left block."
          )
        }

        if (
          any(
            n_by_station$n_observations - i <
              n_per_block
          )
        ) {
          stop(
            "Grid value ", i,
            " leaves fewer than ", n_per_block,
            " observations in the right block for ",
            "at least one station."
          )
        }

        left_block <- data_block |>
          filter(block == "1") |>
          group_by(name) |>
          slice_tail(n = n_per_block) |>
          ungroup()

        right_block <- data_block |>
          filter(block == "2") |>
          group_by(name) |>
          slice_head(n = n_per_block) |>
          ungroup()

        data_block <- bind_rows(
          left_block,
          right_block
        )

        block_info <- list(
          mode = "fixed_rows",
          split_after = i,
          split_before = i + 1L,
          left_years = NULL,
          right_years = NULL,
          n_left = n_per_block,
          n_right = n_per_block
        )
      } else if (cumulative_mode) {
        block_info <- list(
          mode = "cumulative_rows",
          split_after = i,
          split_before = i + 1L,
          left_years = NULL,
          right_years = NULL,
          n_left = i,
          n_right = min(
            n_by_station$n_observations - i
          )
        )
      }
    }

    #### Transform and fit observed blocks ####
    
    data_laplace_block <- trans_fun(
      # data_block,
      select(data_block, X1, X2, name, block),
      n_vars,
      laplace_trans
    )
    
    dist <- suppressMessages(
      calc_dist(
        data_laplace_block,
        start = start,
        ret_dep = ret_dep,
        use_dth = use_dth,
        ...
      )
    )
    
    if (ret_dep) {
      dep <- dist$dep
      dist <- dist$dist
    }

    #### Observed statistics ####

    norm_orig_frob <- compare_blocks(
      dist,
      type = "norm",
      norm_type = "F"
    )

    norm_orig_inf <- compare_blocks(
      dist,
      type = "norm",
      norm_type = "M"
    )

    norm_orig_spec <- compare_blocks(
      dist,
      type = "norm",
      norm_type = "2"
    )

    # norm_orig_ln_spd <- compare_blocks(
    #   dist,
    #   type = "norm",
    #   norm_type = "log SPD"
    # )

    #### Permutation statistics ####

    norm_vals <- j_fun(
      seq_len(n_perm),
      \(j) {
        if (verbose) {
          # message(
          #   "Candidate ", grid_position,
          #   "/", length(grid_vals),
          #   "; permutation ", j,
          #   "/", n_perm
          # )
          # print progress of permutations
          system(sprintf("echo Permutation %s", j))
        }

        # dist_perm <- perm_test_prep(
        #   data_laplace_block,
        #   start = start,
        #   use_dth = use_dth,
        #   ...
        # )
        dist_perm <- perm_test_prep(
          data_laplace_block,
          start = start,
          use_dth = use_dth,
          cond_val = dep_val,
          laplace_sample = laplace_sample,
          ...
        )

        list(
          frob = compare_blocks(
            dist_perm,
            type = "norm",
            norm_type = "F"
          ),
          inf = compare_blocks(
            dist_perm,
            type = "norm",
            norm_type = "M"
          ),
          spec = compare_blocks(
            dist_perm,
            type = "norm",
            norm_type = "2"
          )# ,
          # ln_spd = compare_blocks(
          #   dist_perm,
          #   type = "norm",
          #   norm_type = "log SPD"
          # )
        )
      }
    )

    perm_norms_frob <- vapply(
      norm_vals,
      \(x) x$frob,
      numeric(1)
    )

    perm_norms_inf <- vapply(
      norm_vals,
      \(x) x$inf,
      numeric(1)
    )

    perm_norms_spec <- vapply(
      norm_vals,
      \(x) x$spec,
      numeric(1)
    )

    # perm_norms_ln_spd <- vapply(
    #   norm_vals,
    #   \(x) x$ln_spd,
    #   numeric(1)
    # )

    #### P-values ####

    # p_value_frob <- mean(
    #   perm_norms_frob >= norm_orig_frob
    # )
    p_value_frob <- (
      1 + sum(perm_norms_frob >= norm_orig_frob)
    ) / (length(perm_norms_frob) + 1)

    # p_value_inf <- mean(
    #   perm_norms_inf >= norm_orig_inf
    # )
    p_value_inf <- (
      1 + sum(perm_norms_inf >= norm_orig_inf)
    ) / (length(perm_norms_inf) + 1)

    # p_value_spec <- mean(
    #   perm_norms_spec >= norm_orig_spec
    # )
    p_value_spec <- (
      1 + sum(perm_norms_spec >= norm_orig_spec)
    ) / (length(perm_norms_spec) + 1)

    # p_value_ln_spd <- mean(
    #   perm_norms_ln_spd >= norm_orig_ln_spd
    # )

    if (verbose) {
      message(
        round(
          grid_position / length(grid_vals) * 100,
          1
        ),
        "% of candidate boundaries completed"
      )
    }

    #### Return results ####

    ret_list <- list(
      block_info = block_info,
      p_value_frob = p_value_frob,
      p_value_inf = p_value_inf,
      p_value_spec = p_value_spec,
      # p_value_ln_spd = p_value_ln_spd,
      norm_orig_frob = norm_orig_frob,
      norm_orig_inf = norm_orig_inf,
      norm_orig_spec = norm_orig_spec,
      # norm_orig_ln_spd = norm_orig_ln_spd,
      perm_norms_frob = perm_norms_frob,
      perm_norms_inf = perm_norms_inf,
      perm_norms_spec = perm_norms_spec # ,
      # perm_norms_ln_spd = perm_norms_ln_spd
    )

    if (ret_dep) {
      dep_vals_df <- lapply(dep, coef) |>
        bind_rows(.id = "block")

      ret_list$dep_vals <- dep_vals_df
    }

    ret_list
  }

  i_fun(grid_vals, ret_fun)
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

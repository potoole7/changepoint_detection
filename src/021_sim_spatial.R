#### Spatially varying simulation ####

# TODO Also try for k = 2 for second block? Might have very good performance

# TODO Perform permutation test (and sliding) using evgam to fit CE
# (on spatial simulated data)

#### libs ####

# library(CeCl)
devtools::load_all("../CeCl")
# library(sp)
library(copula)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tidyr)
library(evgam)

#### Metadata ####

seed <- 123
n_vars <- 2 # number of variables
# n_locs <- 60 # number of locations
n_locs <- 120 # number of locations
n_clust <- 3 # number of clusters
# n <- 1000 * (12 / 3) # total number of observations per group
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per location times # locs per group
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # correlation parameters for each group


#### Functions ####

gen_t_loc <- function(rho, n, n_vars = 2, df_t = 3) {
  cop_t <- tCopula(param = rho, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}

#### Method ####

# Function to generate spatial data with varying correlation
gen_spatial_data <- \(
  cor_t_vec,
  n_clust = 3,
  n_locs = 60,
  n_per_loc = 1000,
  sd_loc = 1, # spatial spread around cluster centres
  ell = 0.5, # length scale for spatial smoothing (bigger == smoother)
  eps = 1e-6 # to keep rho inside (-1, 1)
) {
  # generate points for cluster centres
  clust_centres <- data.frame(
    cor = cor_t_vec,
    x = runif(n_clust, 0, 5),
    y = runif(n_clust, 0, 5)
  )

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

  # sd_loc <- 1.0 # spatial spread around each centre (tune this)

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
    x <- gen_t_loc(rho_loc[ell], n_per_loc, n_vars = n_vars, df_t = df_t)
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

# generate data with spatially varying correlation
set.seed(seed)
data_spat <- gen_spatial_data(
  cor_t_vec = cor_t_vec,
  n_clust = n_clust,
  n_locs = n_locs,
  n_per_loc = n_per_loc,
  # sd_loc = 1,
  # sd_loc = 0.5
  sd_loc = 0.8
)

# investigate effect of length scale
sapply(c(0.5, 1, 2, 3, 4, 5), \(x) {
  set.seed(seed)
  gen_spatial_data(
    cor_t_vec = cor_t_vec,
    n_clust = n_clust,
    n_locs = n_locs,
    n_per_loc = n_per_loc,
    ell = x
  )
  return(0)
})
# Use 0.5, as it gives more spatial variability

# plot estimated vs true correlation
data_spat |>
  group_by(name) |>
  summarise(
    cor = cor(X1, X2),
    rho = first(rho),
    x = first(x_loc),
    y = first(y_loc),
    cluster = first(cluster)
  ) |>
  pivot_longer(cols = c(cor, rho), names_to = "type", values_to = "value") |>
  ggplot(aes(x = x, y = y, colour = value, shape = factor(cluster))) +
  geom_point(size = 4) +
  facet_wrap(~type) +
  scale_colour_viridis_c() +
  labs(shape = "Cluster ID") +
  theme_minimal()
# Very similar!

data_spat |>
  group_by(cluster) |>
  summarise(
    cor = cor(X1, X2),
    rho = mean(rho),
    .groups = "drop"
  )

# If using blocks, ...
# TODO Fill in here
diff <- 0.4
cor_t_vec2 <- cor_t_vec
cor_t_vec2[[1]] <- cor_t_vec2[[1]] + diff

set.seed(seed)
# debugonce(gen_spatial_data)
data_spat2 <- gen_spatial_data(
  cor_t_vec = cor_t_vec2,
  n_clust = n_clust,
  n_locs = n_locs,
  n_per_loc = n_per_loc
)

data_spat_block <- bind_rows(
  mutate(data_spat, block = 1),
  mutate(data_spat2, block = 2)
)

data_spat_block |>
  mutate(name = stringr::str_remove(name, "loc_")) |>
  group_by(name, block) |>
  summarise(
    cor = cor(X1, X2),
    rho = first(rho),
    x = first(x_loc),
    y = first(y_loc),
    cluster = first(cluster),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(cor, rho), names_to = "type", values_to = "value") |>
  ggplot(aes(x = x, y = y, colour = value, shape = factor(cluster))) +
  geom_point(size = 4) +
  # also add label of location names
  ggrepel::geom_text_repel(
    aes(label = name),
    size = 3,
    max.overlaps = 10
  ) +
  # facet_wrap(~type) +
  facet_grid(type ~ block) +
  scale_colour_viridis_c() +
  labs(shape = "Cluster ID") +
  theme_minimal()

data_spat_block |>
  group_by(cluster, block) |>
  summarise(
    cor = cor(X1, X2),
    rho = mean(rho),
    .groups = "drop"
  )

#### Fit CE (vanilla + evgam), calculate divergence and compare ####

data_spat_block_lst <- data_spat_block |>
  select(-c(x_loc, y_loc, rho, cluster)) |>
  group_split(block, .keep = FALSE)

# Convert to Laplace

marg_lst <- lapply(data_spat_block_lst, \(x) {
  cecl_marg(
    x,
    thresh_method = "value",
    thresh_args = 0
  )
})

# Fit CE (vanilla)

# condition value for dependence calculation
# (dep_val <- qgpd(0.8, u = 0, xi = -0.05, sigma = 1)) #
dep_val <- qlaplace(p = 0.8)
dep_lst <- lapply(
  marg_lst,
  cecl_dep,
  # cond_prob   = cond_prob,
  cond_val    = dep_val,
  fit_no_keef = TRUE
)

# Calculate divergence
dist <- lapply(seq_along(dep_lst), \(i) {
  cecl_dist(
    dep_lst[[i]],
    marg_lst[[i]],
    seed = seed
  )
})

ggplot(dist[[1]], which = "image") # can see spatial structure
ggplot(dist[[2]], which = "image") # can see lcoation 9 is different, as it's closer to cluster 2

# cluster
clust_lst <- lapply(dist, \(x) {
  cecl_clust(
    x,
    k = 3,
    cluster_mem = sort(rep(unique(data_spat_block$cluster), n_locs / n_clust))
  )
})

# print
clust_lst[[1]]$adj_rand
clust_lst[[2]]$adj_rand
# for n = 60, block 1 = 0.856, block 2 = 0.468
# for n = 120, block 1 = 0.903, block 2 = 0.491

# also try for k = 2
clust_lst_k2 <- lapply(dist, \(x) {
  cecl_clust(
    x,
    k = 2,
    cluster_mem = sort(rep(c(1, 1, 2), times = n_locs / n_clust))
  )
})

clust_lst_k2[[1]]$adj_rand # 1
clust_lst_k2[[2]]$adj_rand # 1


# TODO Now fit with evgam!
# first, put marginal data into form for evgam fitting
marg_join_lst <- lapply(seq_along(marg_lst), \(i) {
  x <- marg_lst[[i]]
  bind_rows(lapply(seq_along(x$transformed), \(j) {
    data.frame(x$transformed[[j]]) |>
      mutate(name = paste0("loc_", j))
  })) |>
    left_join(
      data_spat_block |>
        filter(block == i) |>
        distinct(name, x_loc, y_loc),
      by = "name"
    ) |>
    rename(x = x_loc, y = y_loc) |>
    select(-name)
})

# Function to fit CE model using evgam for a given conditioning variable
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

# debugonce(fit_evgam)
# x1_cond <- fit_evgam(
#   df = marg_join_lst[[1]],
#   dep_val = qlaplace(p = 0.8),
#   var = "X2",
#   cond_var = "X1"
# )

# fit evgam models for each block
x <- marg_join_lst[[1]]
evgam_dep_lst <- lapply(marg_join_lst, \(x) {
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

# Calculate divergence
evgam_dist <- lapply(seq_along(evgam_dep_lst), \(i) {
  cecl_dist(
    dep_lst[[i]],
    marg_lst[[i]],
    seed = seed
  )
})

ggplot(evgam_dist[[1]], which = "image") # can see spatial structure
ggplot(evgam_dist[[2]], which = "image") # can see lcoation 9 is different, as it's closer to cluster 2

# cluster
evgam_clust_lst <- lapply(evgam_dist, \(x) {
  cecl_clust(
    x,
    k = 3,
    cluster_mem = sort(rep(unique(data_spat_block$cluster), n_locs / n_clust))
  )
})

evgam_clust_lst[[1]]$adj_rand
evgam_clust_lst[[2]]$adj_rand
# for n = 60, block 1 = 0.856, block 2 = 0.468
# for n = 120, block 1 = 0.903, block 2 = 0.491, so also the same ...

# also try for k = 2
evgam_clust_lst_k2 <- lapply(evgam_dist, \(x) {
  cecl_clust(
    x,
    k = 2,
    cluster_mem = sort(rep(c(1, 1, 2), times = n_locs / n_clust))
  )
})

evgam_clust_lst_k2[[1]]$adj_rand
evgam_clust_lst_k2[[2]]$adj_rand
# for n = 120, both 1

# TODO Plot alpha, beta estimates for vanilla and evgam fits and compare
# Join coefs to locations, plot on map for each block and method
plot_df <- bind_rows(lapply(seq_along(dep_lst), \(i) {
  coef(dep_lst[[i]]) |>
    mutate(block = i, method = "vanilla") |>
    left_join(
      data_spat_block |>
        filter(block == i) |>
        distinct(name, x_loc, y_loc),
      by = "name"
    )
})) |>
  bind_rows(
    bind_rows(lapply(seq_along(evgam_dep_lst), \(i) {
      coef(evgam_dep_lst[[i]]) |>
        mutate(block = i, method = "evgam") |>
        left_join(
          data_spat_block |>
            filter(block == i) |>
            distinct(name, x_loc, y_loc),
          by = "name"
        )
    }))
  )

plot_df |>
  select(-c("m", "s")) |>
  pivot_longer(cols = c(a, b), names_to = "param", values_to = "value") |>
  mutate(
    block = ifelse(block == 1, "Block 1", "Block 2"),
    method = factor(method, levels = c("vanilla", "evgam"))
  ) |>
  ggplot(aes(x = x_loc, y = y_loc, colour = value)) +
  geom_point(size = 3) +
  facet_grid(param ~ block + method) +
  scale_colour_viridis_c() +
  labs(x = "Longitude", y = "Latitude", colour = "Value") +
  cecl_theme(nejm_pal = FALSE)
# nice, (particularly) alpha and beta estimates do seem a lot more stable/smooth
# using evgam, as desired

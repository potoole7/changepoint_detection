#### Spatially-varying simulation ####

# TODO Allow specification of t-copula correlation for certain "centres"

#### libs ####

library(CeCl)
library(sp)
library(copula)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)

#### Metadata ####

seed <- 123
n_vars <- 2 # number of variables
n_locs <- 60 # number of locations
n_clust <- 3 # number of clusters
# n <- 1000 * (12 / 3) # total number of observations per group
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per location times # locs per group
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # correlation parameters for each group


#### Functions ####

gen_t_loc <- function(rho, n, n_vars = 3, df_t = 3) {
  cop_t <- tCopula(param = rho, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}

#### Generate Data: Option 1, Spatial Kernel ####

# Uses exponential kernel to smooth across space

set.seed(seed)
# generate random spatial coordinates for locations
coords <- cbind(runif(n_locs), runif(n_locs)) # or your real lon/lat (projected)

D <- as.matrix(dist(coords)) # Euclidean distances

# choose centres as points that are maximally distant
centre_loc <- c(1)
K <- n_clust
for (k in 2:K) {
  dists_to_existing <- D[, centre_loc, drop = FALSE]
  min_dists <- apply(dists_to_existing, 1, min)
  next_centre <- which.max(min_dists)
  centre_loc <- c(centre_loc, next_centre)
}

# range <- 0.1 # controls spatial decay
range <- quantile(as.vector(dist(coords)), 0.2) # "local-ish" neighbourhood scale
W <- exp(-D / range) # exponential kernel; W_ij ~ 1 if close, ~  0 if far
diag(W) <- 1 # include self-weight
S <- W / rowSums(W) # normalise to row-stochastic smoothing operator

rho_centres <- cor_t_vec # assign correlation values to centres

# function to map from correlation to latent z space (i.e. invert tanh mapping)
inv_rho_to_z <- function(rho, rho_min, rho_max, a) {
  # map rho in (rho_min, rho_max) -> x in (-1, 1)
  x <- (rho - rho_min) / (rho_max - rho_min)
  x <- 2 * x - 1

  # avoid atanh overflow if rho hits bounds
  eps <- 1e-8
  x <- pmin(pmax(x, -1 + eps), 1 - eps)

  atanh(x) / a
}

# parameters
rho_min <- 0.01
rho_max <- 0.99
# a <- 3 # gain parameter to control steepness of transition
a <- 3

# get latent z values at centres
z_centres <- inv_rho_to_z(rho_centres, rho_min, rho_max, a)

# create smooth spatial random field
set.seed(seed)
z <- rnorm(n_locs) # initialise latent field
z[centre_loc] <- z_centres # impose anchors on latent scale

# n_iter <- 3
# for (t in 1:n_iter) {
#   z <- as.numeric(S %*% z)   # smooth
#   z[centre_loc] <- z_centres # re-impose anchors so they remain exact
# }
# z_smooth <- z

for (t in 1:500) { # try 50–500
  print(t)
  z_old <- z
  z <- as.numeric(S %*% z) # smooth
  z[centre_loc] <- z_centres # re-impose anchors so they remain exact
  if (max(abs(z - z_old)) < 1e-6) break # convergence check
}
scale_factor <- 2.5 # try 1.5, 2, 3, 4
z2 <- z * scale_factor
z2[centre_loc] <- z_centres
z_smooth <- z2

# map real line -> (-1, 1) with tanh, then scale into a practical range (0, 1)
# rho_loc <- rho_min + (rho_max - rho_min) * ((tanh(a * z_smooth) + 1) / 2)
rho_loc <- rho_min + (rho_max - rho_min) * ((tanh(a * z_smooth) + 1) / 2)
rho_loc

cbind(target = rho_centres, achieved = rho_loc[centre_loc])

# plot spatial field of correlations
loc_df <- data.frame(
  loc = paste0("loc_", 1:n_locs),
  x = coords[, 1],
  y = coords[, 2],
  rho = rho_loc,
  centre = ifelse(1:n_locs %in% centre_loc, "centre", "non-centre")
)

ggplot(loc_df, aes(x = x, y = y, colour = rho, shape = centre)) +
  geom_point(size = 4) +
  labs(
    x = "X",
    y = "Y"
  ) +
  cecl_theme(nejm_pal = FALSE) +
  scale_colour_viridis_c() +
  guides(shape = "none") +
  theme(
    legend.text = element_text(angle = 45, hjust = 1)
  )

# check that closer points have more similar rho
rho_diff <- as.vector(dist(rho_loc))
d_vals <- as.vector(dist(coords))
plot(d_vals, rho_diff, xlab = "Distance", ylab = "|rho_i - rho_j|")
abline(lm(rho_diff ~ d_vals), lwd = 2)
summary(lm(rho_diff ~ d_vals))

# generate data for each location
set.seed(seed)
dat_list <- lapply(1:n_locs, function(ell) {
  x <- gen_t_loc(rho_loc[ell], n_per_loc, n_vars = n_vars, df_t = df_t)
  x$name <- paste0("loc_", ell)
  x$rho <- rho_loc[ell]
  x
})
data_spatial <- do.call(rbind, dat_list)

# also plot again with correlations estimated from data to check
data_spatial |>
  group_by(name) |>
  summarise(
    cor = cor(X1, X2),
    rho = first(rho),
    .groups = "drop"
  ) |>
  left_join(select(loc_df, -rho), by = c("name" = "loc")) |>
  ggplot(aes(x = x, y = y, colour = rho)) +
  geom_point(size = 4) +
  labs(
    x = "X",
    y = "Y"
  ) +
  cecl_theme(nejm_pal = FALSE) +
  scale_colour_viridis_c() +
  theme(legend.text = element_text(angle = 45, hjust = 1))

















#### Option 2: spatially varying mixture of cluster correlations ####

K <- 3
# random cluster centers
# centre_loc <- sample(1:n_locs, K)
# centres <- coords[centre_loc, , drop = FALSE] # or fixed points
# choose cluster centres as far away points
# first, calculate distance matrix
D_coords <- as.matrix(dist(coords))
# choose centres as points that are maximally distant
centre_loc <- c(1)
for (k in 2:K) {
  dists_to_existing <- D_coords[, centre_loc, drop = FALSE]
  min_dists <- apply(dists_to_existing, 1, min)
  next_centre <- which.max(min_dists)
  centre_loc <- c(centre_loc, next_centre)
}

dist_to_centre <- sapply(1:K, \(k) sqrt(rowSums((coords - centres[k, ])^2)))
tau <- 0.1 # smaller -> sharper clusters
logw <- -dist_to_centre / tau
w <- exp(logw) / rowSums(exp(logw)) # softmax weights

rho_k <- c(0.1, 0.5, 0.9)
rho_loc <- as.numeric(w %*% rho_k)
rho_loc

# plot spatial field of correlations
loc_df <- data.frame(
  loc = paste0("loc_", 1:n_locs),
  x = coords[, 1],
  y = coords[, 2],
  rho = rho_loc,
  centre = ifelse(1:n_locs %in% centre_loc, "centre", "non-centre")
)

ggplot(loc_df, aes(x = x, y = y, colour = rho, shape = centre)) +
  geom_point(size = 4) +
  labs(
    x = "X",
    y = "Y"
  ) +
  cecl_theme(nejm_pal = FALSE) +
  scale_colour_viridis_c()

#### Spatially varying simulation ####

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

gen_t_loc <- function(rho, n, n_vars = 2, df_t = 3) {
  cop_t <- tCopula(param = rho, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}

#### Method 1 ####

# generate poitns for cluster centres
set.seed(seed)
clust_centres <- data.frame(
  cor = cor_t_vec,
  x = runif(n_clust, 0, 5),
  y = runif(n_clust, 0, 5)
)

# plot
ggplot(clust_centres, aes(x = x, y = y, color = as.factor(cor))) +
  geom_point(size = 4) +
  labs(color = "Correlation") +
  theme_minimal()

# generate locations around cluster centres

# smooth correlation function for these points, while keeping centres fixed

# how many locations per cluster
locs_per_clust <- n_locs / n_clust
stopifnot(locs_per_clust == floor(locs_per_clust))

sd_loc <- 1.0 # spatial spread around each centre (tune this)

set.seed(seed)
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

head(locs)

ggplot() +
  geom_point(data = locs, aes(x = x_loc, y = y_loc), alpha = 0.6) +
  geom_point(data = clust_centres, aes(x = x, y = y, colour = as.factor(cor)), size = 4) +
  theme_minimal() +
  labs(colour = "Centre corr")

# length scale (bigger => smoother, more overlap among centres)
ell <- 2.0

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

# keep rho inside (-1, 1) for copula validity (esp if you later allow negatives)
eps <- 1e-6
rho_loc <- pmin(1 - eps, pmax(-1 + eps, rho_loc))

locs$rho <- rho_loc

ggplot(locs, aes(x = x_loc, y = y_loc, colour = rho)) +
  geom_point(size = 2) +
  scale_colour_viridis_c() +
  theme_minimal() +
  labs(colour = "rho(x,y)")

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

# plot
ggplot(locs2, aes(x = x_loc, y = y_loc, colour = rho, shape = is_centre, size = is_centre)) +
  geom_point() +
  scale_colour_viridis_c() +
  theme_minimal() +
  labs(colour = "rho(x,y)", shape = "Is centre") +
  scale_size_manual(values = c(5, 8))

# for non-centres, tabulate min, max, mean rho
locs |>
  group_by(cluster) |>
  summarise(
    min_rho = min(rho),
    max_rho = max(rho),
    mean_rho = mean(rho)
  )

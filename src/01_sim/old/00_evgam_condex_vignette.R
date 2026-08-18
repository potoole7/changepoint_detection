#### Simulate data ####

set.seed(1)

n_sites <- 40
n_days <- 365 * 5 # 5 years

locations <- data.frame(
  lon = runif(n_sites, -125, -65), # US-ish box
  lat = runif(n_sites, 25, 50)
)

dates <- seq.Date(
  from = as.Date("2010-01-01"),
  by = "day",
  length.out = n_days
)

library(fields) # for rdist

coords <- as.matrix(locations)
D <- rdist(coords)

# Exponential spatial covariance
range_param <- 10
Sigma <- exp(-D / range_param)

# Cholesky for fast simulation
L <- chol(Sigma + 1e-6 * diag(n_sites))

library(evd) # for Laplace transform helpers

day_SFT_max <- matrix(NA, n_days, n_sites)

for (t in 1:n_days) {
  # Smooth seasonal cycle
  seasonal <- 2 * sin(2 * pi * t / 365)

  # Spatially correlated noise
  z <- rnorm(n_sites)
  spatial_field <- L %*% z

  # Gaussian temp
  temp_gauss <- seasonal + spatial_field

  # Convert Gaussian → Uniform → Laplace
  u <- pnorm(temp_gauss)
  temp_laplace <- qlaplace(u)

  day_SFT_max[t, ] <- temp_laplace
}

dir.create("Average surface skin temperature", showWarnings = FALSE)

years <- format(dates, "%Y")

for (yr in unique(years)) {
  idx <- which(years == yr)

  day_SFT_max_year <- day_SFT_max[idx, ]
  dates_year <- dates[idx]

  save(
    dates_year,
    locations,
    day_SFT_max_year,
    file = paste0(
      "Average surface skin temperature/sim_MERRA2_", yr, ".RData"
    )
  )
}

#### Run vignette ####

## ============================================================
## Conditional extremes with evgam — one-go script
## Source: https://byoungman.github.io/evgam/condex/index.html
## Adds: synthetic "MERRA-2-like" data generator if files missing
## ============================================================

## ---- user options ------------------------------------------
# RUN_COMPOSITE_FIT <- FALSE     # set TRUE to run Section 2.3 (can be slow)
RUN_COMPOSITE_FIT <- TRUE
DATA_DIR <- "Average surface skin temperature"
SEED <- 1

## ---- packages ----------------------------------------------
req_pkgs <- c("evgam", "mgcv", "lattice")
to_install <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(evgam)
library(lattice)

## ---- helper: Laplace quantile (same as vignette) -----------
qlaplace_local <- function(p) {
  ifelse(p <= .5, log(2 * p), -log(2 * (1 - p)))
}

## ---- helper: generate synthetic yearly .RData files --------
generate_synthetic_files <- function(
    out_dir = DATA_DIR,
    n_sites = 40,
    start_date = as.Date("2010-01-01"),
    n_days = 365 * 5,
    seed = SEED) {
  set.seed(seed)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # locations: lon/lat matrix like vignette expects
  locations <- cbind(
    lon = runif(n_sites, -125, -65),
    lat = runif(n_sites, 25, 50)
  )

  dates_all <- seq.Date(from = start_date, by = "day", length.out = n_days)

  # build a spatial correlation using an RBF kernel
  coords <- locations
  D <- as.matrix(dist(coords))
  range_param <- 8
  Sigma <- exp(-(D^2) / (2 * range_param^2))
  L <- chol(Sigma + 1e-6 * diag(n_sites))

  # simulate Gaussian field + seasonality, then map to Laplace
  day_SFT_max_all <- matrix(NA_real_, nrow = n_days, ncol = n_sites)
  for (t in seq_len(n_days)) {
    seasonal <- 0.7 * sin(2 * pi * t / 365) + 0.25 * cos(2 * pi * t / 365)
    z <- drop(L %*% rnorm(n_sites))
    g <- seasonal + z

    # Gaussian -> Uniform -> Laplace
    u <- pnorm(g)
    y <- qlaplace_local(u)

    day_SFT_max_all[t, ] <- y
  }

  # write yearly files containing: dates, locations, day_SFT_max
  yrs <- format(dates_all, "%Y")
  for (yr in unique(yrs)) {
    idx <- which(yrs == yr)
    dates <- dates_all[idx]
    # NOTE: vignette later does t(rbind(...)) so storing as days x sites is fine
    day_SFT_max <- day_SFT_max_all[idx, , drop = FALSE]

    save(dates, locations, day_SFT_max,
      file = file.path(out_dir, paste0("sim_MERRA2_", yr, ".RData"))
    )
  }

  invisible(TRUE)
}

## ---- ensure data exists ------------------------------------
fls <- list.files(DATA_DIR, full.names = TRUE)
if (length(fls) == 0) {
  message("No data files found in '", DATA_DIR, "'. Generating synthetic files...")
  generate_synthetic_files(out_dir = DATA_DIR)
  fls <- list.files(DATA_DIR, full.names = TRUE)
}

## ---- load & combine yearly files (as vignette) -------------
lst <- list()
for (i in seq_along(fls)) {
  load(fls[i]) # expects dates, locations, day_SFT_max
  lst[[i]] <- list(dates = dates, locations = locations, day_SFT_max = day_SFT_max)
}

dates <- do.call(c, lapply(lst, "[[", "dates"))
day_SFT_max <- t(do.call(rbind, lapply(lst, "[[", "day_SFT_max")))

## ---- transform to Laplace scale (as vignette) --------------
day_SFT_max_F <- t(apply(day_SFT_max, 1, function(x) (rank(x) - .5) / length(x)))
data_L <- qlaplace_local(day_SFT_max_F)

## ---- metadata in data.frame format (as vignette) -----------
grid_df <- data.frame(lon = locations[, 1], lat = locations[, 2])
d <- nrow(grid_df)

## ============================================================
## 2.2 Single-site conditioning (as vignette)
## ============================================================

cond <- with(as.data.frame(locations), which(lon == median(lon) & lat == median(lat)))
if (length(cond) == 0) cond <- round(nrow(grid_df) / 2) # fallback if no exact median match

lon0 <- locations[cond, 1]
lat0 <- locations[cond, 2]

threshold <- qlaplace_local(0.8)
exceed <- which(data_L[cond, ] > threshold)

cond_mat <- matrix(data_L[cond, ], nrow(data_L), ncol(data_L), byrow = TRUE)
ev_df_mat <- data.frame(lon = grid_df$lon[-cond], lat = grid_df$lat[-cond])
ev_df_mat$y <- data_L[-cond, exceed]
ev_df_mat$x <- cond_mat[-cond, exceed]

fmla_condex <- list(
  y ~ s(lon, lat, k = 10),
  ~1,
  ~ s(lon, lat, k = 10),
  ~ s(lon, lat, k = 10)
)

m1 <- evgam(fmla_condex,
  data = ev_df_mat, family = "condex",
  args = list(x = ev_df_mat$x)
)

pred1 <- grid_df
pred1 <- cbind(pred1, predict(m1, newdata = grid_df, type = "response"))

panel_fn <- function(...) {
  panel.levelplot(...)
  panel.xyplot(lon0, lat0, pch = 19, col = 1, cex = 1.5)
}

l1 <- levelplot(alpha ~ lon * lat, pred1, main = "alpha", panel = panel_fn)
l2 <- levelplot(beta ~ lon * lat, pred1, main = "beta", panel = panel_fn)
l3 <- levelplot(location ~ lon * lat, pred1, main = "mu", panel = panel_fn)
l4 <- levelplot(scale ~ lon * lat, pred1, main = "sigma", panel = panel_fn)

print(l1, split = c(1, 1, 2, 2), more = TRUE)
print(l2, split = c(2, 1, 2, 2), more = TRUE)
print(l3, split = c(1, 2, 2, 2), more = TRUE)
print(l4, split = c(2, 2, 2, 2))

## ============================================================
## 2.3 Composite fitting for multiple conditioning sites (optional)
## ============================================================

if (isTRUE(RUN_COMPOSITE_FIT)) {
  ncond <- 5
  cond_vec <- round(nrow(grid_df) * (1:ncond - .5) / ncond)

  ncol <- max(rowSums(data_L[cond_vec, ] > threshold))
  mat0 <- matrix(NA, nrow(data_L) - 1, ncol)

  ev_df_list <- list()
  for (i in 1:ncond) {
    cond <- cond_vec[i]
    lon0 <- grid_df[cond, "lon"]
    lat0 <- grid_df[cond, "lat"]

    ev_df_i <- data.frame(lon = grid_df$lon[-cond], lat = grid_df$lat[-cond])
    ev_df_i$lon2 <- ev_df_i$lon - lon0
    ev_df_i$lat2 <- ev_df_i$lat - lat0

    matx <- maty <- mat0
    exceed <- which(data_L[cond, ] > threshold)

    maty[, seq_along(exceed)] <- data_L[-cond, exceed]
    ev_df_i$y <- maty

    cond_mat <- matrix(data_L[cond, ], nrow(data_L), ncol(data_L), byrow = TRUE)
    matx[, seq_along(exceed)] <- cond_mat[-cond, exceed]
    ev_df_i$x <- matx

    ev_df_i$weights <- 0 * maty + 1 / ncond
    ev_df_list[[i]] <- ev_df_i
  }

  ev_df_comp <- do.call(rbind, lapply(ev_df_list, function(x) x[, c("lon", "lat", "lon2", "lat2")]))
  ev_df_comp$y <- do.call(rbind, lapply(ev_df_list, function(x) x$y))
  ev_df_comp$x <- do.call(rbind, lapply(ev_df_list, function(x) x$x))
  ev_df_comp$weights <- do.call(rbind, lapply(ev_df_list, function(x) x$weights))

  fmla_condex2 <- list(
    y ~ s(lon2, lat2, k = 10),
    ~1,
    ~ s(lon2, lat2, k = 10),
    ~ s(lon2, lat2, k = 10)
  )

  m2 <- evgam(
    fmla_condex2,
    data = ev_df_comp, family = "condex",
    args = list(x = ev_df_comp$x, weights = ev_df_comp$weights)
  )

  cond <- round(nrow(grid_df) / 2)
  lon0 <- grid_df[cond, "lon"]
  lat0 <- grid_df[cond, "lat"]

  pred2 <- grid_df
  pred2$lon2 <- pred2$lon - lon0
  pred2$lat2 <- pred2$lat - lat0
  pred2 <- cbind(pred2, predict(m2, newdata = pred2, type = "response"))

  panel_fn2 <- function(...) {
    panel.levelplot(...)
    panel.xyplot(0, 0, pch = 19, col = 1, cex = 1.5)
  }

  l1 <- levelplot(alpha ~ lon2 * lat2, pred2, main = "alpha", panel = panel_fn2)
  l2 <- levelplot(beta ~ lon2 * lat2, pred2, main = "beta", panel = panel_fn2)
  l3 <- levelplot(location ~ lon2 * lat2, pred2, main = "mu", panel = panel_fn2)
  l4 <- levelplot(scale ~ lon2 * lat2, pred2, main = "sigma", panel = panel_fn2)

  print(l1, split = c(1, 1, 2, 2), more = TRUE)
  print(l2, split = c(2, 1, 2, 2), more = TRUE)
  print(l3, split = c(1, 2, 2, 2), more = TRUE)
  print(l4, split = c(2, 2, 2, 2))
}

message("Done.")

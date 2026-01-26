#### Simulation example ####

# Questions:
# - Just to check, should the marginals stay the same? Is that realistic?
# - Difference due to randomness, particularly in generating Laplace
#   observations, is quite high. How to reduce this or factor it in?

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr)
library(ggplot2)


#### metadata ####

seed <- 123
n_vars <- 2 # number of variables
# n_locs <- 12 # number of locations
n_locs <- 60
n_clust <- 3 # number of clusters
n_per_loc <- 1000 # total obs per location times # locs per group
n <- n_per_loc * (n_locs / n_clust) # total obs per location times n_locs per group
df_t <- 3 # degrees of freedom for t-copula
cor_t_vec <- c(0.1, 0.5, 0.9) # (initial) correlation parameters for each group


#### Functions ####

# Function to generate multivariate t data with specified correlation
gen_t <- \(cor_t, n) {
  # generate data
  cop_t <- tCopula(param = cor_t, dim = n_vars, df = df_t, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}


#### Scenario 1: Local change ####

# increase lowest correlation by 0.05 until it reaches middle correlation
set.seed(seed)
# start from 0 to have control example (checks difference due to seed)
increments <- c(0, 0:8)
data_loc <- bind_rows(lapply(seq_along(increments), \(i) {
  cor_spec <- cor_t_vec
  cor_spec[[1]] <- cor_spec[[1]] + (0.05 * increments[[i]])
  print(cor_spec)
  # ret <- do.call(rbind, lapply(cor_t_vec, gen_t, n = n))
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



#### Fit CE and cluster ####

# Calculate CE and cluster for all time blocks:

# We know the exact GPD values, so use these to transform to Laplace
locs <- unique(data_loc$name)
vars <- names(data_loc)[1:n_vars]
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
data_laplace_loc <- data_loc |>
  mutate(block = factor(block, levels = unique(block))) |>
  group_split(block) |>
  lapply(\(x) {
    ret <- trans_marg(gpd_vals, x, vars = vars)
    names(ret) <- locs
    ret
  })
names(data_laplace_loc) <- unique(data_loc$block)

# check you've transformed correctly
par(mfrow = c(1, 3))
plot(data_laplace_loc[[1]][[1]])
plot(data_laplace_loc[[2]][[1]])
plot(data_laplace_loc[[10]][[1]])
par(mfrow = c(1, 1))

marg_loc <- lapply(data_laplace_loc, as_cecl_marg)

# fit CE model to each
dep_loc <- lapply(marg_loc, cecl_dep, cond_prob = 0.9, fit_no_keef = TRUE)

# Calculate divergence and cluster
dist_loc <- lapply(seq_along(dep_loc), \(i) {
  cecl_dist(dep_loc[[i]], marg_loc[[i]], seed = seed)
})
clust_loc <- lapply(
  dist_loc,
  cecl_clust,
  k = 3,
  cluster_mem = sort(rep(1:n_clust, each = n_locs / n_clust)),
  seed = seed
)

#### Analysis ####

# plot distance matrices
dist_mat_plots <- lapply(
  dist_loc, ggplot,
  which = "image", col_breaks = seq(0, 0.7, by = 0.1)
)
dist_mat_plots[[1]] # clusters distinct, cluster 1 further from 3 than 2
dist_mat_plots[[5]] # cluster 1 getting closer to 2
dist_mat_plots[[10]] # cluster 1 and 2 indistinguishable

# Check clustering solution accuracy for each
sapply(clust_loc, `[[`, "adj_rand") # Gradually worsens from 1 to 0.5

# Calculate difference in divergence matrices:
# difference between neighbouring div matrices
# dist_diff <- clust_loc[1:(length(clust_loc) - 1)]
dist_diff <- dist_loc[1:(length(clust_loc) - 1)]
for (i in seq_along(dist_diff)) {
  dist_diff[[i]]$dist_mat <- abs(
    dist_loc[[i]]$dist_mat - dist_loc[[i + 1]]$dist_mat
  )
}

# difference between first and subsequent divs
# dist_diff_first <- clust_loc[1:length(clust_loc)]
dist_diff_first <- dist_loc[1:length(dist_loc)]
for (i in seq_along(dist_diff_first)) {
  print(i)
  dist_diff_first[[i]]$dist_mat <- abs(
    dist_loc[[1]]$dist_mat -
      dist_loc[[i]]$dist_mat
  )
}

# plot
dist_diff_plots <- lapply(
  dist_diff,
  ggplot,
  which = "image",
  # show_xlab = TRUE,
  col_breaks = seq(0, 0.4, by = 0.1)
)
# Interpretation:
dist_diff_plots[[1]]
# control: Noticeable differences here due to randomness!
# Interesting given same clustering solution
# Some particularly large differences for 1st vs 3rd cluster comparisons
# TODO Investigate, why is comparison of 1st cluster to others so different?
# Maybe look at these local changes more closely
dist_diff_plots[[2]]
# intermediate change: Doesn't look a whole lot different from last plot!
#
# Maybe slightly larger (colour scale goes up to 0.3 vs 0.25 before)
dist_diff_plots[[9]]
# Final change: Much smaller differences, as percentage increase in
# correlation smaller
# Also no difference for comparing clusters 1 and 2, which is strange
# given that for the first plot there was a noticeable difference there
# which should be just randomness?

# TODO fix colour bar
# debugonce(ggplot.cecl_dist)

dist_diff_first_plots <- lapply(
  dist_diff_first,
  ggplot,
  which = "image",
  col_breaks = seq(0, 0.5, by = 0.1)
)
# Interpretation:
dist_diff_first_plots[[1]] # all 0, as expected
dist_diff_first_plots[[2]]
# control: Surprisingly high diff for what should just be randomness? Max of ~0.25
# TODO Investigate!!
dist_diff_first_plots[[5]]
# Intermediate: Definitely larger than control, max ~0.4
# Differences larger for comparing cluster 1 to others, as expected
dist_diff_first_plots[[10]]
# Max: similar to intermediate, but larger max ~0.5

# Find max-norm for each time step to each subsequent (and first to each)
# Interpret
sapply(dist_diff, \(x) norm(as.matrix(x$dist_mat), type = "F"))
# Gets smaller as percentage increase in correlation smaller
# Largest for first increase in correlation, but still quite large for
# just the control
sapply(dist_diff, \(x) norm(as.matrix(x$dist_mat), type = "I"))
# Seems to capture signal better, larger difference between control and
# first increase in correlation

sapply(dist_diff_first, \(x) norm(as.matrix(x$dist_mat), type = "F"))
sapply(dist_diff_first, \(x) norm(as.matrix(x$dist_mat), type = "I"))
# both increasing, difference between randomness and actual change clearer
# for infinity/max norm

# TODO Quantify difference in clustering solutions via ARI
clust_diff <- vector(length = length(clust_loc) - 1)
for (i in seq_along(clust_diff)) {
  clust_diff[[i]] <- mclust::adjustedRandIndex(
    clust_loc[[i]]$pam$clustering,
    clust_loc[[i + 1]]$pam$clustering
  )
}

# difference between first and subsequent divs
clust_diff_first <- vector(length = length(clust_loc))
for (i in seq_along(clust_diff_first)) {
  clust_diff_first[[i]] <- mclust::adjustedRandIndex(
    clust_loc[[1]]$pam$clustering,
    clust_loc[[i]]$pam$clustering
  )
}

# all decreasing, clustering solution for control is the same so no difference
# there
sapply(clust_loc, `[[`, "adj_rand") # ARI to true clusters
clust_diff
clust_diff_first


#### Scenario 2: Global change ####

# TODO Bring all correlations closer to 1?

#### Analysis ####

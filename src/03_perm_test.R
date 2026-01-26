#### Permutation Test ####

# TODO Add method for computing norms for two distance matrices
# Expand to infinity norm (done)
# TODO Also do for clustering solution
# TODO Change fun to run 1 loop; can parallelise outside (to time runs) (done)

#### libs ####

devtools::load_all("../CeCl")
library(copula) # for generating data
library(dplyr)
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

diff <- 0.4 # difference in correlation for local change
# diff <- 0.05
# diff <- 0.1
# diff <- 0.15
# diff <- 0.16 # successfully picked up with F norm (p-value = 0.047)
diff <- 0
mc_cores <- detectCores() - 1 # number of cores for parallel processing

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

#### Fit CE, calculate divergence and compare ####

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
par(mfrow = c(1, 2))
plot(data_laplace_loc[[1]][[1]])
plot(data_laplace_loc[[2]][[1]])
par(mfrow = c(1, 1))

# convert to cecl_marg objects
marg_loc <- lapply(data_laplace_loc, as_cecl_marg)

# fit CE model to each
# TODO Should cond_prob not be like 0 or something here?
dep_loc <- lapply(marg_loc, cecl_dep, cond_prob = 0.9, fit_no_keef = TRUE)

# Calculate divergence
dist_loc <- lapply(seq_along(dep_loc), \(i) {
  cecl_dist(dep_loc[[i]], marg_loc[[i]], seed = seed)
})
ggplot(dist_loc[[1]], which = "image", col_breaks = seq(0, 0.5, by = 0.1))
ggplot(dist_loc[[2]], which = "image", col_breaks = seq(0, 0.5, by = 0.1))

dist_diff <- dist_loc[[1]]
dist_diff$dist_mat <- abs(dist_loc[[2]]$dist_mat - dist_loc[[1]]$dist_mat)
image(as.matrix(dist_diff$dist_mat))

norm(as.matrix(dist_diff$dist_mat), type = "F")
norm(as.matrix(dist_diff$dist_mat), type = "I")

sort(rowSums(as.matrix(dist_diff$dist_mat)))



# TODO Fix
# ggplot(dist_diff[[1]], which = "image")

# calculate norm of difference between time blocks
(norm_diff_loc_frob <- norm(
  as.matrix(dist_loc[[2]]$dist_mat - dist_loc[[1]]$dist_mat),
  type = "F"
))
(norm_diff_loc_inf <- norm(
  as.matrix(dist_loc[[2]]$dist_mat - dist_loc[[1]]$dist_mat),
  type = "I"
))

# also calculate difference in clustering solutions
clust_loc <- lapply(
  dist_loc,
  cecl_clust,
  k = n_clust,
  cluster_mem = sort(rep(1:n_clust, each = n_locs / n_clust)),
  seed = seed
)

(norm_diff_loc_clust <- mclust::adjustedRandIndex(
  clust_loc[[1]]$pam$clustering,
  clust_loc[[2]]$pam$clustering
))

# Question: Is this norm significantly larger than expected under null
# hypothesis of no change?
# To test this, we can perform a permutation test!

#### Permutation Test ####

# Question:
# If observations within each location were randomly assigned to the two
# blocks, how large would the difference in CE-based dependence matrices
# typically be?

# TODO Parallelise this
n_perm <- 100 # number of permutations

# First, combined transformed data into single data frame (to permute over)
# (As we know the exact marginals, can permutate transformed data directly)
data_laplace_combined <- bind_rows(lapply(seq_along(data_laplace_loc), \(i) {
  bind_rows(lapply(data_laplace_loc[[i]], as.data.frame), .id = "name") %>%
    mutate(block = i)
})) |>
  select(X1, X2, name, block)

perm_test <- \(
  type = c("norm", "clustering"), norm_type = "F"
) {
  type <- match.arg(type, several.ok = FALSE)
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

  # convert to cecl_marg objects
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

  # fit CE models
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

  if (type == "norm") {
    # norm of difference
    return(norm(
      as.matrix(dist_perm[[2]]$dist_mat -
        dist_perm[[1]]$dist_mat),
      type = norm_type
    ))
  } else if (type == "clustering") {
    # clustering solution
    clust_perm <- lapply(
      dist_perm,
      cecl_clust,
      k = n_clust,
      cluster_mem = sort(rep(1:n_clust, each = n_locs / n_clust)),
      seed = seed
    )
    # calculate difference in clustering solutions
    return(mclust::adjustedRandIndex(
      clust_perm[[1]]$pam$clustering,
      clust_perm[[2]]$pam$clustering
    ))
  }
}

# time runs
# microbenchmark::microbenchmark(
#   "norm" = perm_test(type = "norm", norm_type = "F"),
#   "clust" = perm_test(type = "clustering"),
#   times = 100
# )
# Unit: milliseconds
# expr       min       lq     mean   median       uq      max neval cld
# norm  522.4950 552.1078 627.4973 575.6697 630.6949 1323.590   100   a
# clust 510.4279 551.8377 630.4814 570.8510 625.9871 1355.258   100   a

# ~6s per run == ~14/15 minutes for 1000 permutations on 7 cores

# calculate permutation norms
set.seed(seed)
perm_norms_frob <- unlist(mclapply(seq_len(n_perm), \(i) {
  perm_test(type = "norm", norm_type = "F")
}))
perm_norms_inf <- unlist(mclapply(seq_len(n_perm), \(i) {
  perm_test(type = "norm", norm_type = "F")
}))
perm_norms_clust <- unlist(mclapply(seq_len(n_perm), \(i) {
  perm_test(type = "clustering")
}))

par(mfrow = c(3, 1))
hist(perm_norms_frob, xlim = c(0, max(c(norm_diff_loc_frob, perm_norms_frob))))
abline(v = norm_diff_loc_frob, col = "red")
hist(perm_norms_inf, xlim = c(0, max(c(norm_diff_loc_inf, perm_norms_inf))))
abline(v = norm_diff_loc_inf, col = "red")
hist(
  perm_norms_clust,
  xlim = c(0, max(c(norm_diff_loc_clust, perm_norms_clust)))
)
abline(v = norm_diff_loc_clust, col = "red")
par(mfrow = c(1, 1))

# Calculate p-value
(p_value_frob <- mean(perm_norms_frob >= norm_diff_loc_frob))
(p_value_inf <- mean(perm_norms_inf >= norm_diff_loc_inf))
(p_value_clust <- mean(perm_norms_clust >= norm_diff_loc_clust))

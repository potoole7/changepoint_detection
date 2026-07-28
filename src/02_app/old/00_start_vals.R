#### Enable starting values with CeCl ####

# Test:
# (i) How similar parameter values are when we use start values
# (ii) If they are dissimilar to

# library(CeCl)
devtools::load_all("../CeCl")
library(copula)
library(dplyr)
library(ggplot2)

## metadata
n_locs <- 20
dep_val <- qlaplace(0.9) # dependence quantile
seed <- 123

# function to generate multivariate t data with specified correlation
gen_t <- \(cor_t, n_vars = 2, n = 10000) {
  # generate data
  cop_t <- tCopula(param = cor_t, dim = n_vars, df = 3, dispstr = "ex")
  u <- rCopula(n, cop_t)
  data <- data.frame(apply(u, 2, qgpd, xi = -0.05, sigma = 1, u = 0))
}

# generate data for two different correlation settings
set.seed(123)
data <- rbind(
  gen_t(cor_t = 0.2), # first locations have low correlation
  gen_t(cor_t = 0.8) # last locations have high correlation
)
# Add dummy location names
# data$name <- rep(paste0("loc_", 1:10), each = 200)
data <- data |>
  mutate(name = rep(paste0("loc_", 1:n_locs), each = nrow(data) / n_locs))

# check
table(data$name)

# Fit marginal model
marg <- cecl_marg(
  data,
  # thresh_method = "quantile",
  thresh_method = "value",
  # thresh_args = 0.9,
  thresh_args = 0, # don't threshold, since we know all data is extreme
  marg_method = "ecdf"
)

# Fit dependence model
dep <- cecl_dep(obj = marg, cond_val = dep_val)

# Extract starting values for dependence parameters
start <- coef(dep)

# fit again with the same data
# dep_same <- cecl_dep(obj = marg, cond_val = dep_val, start = start)
dep_same <- cecl_dep(obj = marg, cond_prob = 0.9, nruns = 2)

# different with new starting values!
# TODO Improve how nruns works!
coef(dep)[c("a", "b")] - coef(dep_same)[c("a", "b")]

# run again with slightly different data, using start values from before
# (diffeerent seed)
data2 <- rbind(
  gen_t(cor_t = 0.2),
  gen_t(cor_t = 0.8)
) |>
  mutate(name = rep(paste0("loc_", 1:n_locs), each = n() / n_locs))
marg2 <- cecl_marg(
  data,
  thresh_method = "value",
  thresh_args = 0,
  marg_method = "ecdf"
)

# debugonce(cecl_dep)
# devtools::load_all("../CeCl")
dep2 <- cecl_dep(obj = marg2, cond_val = dep_val, start = start)

# check that estimates are fairly similar
# TODO Doesn't seem that small???
# Largest for m and s,
mat_dep <- as.matrix(coef(dep)[, c("a", "b", "m", "s", "ll", "dth")])
mat_dep2 <- as.matrix(coef(dep2)[, c("a", "b", "m", "s", "ll", "dth")])
# norm(as.matrix(coef(dep2)[, c("a", "b", "m", "s", "ll", "dth")] - start[, c("a", "b", "m", "s", "ll", "dth")]))

mat_dep - mat_dep2
image(mat_dep - mat_dep2)

apply(mat_dep - mat_dep2, 2, max)

norm(mat_dep - mat_dep2)
norm(mat_dep[, c("a", "b")] - mat_dep2[, c("a", "b")])

# check that fitting is slightly faster after first time
# (should be, as starting values are closer to optimum)
microbenchmark::microbenchmark(
  cecl_dep(obj = marg, cond_prob = 0.9),
  cecl_dep(obj = marg2, cond_prob = 0.9, start = start),
  times = 10
)
# min       lq     mean   median       uq      max neval cld
# 323.8326 344.9367 388.7484 382.6731 418.3398 716.8430    50  a
# 181.9620 210.0692 232.4247 234.3752 248.7957 289.9436    50   b
# Summary: significantly faster with starting values!

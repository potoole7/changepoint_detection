#### Compare chronology and climate data at Hurricane Ridge ####

# TODO Do outside mean field chronology for better extremes analysis

#### libs ####

library(dplR)
library(dplyr)
library(tidyr)
library(mgcv)
library(evgam)
library(ggplot2)

#### Load Data ####

data(wa082) # example tree-ring width data from dplR package
# data with extremes analysis on chronology
wa082_extremes <- readr::read_csv(
  "data/wa082_chronology_extremes.csv"
)
# Load climate data (GHCND data for airport near Hurricane Ridge)
climate_data <- readr::read_csv(
  "data/climate_data_GHCND_USC00456624_1933_1983.csv"
)

#### Combine Data ####

# first, need to aggregate climate data to annual values
annual_climate <- climate_data %>%
  mutate(
    year = as.numeric(format(as.Date(date), "%Y"))
  ) %>%
  pivot_wider(names_from = datatype, values_from = value) %>%
  janitor::clean_names() %>%
  group_by(year) %>%
  summarise(
    prcp = sum(prcp, na.rm = TRUE),
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    .groups = "drop"
  )

# Now combine with extremes data
combined_data <- wa082_extremes %>%
  inner_join(annual_climate, by = c("year")) |>
  # calculate for lagged variables as well
  mutate(
    across(c(prcp, tmax, tmin), \(x) lag(x, 1), .names = "{col}_lag1"),
    across(c(prcp, tmax, tmin), \(x) lag(x, 2), .names = "{col}_lag2"),
    across(c(prcp, tmax, tmin), \(x) lag(x, 3), .names = "{col}_lag3")
  )

# any years with missing TRI? Nope!
any(is.na(combined_data$index))

# save combined data
readr::write_csv(
  combined_data,
  "data/wa082_chronology_climate_combined.csv"
)

#### Plot ####

# plot cross-correlation between tree-ring index and climate variables
ccf_prcp <- ccf(combined_data$index, combined_data$prcp, lag.max = 5, na.action = na.omit, main = "TRI vs Annual Precipitation")
ccf_tmax <- ccf(combined_data$index, combined_data$tmax, lag.max = 5, na.action = na.omit, main = "TRI vs Average Annual Max Temperature")
ccf_tmin <- ccf(combined_data$index, combined_data$tmin, lag.max = 5, na.action = na.omit, main = "TRI vs Average Annual Min Temperature")
# For precipitaiton, negative ACF, highest in lag 3
# For max temperature, positive ACF, highest in lag 2
# No significant lags for min temperature
# These larger/longer lags seem to be a bit strange?

# Plot simple correlation over years (not lagged)
combined_data %>%
  # select(year, index, prcp, tmax, tmin, contains("lag1"), prcp_lag3, tmax_lag2, tmin_lag2) |>
  select(year, index, prcp, tmax, tmin, contains("lag")) |>
  # pivot_longer(cols = c(prcp, tmax, tmin), names_to = "variable", values_to = "climate") |>
  # pivot_longer(cols = c(prcp, tmax, tmin, prcp_lag1, tmax_lag1, tmin_lag1, prcp_lag3, tmax_lag2, tmin_lag2), names_to = "variable", values_to = "climate") |>
  pivot_longer(cols = c(prcp, tmax, tmin, prcp_lag1, tmax_lag1, tmin_lag1, prcp_lag2, prcp_lag3, tmax_lag2, tmax_lag3, tmin_lag2, tmin_lag2, tmin_lag3), names_to = "variable", values_to = "climate") |>
  ggplot(aes(x = climate, y = index)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  labs(x = "Climate variable", y = "Tree-ring index")
# higher temperature seems to be particularly correlated with a higher TRI
# correlation does seem higher when lagged, matches with ACF plots

# TODO Plot correlations in extremes (chi, etc)
par(mfrow = c(1, 2))
plot(texmex::chi(select(combined_data, index, prcp)))
plot(texmex::chi(select(combined_data, index, prcp_lag1)))
plot(texmex::chi(select(combined_data, index, prcp_lag2))) # actually shows chi!
plot(texmex::chi(select(combined_data, index, tmax)))
plot(texmex::chi(select(combined_data, index, tmax_lag1)))
plot(texmex::chi(select(combined_data, index, tmax_lag2)))
plot(texmex::chi(select(combined_data, index, tmax_lag3))) # same!
plot(texmex::chi(select(combined_data, index, tmin)))
plot(texmex::chi(select(combined_data, index, tmin_lag1)))
plot(texmex::chi(select(combined_data, index, tmin_lag2)))
plot(texmex::chi(select(combined_data, index, tmin_lag3)))


#### Simple modelling ####

summary(lm(index ~ tmax + tmin + year, data = combined_data)) # R^2: 0.1558
summary(lm(index ~ prcp_lag1 + tmax_lag1 + tmin_lag1 + year, data = combined_data)) # R^2: 0.179
summary(lm(index ~ prcp_lag3 + tmax_lag2 + tmin_lag2 + year, data = combined_data)) # R^3: 0.253
# model with lagged versions does best, but not particularly well

# now model precipitation as function of RWI
summary(lm(prcp ~ index + year, data = combined_data)) # R^2: 0.1296
summary(lm(prcp_lag3 ~ index + year, data = combined_data)) # R^2: 0.1878

# try spline
summary(gam(index ~ s(prcp) + s(tmax) + s(tmin) + s(year), data = combined_data)) # R^2: 0.262
summary(gam(index ~ s(prcp_lag3) + s(tmax_lag2) + s(tmin_lag2) + s(year), data = combined_data))
# R^2: 0.253, so not improving on linear model

# try evgam for just modelling extremes
thresh_dat <- \(dat, var, q) {
  var_vals <- dat[[var]]
  thresh <- quantile(var_vals, q, na.rm = TRUE)
  ret <- dat %>%
    filter(!is.na(var_vals) & var_vals >= thresh)
  ret[[var]] <- ret[[var]] - thresh
  ret
}

mod1 <- evgam(
  list(prcp_lag3 ~ s(index) + s(year)),
  data = thresh_dat(combined_data, "prcp_lag3", 0.45)
)
# mostly negative shape
predict(mod1, thresh_dat(combined_data, "prcp_lag3", 0.45))
mod2 <- evgam(
  list(tmax_lag2 ~ s(index) + s(year)),
  data = thresh_dat(combined_data, "tmax_lag2", 0.55)
)
# quite strong positive shape for most
predict(mod2, thresh_dat(combined_data, "tmax_lag2", 0.55))

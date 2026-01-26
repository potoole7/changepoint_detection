#### dplR Workshop Example ####

#### libs ####

library(dplR)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)

#### Load Data ####

# data for raw ring widths for Pacific silver fir Abies amabilis at
# Hurricane Ridge in Washington, USA (23 series covering 286 years)
data(wa082)

#### Metadata ####

# Site_Information
#   Site_Name: Hurricane Ridge
#   Location: Washington
#   Northernmost_Latitude: 47.98
#   Southernmost_Latitude: 47.98
#   Easternmost_Longitude: -123.47
#   Westernmost_Longitude: -123.47
#   Elevation_m: 1550
# Data_Collection
#   Collection_Name: WA082
#   First_Year: 1698
#   Last_Year: 1983
#   Time_Unit: CE
# Species
#   Species_Name: Abies amabilis Douglas ex J. Forbes
#   Common_Name: Pacific silver fir
#   Tree_Species_Code: ABAM

#### Investigate ####

# TODO How is year pulled from this object??
plot(wa082, plot.type = "spag") # spaghetti plot of raw widths

# detrend using age-dependent spline method
# NOTE: May not want to do this for extremes?
# TODO Read other papers to see best practices for extremes
wa082RWI <- detrend(wa082, method = "AgeDepSpline")
class(wa082RWI) <- class(wa082) # set class to rwl
plot(wa082RWI, plot.type = "spag")

# create chronology from detrended ring-width indices
wa082Crn <- chron(wa082RWI)
str(wa082Crn)

# plot chronology with 30-year spline
plot(wa082Crn, add.spline = TRUE, nyrs = 30)

# create chronology with prewhitening
# each series whitened using AR model before chronology is computed
# giving residual chronology/white noise
wa082CrnResid <- chron(wa082RWI, prewhiten = TRUE)
str(wa082CrnResid)
plot(wa082CrnResid)

# truncate to depths > 4
wa082CrnTrunc <- subset(wa082Crn, samp.depth > 4)
plot(wa082CrnTrunc, add.spline = TRUE, nyrs = 30)

# truncate via the subsample signal strength (SSS)
# SSS is a measure of how well the chronology represents the common signal
wa082Ids <- autoread.ids(wa082)
sssThresh <- 0.85
wa082SSS <- sss(wa082RWI, wa082Ids)
yrs <- time(wa082)
yrCutoff <- max(yrs[wa082SSS < sssThresh])
ggplot() +
  geom_rect(
    aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = yrCutoff),
    fill = "darkred", alpha = 0.5
  ) +
  annotate(geom = "text", y = 1.5, x = 1725, label = "SSS < 0.85") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = yrs, y = wa082Crn$std)) +
  labs(x = "Year", y = "RWI") +
  theme_minimal()

# cutoff the raw data to only those years with SSS > 0.85, redo chronology
wa082RwlSSS <- wa082[wa082SSS > sssThresh, ]
wa082RwiSSS <- detrend(wa082RwlSSS, method = "AgeDepSpl")

wa082CrnSSS <- chron(wa082RwiSSS)
ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = time(wa082CrnSSS), y = wa082CrnSSS$std)) +
  geom_line(
    aes(
      x = time(wa082CrnSSS),
      y = caps(wa082CrnSSS$std, nyrs = 30)
    ),
    color = "darkred"
  ) +
  labs(x = "Year", y = "RWI") +
  theme_minimal()

#### Extremes analysis ####

# wa082_rwi <- detrend(rwl = wa082, method = "ModNegExp")
wa082_rwi <- detrend(wa082, method = "AgeDepSpline") # RWI = ring-width indices

# create chronology from detrended ring-width indices
wa082_crn <- chron(wa082_rwi)
plot(wa082_crn, add.spline = TRUE, nyrs = 30)

site_chron_df <- data.frame(
  year = as.numeric(row.names(wa082_crn)),
  index = wa082_crn$std, # use 'std' column
  samp.depth = wa082_crn$samp.depth # number of trees contributing
)

# Define top/bottom 5% as extremes
upper_thresh <- quantile(site_chron_df$index, 0.95, na.rm = TRUE)
lower_thresh <- quantile(site_chron_df$index, 0.05, na.rm = TRUE)

site_chron_df <- site_chron_df %>%
  mutate(
    extreme_high = index >= upper_thresh,
    extreme_low  = index <= lower_thresh
  )

# View extreme years
site_chron_df %>%
  filter(extreme_high | extreme_low)

# Plot chronology with extremes highlighted
ggplot(site_chron_df, aes(x = year, y = index)) +
  geom_line() +
  geom_point(
    data = subset(site_chron_df, extreme_high),
    aes(x = year, y = index),
    color = "red", size = 2
  ) +
  geom_point(
    data = subset(site_chron_df, extreme_low),
    aes(x = year, y = index),
    color = "blue", size = 2
  ) +
  labs(
    title = "Tree-Ring Chronology with Extreme Years Highlighted",
    x = "Year",
    y = "Ring-Width Index"
  ) +
  theme_bw()

# Function to flag extremes per core
flag_extremes <- function(x) {
  upper <- quantile(x, 0.95, na.rm = TRUE)
  lower <- quantile(x, 0.05, na.rm = TRUE)
  data.frame(
    extreme_high = x >= upper,
    extreme_low  = x <= lower
  )
}

# Apply to all cores
extremes_list <- lapply(wa082_rwi, flag_extremes)

# Count how many trees are extreme each year
extreme_counts <- data.frame(
  year = as.numeric(row.names(wa082_rwi)),
  high_count = rowSums(sapply(extremes_list, function(df) df$extreme_high), na.rm = TRUE),
  low_count = rowSums(sapply(extremes_list, function(df) df$extreme_low), na.rm = TRUE)
)

# plot histogram
extreme_counts |>
  tidyr::pivot_longer(
    cols = c(high_count, low_count),
    names_to = "extreme_type",
    values_to = "count"
  ) |>
  ggplot(aes(x = year, y = count, fill = extreme_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("high_count" = "red", "low_count" = "blue"),
    labels = c("High Extremes", "Low Extremes")
  ) +
  labs(
    title = "Count of Extreme Ring-Width Indices per Year",
    x = "Year",
    y = "Number of Trees with Extreme Indices",
    fill = "Extreme Type"
  ) +
  theme_minimal()

# climate_annual: year, TAVG, PRCP, etc.
combined_df <- site_chron_df %>%
  left_join(extreme_counts, by = "year")

# save
readr::write_csv(
  combined_df,
  file = "data/wa082_chronology_extremes.csv"
)

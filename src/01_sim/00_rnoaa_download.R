#### Download NOAA Climate Data for a Given Station and Years ####

library(rnoaa)
library(lubridate)
library(dplyr, warn.conflicts = FALSE)

# Choose a station; COOP Port Angeles example:
station_id <- "GHCND:USC00456624" # Port Angeles COOP
# or airport:
# station_id <- "GHCND:USW00094266"

# find data available for this station
(data_avail <- ncdc_datasets(stationid = station_id)$data)
# available: GHCND, GSOM, GSOY, NORMAL_ANN, NORMAL_MLY, NORMAL_DLY
"GHCND" %in% data_avail$id

# get station metadata
(station <- ncdc_stations(stationid = station_id)$data)
# data available from 1933!

# loop through years to 1983
# round to nearest year
# temp: use fixed start date
# start_date <- as.Date("1980-01-01")
start_date <- lubridate::round_date(as.Date(station$mindate), unit = "year")
end_date <- as.Date("1983-12-31")
years <- seq(from = year(start_date), to = year(end_date), by = 1)

# function to download data for given year (max allowed in API call)
# download_data_year <- \(year, wait_time = 0) {
#   # to avoid hitting rate limits
#   if (wait_time > 0) {
#     Sys.sleep(wait_time)
#   }
#   tryCatch({
#     ncdc(
#       datasetid = "GHCND",
#       stationid = station_id,
#       datatypeid = c("TMAX", "TMIN", "PRCP"),
#       startdate = paste0(year, "-01-01"),
#       enddate   = paste0(year, "-12-31"),
#       limit = 5000
#     )$data
#   }, error = function(e){
#     message("Error for year ", year, ": ", e$message)
#     return(NULL)
#   })
# }
download_data_year <- \(year, datatypeid = c("TMAX", "TMIN", "PRCP"), wait_time = 1) {
  Sys.sleep(wait_time)
  all_data <- list()
  offset <- 1
  repeat {
    tmp <- tryCatch(
      {
        ncdc(
          datasetid = "GHCND",
          stationid = station_id,
          startdate = paste0(year, "-01-01"),
          enddate = paste0(year, "-12-31"),
          datatypeid = datatypeid,
          limit = 1000,
          offset = offset
        )$data
      },
      error = function(e) {
        message("Error fetching year ", year, " at offset ", offset, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(tmp) || nrow(tmp) == 0) break
    all_data[[length(all_data) + 1]] <- tmp
    if (nrow(tmp) < 1000) break # last page
    offset <- offset + 1000
  }
  bind_rows(all_data)
}

all_years <- lapply(years, function(y) {
  message("Downloading year ", y)
  download_data_year(y)
})

df_all_years <- bind_rows(all_years)
# check that all years are present
missing_years <- setdiff(
  as.character(c(years, last(years) + 1)),
  unique(substr(lubridate::round_date(as.Date(df_all_years$date), unit = "year"), 0, 4))
)
df_all_years2 <- df_all_years
while (length(missing_years) > 0) {
  df_missing <- lapply(missing_years, function(y) {
    message("Retrying download for missing year ", y)
    download_data_year(as.integer(y))
  })

  df_all_years2 <- bind_rows(df_all_years, bind_rows(df_missing)) |>
    arrange(date, datatype)
  (missing_years <- setdiff(
    as.character(c(years, last(years) + 1)),
    unique(substr(lubridate::round_date(as.Date(df_all_years2$date), unit = "year"), 0, 4))
  ))
}

# Transform to regular units (GHCND stores T in tenths of degrees C, P in tenths of mm)
df_clean <- df_all_years2 %>%
  mutate(
    value = case_when(
      datatype %in% c("TMAX", "TMIN") ~ value / 10, # to degrees C
      datatype == "PRCP" ~ value / 10, # to mm
      TRUE ~ value
    )
  ) %>%
  select(date, datatype, value)

# save
readr::write_csv(
  df_clean,
  file = paste0("data/climate_data_", gsub(":", "_", station_id), "_1933_1983.csv")
)

# to check against dendrochronology data:
# - visual comparison of time series
# - cross-correlation
# - linear model
# - mgcv analysis
# - evgam analysis
# - What else did other papers do???

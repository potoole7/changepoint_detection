#### Data Exploration ####


# TODO Look at extreme drought (i.e. lack of precipitation!)

# TODO Investigate missing data, especially for wind speeds

# 1. Plot means for 10 year intervals (done)
# Also plot differences from most recent decade (done)
# Do for Summer and Winter (done)
# 2. Do the same for 99th quantile (done)
# 3. Plot chi, chibar for temp vs rain & rain vs ws over entire period (done)
# Plot chi for temp vs rain & rain vs ws (over 10 year intervals) (done)
# plot differences (done)

# TODO Interpret!!!

# Other ideas:
# TODO - PACF/ACF plots? Or extremogram?

#### Libs ####

library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(lubridate)
library(parallel)
library(data.table)
library(sf)
library(ggplot2)
library(patchwork)
library(extRemes)
devtools::load_all("../CeCl")

#### Functions ####

# Function to plot various variables on map
map_plot <- \(
  x,
  areas,
  var,
  season = NULL,
  col_lab = NULL,
  title = NULL,
  fill_fun = scale_fill_viridis_c
) {
  p_dat <- x |>
    mutate(decade = paste0(decade, "s"))

  if (!is.null(season)) {
    p_dat <- filter(p_dat, season == !!season)
  }
  p <- p_dat |>
    ggplot() +
    geom_sf(data = areas, fill = NA, colour = "black") +
    geom_sf(aes(fill = .data[[var]]), shape = 21, size = 5, colour = "black") +
    facet_wrap(~decade, nrow = 2)

  if (!is.null(col_lab)) {
    p <- p + labs(fill = col_lab)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(fill_fun)) {
    p <- p + fill_fun()
  }
  p +
    CeCl::cecl_theme(legend.position = "right", nejm_pal = FALSE)
}

# Function to calculate differences from first decade (1960s)
calc_diff <- \(x, cols, decade = 1960) {
  baseline_cols <- paste0(cols, "_1960")

  x %>%
    left_join(
      x %>%
        filter(decade == !!decade) |>
        select(-decade) |>
        rename_with(
          ~ paste0(.x, "_1960"),
          all_of(cols)
        ) |>
        st_drop_geometry()
    ) %>%
    filter(decade != first(decade)) %>%
    mutate(
      across(
        all_of(cols),
        \(x) x - pick(paste0(cur_column(), "_1960"))[[1]],
        .names = "{.col}_diff"
      )
    )
}

# Function to plot, for both Seasons, ordinary and differenced values
plot_pdf <- \(x, x_diff, areas, var, col_lab = NULL) {
  # plot for both seasons
  season <- list("Summer", "Winter", NULL)
  season_lab <- c(unlist(season[1:2]), "Year")
  p_lst <- lapply(seq_along(season), \(i) {
    map_plot(x, areas, var, season[[i]], col_lab, season_lab[[i]])
  })
  # plot differences
  p_diff_lst <- lapply(seq_along(season), \(i) {
    map_plot(
      x_diff, areas, paste0(var, "_diff"), season[[i]],
      col_lab, season_lab[[i]],
      fill_fun = diff_fill_fun
    )
  })
  return(c(p_lst, p_diff_lst))
}

# Function to calc chi, chibar at quantile q (or closest q) at each  site
calc_chi <- \(data, var1, var2, chi_q) {
  stations <- unique(data$station_id)
  chi_95_df <- bind_rows(lapply(stations, \(x) {
    chi <- data %>%
      filter(station_id == x) |>
      dplyr::select(!!var1, !!var2) %>%
      texmex::chi()

    # whether to show chi or not, based on whether chibar upper extend crosses 1
    show_chi <- !prod(tail(chi$chibar[, 3]) < 1)

    loc <- which.min(abs(chi$quantile - chi_q))
    return(data.frame(
      "station_id" = x,
      "chi" = chi$chi[loc, 2, drop = TRUE],
      "chibar" = chi$chibar[loc, 2, drop = TRUE],
      "show_chi" = show_chi
    ))
  }))
  rownames(chi_95_df) <- NULL

  # join in area statistics
  chi_95_df %>%
    pivot_longer(c("chi", "chibar"), names_to = "var") %>%
    # always show chibar plot
    mutate(show_chi = ifelse(value == "chibar", TRUE, show_chi)) %>%
    left_join(
      distinct(data, station_id, station_name, lon, lat)
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas), remove = FALSE)
}

# Function to plot chi and chibar statistics on map
chi_map_plot <- \(
  chi_95_sf,
  areas,
  var = c("chi", "chibar"),
  # scales = seq(-0.1, 0.6, by = 0.1),
  # point_ranges = c(2, 6),
  spec_locs = NULL,
  locs_nudge_x = 0.1,
  locs_nudge_y = 0.1,
  rm_axis = TRUE
) {
  # lab <- ifelse(var == "chi", expression(chi(u)), expression(bar(chi)(u)))
  lab <- ifelse(
    var == "chi",
    expression(chi(0.95)),
    expression(bar(chi)(0.95))
  )

  plot_data <- filter(chi_95_sf, var == !!var)
  if (var == "chi") {
    plot_data <- plot_data %>%
      mutate(show_chi = factor(ifelse(show_chi == TRUE, "yes", "no")))
  }

  p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      # geom_sf_pattern(
      data = plot_data,
      aes(fill = value, size = value),
      # aes(fill = value, size = value, pattern = show_chi),
      pch = 21,
      stroke = 1
    )

  # add numbers for specific locations, if desired
  if (!is.null(spec_locs)) {
    # plot_data <- plot_data |>
    #   mutate(spec_name = ifelse(name %in% spec_locs, TRUE, FALSE))

    plot_data_num <- plot_data %>%
      filter(name %in% spec_locs) |>
      # arrange by inputted spec_locs
      slice(match(spec_locs, name)) %>%
      mutate(label_num = row_number())

    p <- p +
      geom_sf(
        data = plot_data_num,
        aes(fill = value, size = value),
        colour = "red",
        pch = 21,
        stroke = 1,
        show.legend = FALSE
      ) +
      # add numbers beside points
      geom_sf_text(
        data = plot_data_num,
        aes(label = label_num),
        nudge_x = locs_nudge_x,
        nudge_y = locs_nudge_y,
        # size = 5,
        size = 6,
        fontface = "bold",
        colour = "red",
        show.legend = FALSE
      )
  }

  p <- p +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # add colour scheme afterwards
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   breaks = scales_chi,
    #   labels = as.character(scales_chi),
    #   guide = "legend"
    # ) +
    # scale_size_continuous(
    #   range  = point_ranges,
    #   breaks = scales,
    #   labels = as.character(scales),
    #   guide  = "legend"
    # ) +
    # scale_pattern_manual(values = c("stripe", "none")) +
    # maximise plot within frame
    coord_sf(expand = FALSE) +
    labs(fill = lab, size = lab, x = "", y = "") +
    guides(fill = guide_legend(), size = guide_legend(), pattern = "none") +
    # CeCl::cecl_theme(legend.position = "right") +
    CeCl::cecl_theme() +
    theme(legend.key = element_blank())

  # remove axis text and ticks if required
  if (rm_axis == TRUE) {
    p <- p +
      theme(
        axis.text  = element_blank(),
        axis.ticks = element_blank()
      )
  }

  return(p)
}

# Function to join both chi and chibar plots together
join_chi_plots <- \(chi_95_sf, areas, fill_fun1, fill_fun2) {
  # plot chibar and chi
  chibar_p <- chi_map_plot(chi_95_sf, areas, "chibar", rm_axis = FALSE) +
    # scale_fill_gradientn(
    #   colours = rev(heat.colors(7)),
    #   # breaks = scales,
    #   # labels = as.character(scales),
    #   guide = "legend"
    # )
    diff_fill_fun()

  # chi_p <- chi_map_plot(chi_95_sf, "chi") +
  chi_p <- chi_map_plot(chi_95_sf, areas, "chi", rm_axis = FALSE) +
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   # breaks = scales,
    #   # labels = as.character(scales),
    #   guide = "legend"
    # ) +
    diff_fill_fun()

  chibar_p + chi_p
}


#### Metadata ####

q <- 0.99 # plot 99th quantile to explore extremes (max is quite noisy)
chi_q <- 0.95 # ???

decades <- seq(1960, 2010, by = 10) # complete decades in record only

# divergent colour scheme (for difference plots)
diff_fill_fun <- \(...) {
  scale_fill_gradient2(
    low = "blue3", high = "red3", na.value = "grey", guide = "legend", ...
  )
}

#### Load Data ####

# load ECAD data
data <- readr::read_csv("data/02_app/ecad_clean.csv.gz")

# load map of continental Spain to use as background in plots
areas <- mapSpain::esp_get_munic_siane(moveCAN = TRUE) |>
  filter(!ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla"))


#### Initial Calculations ####

# add decade and season labels
data <- data %>%
  mutate(
    decade = floor(year(date) / 10) * 10,
    season = case_when(
      month(date) %in% c(1:3, 10:12) ~ "Winter",
      TRUE ~ "Summer"
    )
  ) |>
  filter(decade %in% decades)

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# change to sf object
data_sf <- st_as_sf(
  data,
  coords = c("lon", "lat"),
  crs = st_crs(areas_ccaa)
)

readr::write_csv(
  data,
  "data/02_app/ecad_clean_full.csv.gz"
)


#### Means ####

# calculate means
data_mean <- data_sf |>
  group_by(station_id, station_name, decade, season) |>
  summarise(
    across(c(temp, rain, wind_speed, contains("drought")), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )


# calculate differences from first decade (1960s)
data_mean_diff <- calc_diff(
  data_mean,
  c("temp", "rain", "wind_speed", "drought_global", "drought_local"),
  first(decades)
)

# create PDFs of saved plots for both Seasons, as well as across the whole Year
pdf("plots/02_app/01_temp_means.pdf", width = 10, height = 8)
plot_pdf(data_mean, data_mean_diff, areas_ccaa, "temp", "°C")
dev.off()

pdf("plots/02_app/01_rain_means.pdf", width = 10, height = 8)
plot_pdf(data_mean, data_mean_diff, areas_ccaa, "rain", "mm")
dev.off()

pdf("plots/02_app/01_drought_global_means.pdf", width = 10, height = 8)
plot_pdf(data_mean, data_mean_diff, areas_ccaa, "drought_local", "Local Drought Index")
dev.off()

pdf("plots/02_app/01_drought_local_means.pdf", width = 10, height = 8)
plot_pdf(data_mean, data_mean_diff, areas_ccaa, "drought_global", "Global Drought Index")
dev.off()

pdf("plots/02_app/01_ws_means.pdf", width = 10, height = 8)
plot_pdf(data_mean, data_mean_diff, areas_ccaa, "wind_speed", "m/s")
dev.off()


#### Extremes ####

# calculate 99th quantile, to look at extremes
data_q <- data_sf |>
  group_by(station_id, station_name, decade, season) |>
  summarise(
    across(c(temp, rain, wind_speed, contains("drought")), \(x) quantile(x, q, na.rm = TRUE)),
    .groups = "drop"
  )

data_q_diff <- calc_diff(
  data_q,
  c("temp", "rain", "wind_speed", "drought_global", "drought_local"),
  first(decades)
)

# create PDFs of saved plots for both Seasons, as well as across the whole Year
pdf("plots/02_app/01a_temp_q.pdf", width = 10, height = 8)
plot_pdf(data_q, data_q_diff, areas_ccaa, "temp", "°C")
dev.off()

pdf("plots/02_app/01a_rain_q.pdf", width = 10, height = 8)
plot_pdf(data_q, data_q_diff, areas_ccaa, "rain", "mm")
dev.off()

pdf("plots/02_app/01a_drought_global_q.pdf", width = 10, height = 8)
plot_pdf(data_q, data_q_diff, areas_ccaa, "drought_local", "Local Drought Index")
dev.off()

pdf("plots/02_app/01a_drought_local_q.pdf", width = 10, height = 8)
plot_pdf(data_q, data_q_diff, areas_ccaa, "drought_global", "Global Drought Index")
dev.off()

pdf("plots/02_app/01a_ws_q.pdf", width = 10, height = 8)
plot_pdf(data_q, data_q_diff, areas_ccaa, "wind_speed", "m/s")
dev.off()


#### Chi ####

# Plot chi and chi-squared in Summer for one location
spec_loc <- "Valencia"
chi_spec <- data %>%
  filter(
    station_name == spec_loc,
    # decade == last(decades),
    season == "Summer"
  ) |>
  # dplyr::select(temp, rain) %>%
  dplyr::select(temp, drought_global) %>%
  # texmex::chi(qlim = c(0.1, 0.82)) # highest and lowest quantiles available
  texmex::chi()

(chi_plot_spec <- chi_spec %>%
  # use texmex plotting method for chi and chibar
  ggplot(main = c("ChiBar" = "", "Chi" = ""), plot. = FALSE) |>
  lapply(\(x) x + CeCl::cecl_theme()) %>%
  wrap_plots() +
  # `[[`(2) +
  # add centred title through patchwork
  patchwork::plot_annotation(
    title = paste0("Tail Dependence, ", spec_loc),
    theme = CeCl::cecl_theme()
  ))

# calculate chi, chibar at 95th quantile for all sites
# In Summer, for temperature and rain
chi_95_sf <- calc_chi(
  filter(data, season == "Winter"),
  var1 = "temp", var2 = "drought_global", chi_q = chi_q
)

# plot for every setting
join_chi_plots(chi_95_sf, areas_ccaa, diff_fill_fun, diff_fill_fun) +
  plot_annotation(title = "Temperature vs Global Drought Index, Winter")

chi_settings_df <- tidyr::crossing(
  "season" = c("Winter", "Summer", NA),
  # "combo"  = c("temp - rain", "rain - wind_speed"),
  "combo"  = c("temp - drought_global", "temp - drought_local", "rain - wind_speed"),
  "decade" = c(NA, decades)
) |>
  separate_wider_delim(combo, delim = " - ", names = c("var1", "var2")) |>
  arrange(season, var1, var2, desc(is.na(decade)), decade)
chi_settings_df

# TODO Functionalise
plot_list <- mclapply(seq_len(nrow(chi_settings_df)), \(i) {
  row <- chi_settings_df[i, ]
  season <- row[["season"]]
  var1 <- row[["var1"]]
  var2 <- row[["var2"]]
  decade <- row[["decade"]]

  # var1_lab <- ifelse(var1 == "rain", "Precipitation", "Temperature")
  # var2_lab <- ifelse(var2 == "rain", "Precipitation", "Wind Speed")

  var_lab <- \(x) {
    case_when(
      x == "rain" ~ "Precipitation",
      x == "temp" ~ "Temperature",
      x == "wind_speed" ~ "Wind Speed",
      x == "drought_local" ~ "Local Drought Index",
      TRUE ~ "Global Drought Index"
    )
  }
  var1_lab <- var_lab(var1)
  var2_lab <- var_lab(var2)

  data_spec <- data
  season_lab <- ", Year"
  if (!is.na(season)) {
    data_spec <- filter(data_spec, season == !!season)
    season_lab <- paste(",", season)
  }

  decade_lab <- ""
  if (!is.na(decade)) {
    data_spec <- filter(data_spec, decade == !!decade)
    decade_lab <- paste0(", ", decade, "s")
  }

  # finally, check for locations with all missing data for one variable, and
  # ignore
  na_locs <- data_spec |>
    group_by(station_id, station_name) |>
    summarise(
      n = n(),
      # n_na_var1 = sum(is.na(!!var1)),
      n_na_var1 = sum(is.na(.data[[var1]])),
      n_na_var2 = sum(is.na(.data[[var2]])),
      .groups = "drop"
    ) |>
    # where there are all NAs for a variable
    filter(n_na_var1 == n | n_na_var2 == n)

  data_spec <- anti_join(data_spec, na_locs)

  chi_95_sf_spec <- calc_chi(
    data_spec,
    var1 = var1, var2 = var2, chi_q = chi_q
  ) |>
    mutate(var1 = !!var1, var2 = !!var2, season = !!season, decade = !!decade)

  p_spec <- join_chi_plots(chi_95_sf, areas_ccaa, diff_fill_fun, diff_fill_fun) +
    # plot_annotation(title = "Temperature vs Precipitation, Winter")
    plot_annotation(title = paste0(
      var1_lab, " vs ", var2_lab, season_lab, decade_lab
    ))

  return(list("data" = chi_95_sf_spec, "plot" = p_spec))
})
names(plot_list) <- with(
  chi_settings_df,
  paste(
    var1,
    var2,
    ifelse(is.na(season), "Year", season),
    ifelse(is.na(decade), "All", decade),
    sep = " - "
  )
)

# save plots
# pdf("plots/02_app/01b_chi_maps_temp_rain.pdf", width = 10, height = 8)
# lapply(plot_list[grepl("temp", names(plot_list))], `[[`, "plot")
# dev.off()

pdf("plots/02_app/01b_chi_maps_temp_drought_local.pdf", width = 10, height = 8)
lapply(plot_list[grepl("drought_local", names(plot_list))], `[[`, "plot")
dev.off()

pdf("plots/02_app/01b_chi_maps_temp_drought_global.pdf", width = 10, height = 8)
lapply(plot_list[grepl("drought_global", names(plot_list))], `[[`, "plot")
dev.off()

pdf("plots/02_app/01b_chi_maps_rain_ws.pdf", width = 10, height = 8)
lapply(plot_list[grepl("wind_speed", names(plot_list))], `[[`, "plot")
dev.off()

# Calculate differences and plot
plot_list_1960 <- plot_list[grepl("1960", names(plot_list))]

df_1960 <- lapply(plot_list_1960, `[[`, "data") |>
  bind_rows() |>
  rename(value_1960 = value) |>
  mutate(season = ifelse(is.na(season), "Year", season)) |>
  select(-decade) |>
  st_drop_geometry()

plot_list_diff <- plot_list[!grepl("All|1960", names(plot_list), fixed = FALSE)]

plot_list_diff <- lapply(plot_list_diff, \(x) {
  chi_95_sf_diff <- x$data |>
    mutate(season = ifelse(is.na(season), "Year", season)) |>
    left_join(df_1960) |>
    mutate(value = value - value_1960) |>
    select(-value_1960)

  p_spec <- join_chi_plots(
    chi_95_sf_diff, areas_ccaa, diff_fill_fun, diff_fill_fun
  ) +
    plot_annotation(title = paste("Diff,", x$plot$patches$annotation$title))

  return(list("data" = chi_95_sf_diff, "plot" = p_spec))
})

plot_list_diff$`temp - drought_global - Summer - 1990`$plot

# save plots
# pdf("plots/02_app/01c_chi_maps_diff_temp_rain.pdf", width = 10, height = 8)
# lapply(plot_list_diff[grepl("temp", names(plot_list_diff))], `[[`, "plot")
# dev.off()

pdf("plots/02_app/01c_chi_maps_diff_temp_drought_local.pdf", width = 10, height = 8)
lapply(plot_list_diff[grepl("drought_local", names(plot_list_diff))], `[[`, "plot")
dev.off()

pdf("plots/02_app/01c_chi_maps_diff_temp_drought_global.pdf", width = 10, height = 8)
lapply(plot_list_diff[grepl("drought_global", names(plot_list_diff))], `[[`, "plot")
dev.off()

pdf("plots/02_app/01c_chi_maps_diff_rain_ws.pdf", width = 10, height = 8)
lapply(plot_list_diff[grepl("wind_speed", names(plot_list_diff))], `[[`, "plot")
dev.off()


#### Extremal autocorrelation ####

# Empirical lagged extremogram:
# P(X_{t+h} > u | X_t > u)
lag_extremogram <- \(x, q = 0.95, maxlag = 20) {
  u <- quantile(x, q, na.rm = TRUE)
  I <- x > u

  tibble(
    lag = 1:maxlag,
    chi = sapply(1:maxlag, function(h) {
      idx <- seq_len(length(I) - h)

      denom <- sum(I[idx], na.rm = TRUE)

      if (denom == 0) {
        return(NA_real_)
      }

      sum(I[idx] & I[idx + h], na.rm = TRUE) / denom
    })
  )
}


# Rolling version for one station and one season
roll_ext_station <- \(data_spec,
  var = "temp",
  q = 0.95,
  years = 15,
  maxlag = 20) {
  data_spec <- data_spec |>
    arrange(date)

  rows_per_year <- data_spec |>
    mutate(year = lubridate::year(date)) |>
    count(year) |>
    summarise(n = median(n)) |>
    pull(n)

  window <- round(years * rows_per_year)
  step <- round(rows_per_year)

  starts <- seq(
    1,
    nrow(data_spec) - window + 1,
    by = step
  )

  purrr::map_dfr(
    starts,
    function(s) {
      idx <- s:(s + window - 1)

      lag_extremogram(
        x = data_spec[[var]][idx],
        q = q,
        maxlag = maxlag
      ) |>
        mutate(
          start = data_spec$date[min(idx)],
          end = data_spec$date[max(idx)],
          mid = start + (end - start) / 2
        )
    }
  )
}

pdf(
  "plots/02_app/02_rolling_extremogram_temp.pdf",
  width = 12,
  height = 8
)

for (st in sort(unique(data$station_name))) {
  roll_df <- purrr::map_dfr(
    c("Winter", "Summer"),
    \(ss) {
      d <- data |>
        filter(
          station_name == st,
          season == ss
        )

      if (nrow(d) == 0) {
        return(tibble())
      }

      roll_ext_station(
        data_spec = d,
        var = "temp",
        q = 0.95,
        years = 15,
        # maxlag = 20
        maxlag = 10
      ) |>
        mutate(season = ss)
    }
  )

  p <- ggplot(
    roll_df,
    aes(lag, mid, fill = chi)
  ) +
    geom_tile(height = 365, width = 1) +
    # geom_raster() +
    facet_wrap(~season, ncol = 2) +
    scale_y_date(
      date_breaks = "10 years",
      date_labels = "%Y"
    ) +
    scale_fill_viridis_c(
      limits = c(0, 1),
      na.value = "grey80"
    ) +
    labs(
      title = st,
      x = "Lag (weeks)",
      y = "Window midpoint",
      fill = expression(
        P(X[t + h] > u ~ "|" ~ X[t] > u)
      )
    ) +
    cecl_theme(nejm_pal = FALSE) +
    theme(legend.text = element_text(angle = 45, hjust = 1))

  print(p)
}

dev.off()

# TODO functionalise
pdf(
  "plots/02_app/02_rolling_extremogram_rain.pdf",
  width = 12,
  height = 8
)

for (st in sort(unique(data$station_name))) {
  roll_df <- purrr::map_dfr(
    c("Winter", "Summer"),
    \(ss) {
      d <- data |>
        filter(
          station_name == st,
          season == ss
        )

      if (nrow(d) == 0) {
        return(tibble())
      }

      roll_ext_station(
        data_spec = d,
        var = "rain",
        q = 0.95,
        years = 15,
        # maxlag = 20
        maxlag = 10
      ) |>
        mutate(season = ss)
    }
  )

  p <- ggplot(
    roll_df,
    aes(lag, mid, fill = chi)
  ) +
    geom_tile(height = 365, width = 1) +
    # geom_raster() +
    facet_wrap(~season, ncol = 2) +
    scale_y_date(
      date_breaks = "10 years",
      date_labels = "%Y"
    ) +
    scale_fill_viridis_c(
      limits = c(0, 1),
      na.value = "grey80"
    ) +
    labs(
      title = st,
      x = "Lag (weeks)",
      y = "Window midpoint",
      fill = expression(
        P(X[t + h] > u ~ "|" ~ X[t] > u)
      )
    ) +
    cecl_theme(nejm_pal = FALSE) +
    theme(legend.text = element_text(angle = 45, hjust = 1))

  print(p)
}

dev.off()

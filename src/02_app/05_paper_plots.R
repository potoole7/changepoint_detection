#### Paper plots ####

# Plots for Application section of paper

# TODO What about matrices for individual variables? Right now we're only
# looking at aggregated matrices!!! Would need to do entire analysis again ...

# TODO Also need to do for DQU = 0.85 ...

# TODOs for intro plot
# Try make elev_df smaller so it's easier to plot (done)
# Fix map so that weird very low down location isn't present!!! (done)
# Fill in location of two spec_locs on map! (done)
# TODO Can we actually use 1950s data for some locations???
# TODO Do I need to label SPI as negative SPI? Or can this be introduced in
# the paper?
# Add elevation to plot of country (done)
# Clip significant whitespace away from plot (done)

# TODOs for CE plots
# TODO Should we use the same Laplace sample for each clustering, like we did
# with the changepoint algorithm?
# TODO For Winter, should I use all data pre and post 1998 to cluster?
# Will obviously be different number of years for both
# TODO Some k values in elbow plots suggest different values than k = 3! plot
# for other k values, even just to show Christian ...
# Add elevations (done)
# Ensure clusterings don't change colours!! (done)

#### Libs ####

devtools::load_all("../CeCl")
library(grid)
library(lubridate)
library(stringr)
library(RColorBrewer)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(purrr)
library(parallel)
library(evgam)
library(ggridges)
library(sf)
library(patchwork)

# source custom functions
source("src/00_functions.R")


#### Metadata ####

# variables of interest
dep_var <- "drought_local_rev"
temp_var <- "temp_max"

seasons <- c("Winter", "Spring", "Summer", "Autumn")

# chi_q <- 0.95
spec_years <- 1960:2020

# divergent colour scheme (for difference plots)
# TODO Needed?
diff_fill_fun <- \(...) {
  scale_fill_gradient2(
    low = "blue3", high = "red3", na.value = "grey", guide = "legend", ...
  )
}


# Identified changepoint years
change_year_winter <- 1998
change_year_summer <- 1990

n_perm <- 200
# dqu <- 0.8
dqu <- 0.85


#### Functions ####


#### Load Data ####

# data
data <- readr::read_csv(
  "data/02_app/ecad_clean.csv.gz"
) |>
  # mutate(decade = factor(floor(year(date) / 10) * 10, levels = decades)) |>
  mutate(
    year = as.numeric(format(as.Date(date), "%Y"))
  ) |>
  rename(name = station_name) |>
  identity()


if (dep_var == "rain") {
  data <- data |>
    filter(rain > 0)
}

# if specified, use maximum temperature rather than 90th quantile
if (temp_var == "temp_max") {
  data <- data |>
    mutate(
      temp_max = ifelse(is.infinite(temp_max), NA, temp_max),
      temp     = temp_max
    ) |>
    filter(!is.na(temp))
  temp_var <- "temp"
}

# reverse drought_local variable to give positive alpha values, if desired
if (dep_var == "drought_local_rev") {
  data <- data |>
    mutate(drought_local = -drought_local)
  dep_var <- c("drought_local")
}

# check percentage of dates
station_count <- n_distinct(data$name)

date_coverage <- bind_rows(lapply(seasons, \(s) {
  season_data <- data |>
    filter(.data$season == s)

  all_season_dates <- season_data |>
    distinct(.data$date) |>
    pull(.data$date)

  valid_dates <- season_data |>
    # Add filters for missing measurements here if required
    distinct(.data$date, .data$name) |>
    count(.data$date, name = "n_stations") |>
    filter(.data$n_stations == station_count) |>
    pull(.data$date)

  tibble(
    season        = s,
    n_valid_dates = length(valid_dates),
    n_dates       = length(all_season_dates),
    perc          = n_valid_dates / n_dates
  )
}))

date_coverage # 94% - 80%: Fine!

# only keep dates available at every station (by season)
valid_dates <- data %>%
  distinct(season, date, name) %>%
  group_by(season) %>%
  mutate(n_stations_in_season = n_distinct(name)) %>%
  group_by(season, date) %>%
  filter(n_distinct(name) == first(n_stations_in_season)) %>%
  distinct(season, date)

data <- data %>%
  semi_join(valid_dates, by = c("season", "date"))


# # marginal fits
# marg_season <- readRDS(file)
marg_season <- readRDS("data/02_app/marg_season_roll_emp.rds")
# if (!is.nkull(names(marg_season))) {
#   seasons <- names(marg_season)
# }

# load map of continental Spain to use as background in plots
areas <- read_sf("data/02_app/spain_shapefile.geojson") |>
  filter(
    !ine.ccaa.name %in% c("Canarias", "Balears, Illes", "Ceuta", "Melilla")
  )

# simplify areas into autonomous communities/provinces
areas_ccaa <- areas %>%
  group_by(ine.ccaa.name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ensure areas_ccaa has the same bounding box as elev_df
areas_proj <- st_transform(areas_ccaa, crs = 25830)

# extract point location of each station for plotting on map
# pts <- data %>%
#   distinct(name, lon, lat) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = st_crs(areas_proj))
pts <- data %>%
  distinct(name, lon, lat) %>%
  st_as_sf(
    coords = c("lon", "lat"),
    crs = st_crs(areas_ccaa) # ETRS89, so handles lon/lat in degrees
  ) %>%
  st_transform(st_crs(areas_proj)) # EPSG:25830: UTM coordinates in metres

station_names <- unique(data$name)

# TODO Try make smaller, lots of rows, don't need it to be that granular!
# elev_df <- readr::read_csv("data/spain_elev_binned.csv")
elev_df <- readr::read_csv("data/spain_elev_binned_reduced.csv")
# ensure bounding box is same as for areas
elev_df <- elev_df %>%
  filter(
    x >= st_bbox(areas_proj)$xmin,
    x <= st_bbox(areas_proj)$xmax,
    y >= st_bbox(areas_proj)$ymin,
    y <= st_bbox(areas_proj)$ymax
  )


#### Introduction plot (scatter plots) ####


# find names and location of (optionally, specific) stations
# plot_locs <- \(pts, areas, add_names = FALSE, spec_locs = NULL, elev_df = NULL, legend.position = "bottom") {
#   # browser()
#   pts_plt <- pts
#   # filter to specific locations, if desired
#   if (!is.null(spec_locs)) {
#     pts_plt <- filter(pts, name %in% spec_locs)
#   }
#   # initialise plot
#   p <- pts_plt |>
#     ggplot()
#
#   # optionally, add elevation raster to background of plot
#   if (!is.null(elev_df)) {
#     p <- p +
#       geom_tile(
#         # ensure order is correct
#         data = elev_df |>
#           mutate(
#             elev_bin = factor(elev_bin, levels = unique(.data[["elev_bin"]]))
#           ),
#         aes(x = x, y = y, fill = elev_bin),
#         width = diff(range(elev_df$x)) / length(unique(elev_df$x)),
#         height = diff(range(elev_df$y)) / length(unique(elev_df$y)),
#       ) +
#       cecl_theme(nejm_pal = FALSE, legend.position = legend.position) +
#       labs(x = "", y = "") +
#       scale_fill_manual(
#         values = c(
#           "<200"      = "#2E8B57",
#           "200–600" = "#9ACD66",
#           "600–1000" = "#F7E56B",
#           "1000–1500" = "#F29E38",
#           "1500–2000" = "#C46A3A",
#           ">2000"     = "#7A3B1E"
#         ),
#         drop = FALSE,
#         name = "Elevation (m)"
#       )
#
#     # right-justify facet labels, if desired
#     if (legend.position %in% c("left", "right")) {
#       p <- p +
#         theme(
#           legend.text = element_text(hjust = 1)
#         )
#     }
#   } else {
#     p <- p +
#       cecl_theme(legend.position = legend.position)
#   }
#
#   # plot station locations
#   p <- p +
#     geom_sf(data = areas, colour = "black", fill = NA) + # region boundaries
#     # TODO Replace with hollow dots
#     # geom_sf(size = 5, alpha = 0.9)
#     geom_sf(
#       size = 5,
#       shape = 21,
#       fill = NA,
#       colour = "black",
#       # stroke = 1,
#       # stroke = 1.2,
#       stroke = 1.5,
#       alpha = 0.9
#     )
#
#   # optionally, add station names to map
#   if (add_names == TRUE) {
#     p <- p +
#       ggrepel::geom_text_repel(
#       aes(label = name, geometry = geometry),
#       stat = "sf_coordinates"
#     )
#   }
#
#   return(p)
# }
plot_locs <- \(
  pts,
  areas,
  add_names = FALSE,
  spec_locs = NULL,
  highlight_locs = NULL,
  elev_df = NULL,
  legend.position = "bottom"
) {
  pts_plt <- pts

  if (!is.null(spec_locs)) {
    pts_plt <- pts |>
      dplyr::filter(name %in% spec_locs)
  }

  p <- ggplot2::ggplot()

  if (!is.null(elev_df)) {
    p <- p +
      ggplot2::geom_tile(
        data = elev_df |>
          dplyr::mutate(
            elev_bin = factor(
              elev_bin,
              levels = unique(.data[["elev_bin"]])
            )
          ),
        ggplot2::aes(x = x, y = y, fill = elev_bin),
        width = diff(range(elev_df$x)) /
          length(unique(elev_df$x)),
        height = diff(range(elev_df$y)) /
          length(unique(elev_df$y))
      ) +
      cecl_theme(
        nejm_pal = FALSE,
        legend.position = legend.position
      ) +
      ggplot2::labs(x = "", y = "") +
      ggplot2::scale_fill_manual(
        values = c(
          "<200"      = "#2E8B57",
          "200–600"   = "#9ACD66",
          "600–1000"  = "#F7E56B",
          "1000–1500" = "#F29E38",
          "1500–2000" = "#C46A3A",
          ">2000"     = "#7A3B1E"
        ),
        drop = FALSE,
        name = "Elevation (m)"
      )

    if (legend.position %in% c("left", "right")) {
      p <- p +
        ggplot2::theme(
          legend.text = ggplot2::element_text(hjust = 1)
        )
    }
  } else {
    p <- p +
      cecl_theme(legend.position = legend.position)
  }

  # Region boundaries
  p <- p +
    ggplot2::geom_sf(
      data = areas,
      colour = "black",
      fill = NA
    )

  # All displayed stations
  p <- p +
    ggplot2::geom_sf(
      data = pts_plt,
      size = 5,
      shape = 21,
      fill = NA,
      colour = "black",
      stroke = 1.5,
      alpha = 0.9
    )

  # Optionally highlight selected stations
  if (!is.null(highlight_locs)) {
    highlight_pts <- pts |>
      dplyr::filter(name %in% highlight_locs) |>
      dplyr::mutate(
        plot_name = stringr::str_remove(name, " Aeropuerto$")
      )

    p <- p +
      ggplot2::geom_sf(
        data = highlight_pts,
        shape = 21,
        size = 7,
        # fill = "#F7E56B",
        fill = "black",
        colour = "black",
        stroke = 1.5,
        alpha = 1
      ) +
      # ggrepel::geom_text_repel(
      #   data = highlight_pts,
      #   ggplot2::aes(
      #     label = plot_name,
      #     geometry = geometry
      #   ),
      #   stat = "sf_coordinates",
      #   fontface = "bold",
      #   size = 4,
      #   box.padding = 0.5,
      #   point.padding = 0.5,
      #   min.segment.length = 0
      # ) +
      NULL
  }

  # Optionally label all displayed stations
  if (isTRUE(add_names)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = pts_plt,
        ggplot2::aes(
          label = name,
          geometry = geometry
        ),
        stat = "sf_coordinates"
      )
  }

  p
}

# find (pearson) correlation between temp and drought_local at each location
# pick locations with highest and lowest correlation
spec_locs <- data |>
  filter(year %in% spec_years, !is.na(temp), !is.na(drought_local)) |>
  group_by(name) |>
  # group_by(name, season) |>
  summarise(cor = cor(temp, drought_local), .groups = "drop") |>
  arrange(desc(cor)) |>
  slice(c(1, n())) |>
  pull(name)

# (p_loc_names <- plot_locs(pts, areas_proj))
# ggsave(plot = p_loc_names, "plots/02_app/p_loc_names.png", width = 14, height = 10)
p_loc_names_elev <- plot_locs(
  pts,
  areas_proj,
  elev_df = elev_df,
  highlight_locs = spec_locs,
  legend.position = "right"
  # legend.position = "bottom"
)
p_loc_names_elev
ggsave(plot = p_loc_names_elev, "plots/02_app/p_loc_names_elev.png", width = 14, height = 10)

# pick two locations in South and North of Spain, and compare
# spec_locs <- c("Malaga", "Bilbao Aeropuerto")
# spec_locs <- c("Malaga", "Santander")

data_laplace <- bind_rows(lapply(seq_along(marg_season), \(i) {
  x <- marg_season[[i]]
  bind_rows(lapply(seq_along(x$transformed), \(j) {
    y <- as.data.frame(x$transformed[[j]])
    y$name <- names(x$transformed)[[j]]
    y$season <- names(marg_season)[[i]]
    # TODO Need to add date?

    y
  }))
}))

# TODO Replace with Laplace transformed data!
# data_intro <- data |>
#   filter(name %in% spec_locs, year %in% spec_years) |>
#   mutate(name = ifelse(
#     grepl("Aeropuerto", name),
#     str_remove(name, " Aeropuerto"),
#     name
#   ))

data_intro <- data_laplace |>
  filter(name %in% spec_locs) |>
  mutate(name = ifelse(
    grepl("Aeropuerto", name),
    str_remove(name, " Aeropuerto"),
    name
  ))


# first, plot scatter plots of variables against eachother
p_scatter <- data_intro |>
  ggplot(aes(x = drought_local, y = temp, colour = season)) +
  geom_point() +
  facet_wrap(~name) +
  labs(y = "Temperature", x = "SPI", colour = "Season") +
  cecl_theme(legend.position = "right") +
  # cecl_theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 6)))
p_scatter

# now, plot time series for each
# data_intro |>
#   pivot_longer(c(temp, drought_local), names_to = "var") |>
#   mutate(ind = paste0(name, " - ", var)) |>
#   ggplot(aes(x = season_year, y = value, colour = season)) +
#   geom_point() +
#   facet_wrap(~ind) +
#   cecl_theme()
#
# data_intro |>
#   pivot_longer(c(temp, drought_local), names_to = "var") |>
#   mutate(ind = paste0(name, " - ", var)) |>
#   group_by(season, ind) |>
#   filter(value >= quantile(value, 0.95, na.rm = TRUE)) |>

# p_intro <- p_loc_names_elev + p_scatter
p_intro <- wrap_plots(list(p_loc_names_elev, p_scatter))

# ggsave(plot = p_intro, "latex/plots/intro.png", width = 12, height = 8)
ggsave(plot = p_intro, "latex/plots/intro_laplace.png", width = 12, height = 8)


#### Screening plot ####

screen_res_df <- readr::read_csv(paste0("data/02_app/screen_res_dqu_", dqu, ".csv.gz"))

screen_res_df_plt <- screen_res_df |>
  arrange(
    season,
    n_years_per_block,
    change_after_year
  ) |>
  group_by(
    season,
    n_years_per_block
  ) |>
  mutate(
    # Separate successful line segments at failed candidates.
    success_run = cumsum(!success),

    # Identify Frobenius local peaks.
    local_peak_frob = (
      success &
        lag(success, default = FALSE) &
        lead(success, default = FALSE) &
        frob > lag(frob) &
        frob >= lead(frob)
    ),
    local_peak_frob = replace_na(
      local_peak_frob,
      FALSE
    )
  ) |>
  ungroup() |>
  mutate(
    # setting = paste0(
    #   season,
    #   ", ",
    #   n_years_per_block,
    #   " years"
    # )
    setting = season
  )

# Preserve the desired facet ordering.
# setting_levels <- screen_res_df_plt |>
#   distinct(
#     season,
#     n_years_per_block,
#     setting
#   ) |>
#   arrange(
#     season,
#     n_years_per_block
#   ) |>
#   pull(setting) |>
#   unique()
setting_levels <- seasons

screen_res_df_plt <- screen_res_df_plt |>
  mutate(
    setting = factor(
      setting,
      levels = setting_levels
    )
  )

# plot norm peaks
year_breaks <- seq(
  floor(
    min(
      screen_res_df_plt$change_after_year,
      na.rm = TRUE
    ) / 2
  ) * 2,
  ceiling(
    max(
      screen_res_df_plt$change_after_year,
      na.rm = TRUE
    ) / 2
  ) * 2,
  by = 2
  # by = 1
)

# Convert all norms to long format
screen_res_long <- screen_res_df_plt |>
  pivot_longer(
    cols = c(
      frob,
      inf,
      spec
    ),
    names_to = "norm",
    values_to = "value"
  ) |>
  mutate(
    norm = recode(
      norm,
      frob = "Frobenius",
      inf = "Infinity",
      spec = "Spectral"
    ),
    norm = factor(
      norm,
      levels = c(
        "Frobenius",
        "Infinity",
        "Spectral"
      )
    )
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    change_after_year
  ) |>
  group_by(
    season,
    n_years_per_block,
    norm
  ) |>
  mutate(
    # Peaks are calculated separately for each norm.
    local_peak = (
      success &
        lag(success, default = FALSE) &
        lead(success, default = FALSE) &
        value > lag(value) &
        value >= lead(value)
    ),
    local_peak = replace_na(
      local_peak,
      FALSE
    )
  ) |>
  ungroup()


# Peak summaries for every norm
candidate_peaks_all <- screen_res_long |>
  filter(
    success,
    local_peak
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    desc(value)
  )

top_local_peaks_all <- candidate_peaks_all |>
  group_by(
    season,
    n_years_per_block,
    norm
  ) |>
  slice_max(
    order_by = value,
    n = 3,
    with_ties = FALSE
  ) |>
  arrange(
    season,
    n_years_per_block,
    norm,
    desc(value)
  ) |>
  ungroup()

# # Combined plot for all norms
# plot_all_norms <- \(df, spec_year = NULL) {
#   df_plot <- df
#   if (!is.null(spec_year)) {
#     df_plot <- df_plot |>
#       filter(n_years_per_block == spec_year)
#   }
#   p <- df_plot |>
#     ggplot(
#       aes(
#         x = change_after_year,
#         y = value,
#         colour = norm
#       )
#     ) +
#     geom_line(
#       data = \(x) filter(x, success),
#       aes(
#         group = interaction(
#           setting,
#           norm,
#           success_run
#         )
#       ),
#       show.legend = FALSE
#     ) +
#     geom_point(
#       data = \(x) filter(x, success)
#     ) +
#     geom_point(
#       data = \(x) filter(x, local_peak),
#       colour = "red",
#       size = 3
#     ) +
#     facet_wrap(
#       ~setting,
#       scales = "free_y"
#     ) +
#     scale_x_continuous(
#       breaks = year_breaks
#     ) +
#     scale_colour_brewer(
#       palette = "Dark2"
#     ) +
#     labs(
#       x = "Season Year",
#       y = expression(D),
#       colour = "Norm",
#     ) +
#     cecl_theme() +
#     theme(
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1
#       ),
#       legend.position = "bottom"
#     ) +
#     guides(colour = guide_legend(override.aes = list(size = 6)))
# }

# Combined plot for all norms
df <- screen_res_long
spec_year <- 25
plot_all_norms <- \(df, spec_year = NULL) {
  df_plot <- df
  if (!is.null(spec_year)) {
    df_plot <- df_plot |>
      filter(n_years_per_block == spec_year)
  }
  p <- df_plot |>
    group_by(norm, n_years_per_block) |>
    # scale between 0 and 1
    # mutate(value = scale(value, center = TRUE, scale = TRUE)) |>
    mutate(value = boot::inv.logit(scale(value))) |>
    filter(norm != "Spectral") |>
    ggplot(
      aes(
        x = change_after_year,
        y = value,
        # colour = norm
        colour = season
      )
    ) +
    geom_line(
      data = \(x) filter(x, success),
      aes(
        group = interaction(
          setting,
          norm,
          success_run
        )
      ),
      show.legend = FALSE
    ) +
    geom_point(
      data = \(x) filter(x, success)
    ) +
    geom_point(
      data = \(x) filter(x, local_peak),
      # colour = "red",
      colour = "black",
      shape = 4,
      size = 5
    ) +
    # facet_wrap(
    #   ~setting,
    #   scales = "free_y"
    # ) +
    facet_wrap(
      ~norm,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = year_breaks
    ) +
    scale_colour_brewer(
      palette = "Dark2"
    ) +
    labs(
      x = "Season Year",
      y = expression(D),
      # colour = "Norm",
      colour = "Season",
    ) +
    cecl_theme() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 6)))
}


(p_all_norms_25 <- plot_all_norms(screen_res_long, spec_year = 25))

# ggsave(paste0("plots/02_app/p_all_norms_25_dqu_", dqu, ".png"), p_all_norms_25, width = 12, height = 8)
# ggsave("latex/plots/screen.png", p_all_norms_25, width = 12, height = 8)
ggsave(paste0("latex/plots/screen_dqu_", dqu, ".png"), p_all_norms_25, width = 12, height = 8)


#### Changepoint plot ####

permutation_scan_results <- readRDS(
  paste0(
    "data/02_app/",
    "spain_perm_test_",
    "n_perm_",
    # 200,
    n_perm,
    "_dqu_",
    dqu,
    ".rds"
  )
)

# preprocess every season
permutation_scan_tidy <- lapply(
  permutation_scan_results,
  preprocess_permutation_scan
)

permutation_summary_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$summary
  )
)

permutation_values_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$permutations
  )
)

permutation_dependence_df <- bind_rows(
  lapply(
    permutation_scan_tidy,
    \(x) x$dependence
  )
)

# check failures
permutation_summary_df |>
  distinct(
    season,
    change_after_year,
    success,
    n_attempted,
    n_successful,
    n_failed,
    failure_rate,
    error_stage,
    error_message
  ) |>
  arrange(
    season,
    change_after_year
  )

# create plots for every season
permutation_scan_plots <- lapply(
  permutation_scan_tidy,
  plot_permutation_scan,
  plot_ce_parameters = TRUE,
  variable_names = c(
    X1 = "Maximum temperature",
    X2 = "Drought"
  )
)

# extract data for plotting for each season so we can facet properly
data_p <- bind_rows(lapply(names(permutation_scan_plots), \(x) {
  permutation_scan_plots[[x]]$p_value_profile$data
}))

year_breaks <- sort(
    unique(data_p$change_after_year)
)

# plot with faceted
p_change <- data_p |>
  filter(norm != "Spectral") |>
  mutate(season = factor(season, levels = seasons)) |>
  ggplot(
    aes(
      x = change_after_year,
      y = p_value,
      colour = norm
    )
  ) +
  geom_hline(
    yintercept = 0.05,
    colour = "grey40",
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0.10,
    colour = "gold4",
    linetype = "dashed"
  ) +
  geom_line(
    data = \(x) filter(x, success),
    aes(group = norm),
    show.legend = FALSE
  ) +
  geom_point(
    data = \(x) filter(x, success),
    size = 2
  ) +
  facet_wrap(~ season) +
  scale_x_continuous(
    breaks = year_breaks
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    x = "Season Year",
    y = "p-value",
    colour = "Norm"
  ) +
  cecl_theme() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 6)))
p_change

# ggsave("latex/plots/change.png", p_change, width = 12, height = 8)
ggsave(paste0("latex/plots/change_dqu_", dqu, ".png"), p_change, width = 12, height = 8)


#### Cluster maps before/after changepoint ####

# load marginal models
marg_season <- readRDS("data/02_app/marg_season_roll_emp.rds")

# join dates to transformed data
data_laplace_season <- lapply(marg_season, \(x) {
  bind_rows(Map(
    \(transformed, dates, station) {
      stopifnot(nrow(transformed) == length(dates))

      as.data.frame(transformed) |>
        mutate(
          date = as.Date(dates),
          name = station,
          .before = 1
        )
    },
    x$transformed,
    x$dates,
    names(x$transformed)
  ))
})

# add season_date and join data
data_laplace <- bind_rows(lapply(
  names(data_laplace_season),
  \(s) {
    data_laplace_season[[s]] |>
      mutate(
        season = s,
        season_year = if_else(
          season == "Winter" & month(date) == 12L,
          year(date) + 1L,
          year(date)
        )
      )
  }
))

# Split by season and changepoints (in the case of Summer and Winter)
marg_laplace_cp <- data_laplace |>
  mutate(
    group = case_when(
      season %in% c("Spring", "Autumn")        ~ 1,
      season == "Winter" & season_year <= 1998 ~ 1,
      season == "Summer" & season_year <= 1990 ~ 1,
      TRUE                                     ~ 2
    )
  ) |>
  mutate(
    season = factor(season, levels = seasons),
    group  = factor(group, levels = c(1, 2)),
    ind = paste0(season, " - ", group),
    ind = case_when(
      ind == "Winter - 1" ~ "Winter - 1960-1998",
      ind == "Winter - 2" ~ "Winter - 1999-2020",
      ind == "Summer - 1" ~ "Summer - 1960-1990",
      ind == "Summer - 2" ~ "Summer - 1991-2020",
      TRUE                ~ season
    ),
    ind = factor(
      ind,
      levels = c(
        "Winter - 1960-1998", "Winter - 1999-2020",
        "Spring",
        "Summer - 1960-1990", "Summer - 1991-2020",
        "Autumn"
      )
    )
  ) |>
  group_split(ind, .keep = TRUE) |>
  lapply(\(x) {
    ret <- x |>
      select(name, temp, drought_local) |>
      as_cecl_marg()

    ret$date <- x$date
    ret$season_year <- x$season_year
    ret$season <- unique(x$season)
    ret$group <- unique(x$group)
    ret$ind <- unique(x$ind)

    ret
  })

# label appropriately
names(marg_laplace_cp) <- vapply(marg_laplace_cp, \(x) {
  # paste0(x$season, " - ", x$group)
  paste0(x$ind)
}, character(1))


# fit CE model for each
ce_cp <- lapply(marg_laplace_cp, \(x) {
  cecl_dep(
    x,
    vars      = c("temp", dep_var),
    cond_vars = c("temp", dep_var),
    cond_prob = dqu,
    nruns     = 3,
    ncores    = getOption("mc_cores", 2L)
  )
})


# calculate distances
set.seed(123)
dist_cp <- lapply(seq_along(ce_cp), \(i) {
  cecl_dist(
    # dep_var_season[[i]],
    # marg_season[[name_df$season[[i]]]],
    ce_cp[[i]],
    marg_laplace_cp[[i]],
    ncores = getOption("mc_cores", 2L)
  )
})
names(dist_cp) <- names(ce_cp)

# # function to label years from group and season
# year_from_season_group <- \(season, group) {
#  return(case_when(
#     # season %in% c("Spring", "Autumn") ~ NA,
#     season == "Summer" & group == 1   ~ "1960-1990",
#     season == "Summer" & group == 2   ~ "1991-2020",
#     season == "Winter" & group == 1   ~ "1960-1998",
#     season == "Winter" & group == 2   ~ "1999-2020",
#     TRUE                              ~ NA
#   ))
# }

# find k via Elbow plot
twgss_vals <- bind_rows(lapply(seq_along(dist_cp), \(i) {

  season <- as.character(marg_laplace_cp[[i]]$season)
  # group <- as.numeric(marg_laplace_cp[[i]]$group)
  # years <- year_from_season_group(season, group)
  years <- ifelse(
    season == names(dist_cp)[[i]],
    NA,
    str_split(names(dist_cp)[[i]], " - ")[[1]][[2]]
  )

  data.frame(
    "twgss"    = plot_scree(dist_cp[[i]]$dist_mat, 1:10, "none"),
    "k"        = 1:10, # TODO Add to metadata
    "season"   = season,
    "years"    = years
  )
}))

twgss_vals |>
  arrange(factor(season, levels = seasons)) |>
  mutate(
    ind = ifelse(
      !is.na(years), paste0(season, " - ", years), season
    ),
    ind = factor(
      ind,
      levels = c(
        "Winter - 1960-1998", "Winter - 1999-2020",
        "Spring",
        "Summer - 1960-1990", "Summer - 1991-2020",
        "Autumn"
      )
    )
  ) |>
  ggplot(aes(x = k, y = twgss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = unique(twgss_vals$k)) +
  facet_wrap(~ind, scale = "free_y") +
  labs(y = "TWD") +
  cecl_theme()

# TODO Fit for k = 2 and k = 4 as well, even just so as to save
k <- 3 # Do for now, but might need to test for others later ...

# cluster
clust_cp <- lapply(dist_cp, \(x) {
  cecl_clust(
      x,
      k      = k,
      ncores = getOption("mc_cores", 2L)
    )
})

# function to ensure colours in map plots are always aligned
align_pam_labels <- \(reference_pam, target_pam) {
  reference <- reference_pam$clustering
  target <- target_pam$clustering

  if (is.null(names(reference)) || is.null(names(target))) {
    stop(
      "Both cluster-membership vectors must be named using station names."
    )
  }

  common_names <- intersect(names(reference), names(target))

  if (length(common_names) == 0L) {
    stop(
      "Reference and target clusterings have no station names in common."
    )
  }

  reference_levels <- sort(unique(reference))
  target_levels <- sort(unique(target))

  if (length(reference_levels) != length(target_levels)) {
    stop(
      "Reference and target solutions have different numbers of clusters."
    )
  }

  overlap <- table(
    reference = factor(
      reference[common_names],
      levels = reference_levels
    ),
    target = factor(
      target[common_names],
      levels = target_levels
    )
  )

  assignment <- clue::solve_LSAP(
    overlap,
    maximum = TRUE
  )

  # Names are target labels; values are corresponding reference labels
  label_map <- stats::setNames(
    reference_levels,
    target_levels[as.integer(assignment)]
  )

  aligned_labels <- label_map[as.character(target)]
  aligned_labels <- as.integer(aligned_labels)

  # Restore station names after integer coercion
  names(aligned_labels) <- names(target)

  aligned_pam <- target_pam
  aligned_pam$clustering <- aligned_labels

  attr(aligned_pam, "cluster_label_map") <- label_map
  attr(aligned_pam, "cluster_overlap") <- overlap

  aligned_pam
}

# align colours
clust_cp_aligned <- clust_cp
for (i in seq_along(clust_cp_aligned)) {
  if (i != reference_i) {
    clust_cp_aligned[[i]]$pam <- align_pam_labels(
      reference_pam = clust_cp[[reference_i]]$pam,
      target_pam = clust_cp[[i]]$pam
    )
  }
}

# assign distinct colours for each cluster
cluster_cols <- stats::setNames(
  ggsci::pal_nejm()(k),
  as.character(seq_len(k))
)

# map_plots_cp <- mclapply(seq_along(clust_cp_aligned), \(i) {
map_plots_cp <- lapply(seq_along(clust_cp_aligned), \(i) {
  season <- as.character(marg_laplace_cp[[i]]$season)

  years <- ifelse(
    season == names(dist_cp)[[i]],
    NA,
    str_split(names(dist_cp)[[i]], " - ")[[1]][[2]]
  )

  ind <- ifelse(
    !is.na(years),
    paste0(season, " - ", years),
    season
  )

  p <- plt_clust_map(
    pts = pts,
    areas = areas_proj,
    clust_obj = clust_cp_aligned[[i]]$pam,
    plot_medoids = FALSE,
    pt_size = 6,
    elev_df = elev_df,
    point_cols = cluster_cols,
    alpha = 1
  ) +
    ggtitle(ind)

  if (i %in% c(1, 2, 3)) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  if (i %in% c(2, 3, 5, 6)) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  p
})

names(map_plots_cp) <- names(clust_cp_aligned)

map_plots_join <- wrap_plots(map_plots_cp)

map_plots_join

# ggsave(plot = map_plots_join, "latex/plots/clust_map.png", width = 12, height = 8)
ggsave(plot = map_plots_join, paste0("latex/plots/clust_map_dqu_", dqu, ".png"), width = 12, height = 8)


#### Results for different dissimilarity matrices?? ####
#### DQU Sensitivity Plot (for 80 and 85%) ####

#### Supplementary Materials ####
#### TODO Parameter maps before/after changepoint ####
#### TODO Gaussian QQ plots ####



#### Old ####
#### Chi plot ####

# TODO Chi plot:
# TODO - Functionalise plotting code
# TODO - Only have chi plot (not chibar)
# TODO - Add ability to specify years we want to include
# TODO - Add elevation

chi_settings_df <- tidyr::crossing(
  "season" = c("Winter", "Spring", "Summer", "Autumn"),
  "combo"  = c("temp - drought_local"),
) |>
  separate_wider_delim(combo, delim = " - ", names = c("var1", "var2")) |>
  arrange(season, var1, var2)
chi_settings_df

# TODO Functionalise
plot_list <- mclapply(seq_len(nrow(chi_settings_df)), \(i) {
  row <- chi_settings_df[i, ]
  season <- row[["season"]]
  var1 <- row[["var1"]]
  var2 <- row[["var2"]]
  # decade <- row[["decade"]]

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

  # decade_lab <- ""
  # if (!is.na(decade)) {
  #   data_spec <- filter(data_spec, decade == !!decade)
  #   decade_lab <- paste0(", ", decade, "s")
  # }

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
    # mutate(var1 = !!var1, var2 = !!var2, season = !!season, decade = !!decade)
    mutate(var1 = !!var1, var2 = !!var2, season = !!season)

  p_spec <- join_chi_plots(chi_95_sf, areas_ccaa, diff_fill_fun, diff_fill_fun) +
    # plot_annotation(title = "Temperature vs Precipitation, Winter")
    plot_annotation(title = paste0(
      var1_lab, " vs ", var2_lab, season_lab # , decade_lab
    ))

  return(list("data" = chi_95_sf_spec, "plot" = p_spec))
})
names(plot_list) <- with(
  chi_settings_df,
  paste(
    var1,
    var2,
    ifelse(is.na(season), "Year", season),
    # ifelse(is.na(decade), "All", decade),
    sep = " - "
  )
)

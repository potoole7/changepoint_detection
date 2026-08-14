# TODO Make this into stand alone script! Missing code to load areas

# get raster of elevations from areas file
get_elev_rast <- \(areas, z = 9, bins = NULL, labels = NULL, reduce = TRUE) {
  # convert to projected CRS  # areas_proj <- st_transform(areas, IrishGrid = 29902)
  # for Spain this time, not Ireland!
  areas_proj <- st_transform(areas, crs = 25830)

  # pull elevation data from Open Elevation API (DEM = Digital Elevation Model)
  dem <- terra::rast(elevatr::get_elev_raster(
    locations = areas_proj, z = z, clip = "locations"
  ))

  # Mask DEM to the exact MULTIPOLYGON footprint
  dem_masked <- terra::mask(dem, terra::vect(areas_proj))

  # If specified, 8 × 8 cells become one cell: approximately 64× fewer rows
  dem_reduced <- dem_masked
  if (reduce) {
    dem_reduced <- terra::aggregate(
      dem_masked,
      fact = 8,
      fun = mean,
      na.rm = TRUE
    )
  }

  # Convert raster to a data.frame for ggplot
  df_dem <- as.data.frame(dem_reduced, xy = TRUE, na.rm = TRUE)
  names(df_dem) <- c("x", "y", "elevation")
  # filter out negative elevations
  df_dem <- filter(df_dem, elevation >= 0)

  if (!is.null(bins)) {
    df_dem <- df_dem %>%
      mutate(
        elev_bin = factor(cut(
          elevation,
          # breaks = c(seq(0, 500, by = 100), Inf),
          breaks = bins,
          # labels = c(
          #   "200–300 m", "300–400 m", "400–500 m", "> 500 m"
          # ),
          labels = labels,
          right  = FALSE
        ), levels = labels)
      )
  }

  return(df_dem)
}

elev_df <- get_elev_rast(areas_ccaa)

# ensure elev_df has the same bounding box as areas_ccaa
areas_proj <- st_transform(areas_ccaa, crs = 25830)
elev_df_crop <- elev_df %>%
  filter(
    x >= st_bbox(areas_proj)$xmin,
    x <= st_bbox(areas_proj)$xmax,
    y >= st_bbox(areas_proj)$ymin,
    y <= st_bbox(areas_proj)$ymax
  )

# areas_proj |>
#   ggplot() +
#   geom_sf(fill = NA, colour = "black")

p <- areas_proj |>
  ggplot() +
  geom_sf(fill = NA, colour = "black") +
  geom_tile(data = elev_df_crop, aes(x = x, y = y, fill = elevation)) +
  cecl_theme(nejm_pal = FALSE, legend.position = "right")

ggsave(plot = p, "spain_elevation.png", width = 10, height = 8)

summary(elev_df_crop$elevation)

elev_df_bin <- get_elev_rast(
  areas_ccaa,
  bins = c(0, 200, 600, 1000, 1500, 2000, Inf),
  labels = c(
    "<200",
    "200–600",
    "600–1000",
    "1000–1500",
    "1500–2000",
    ">2000"
  )
)

elev_df_bin_crop <- elev_df_bin %>%
  filter(
    x >= st_bbox(areas_proj)$xmin,
    x <= st_bbox(areas_proj)$xmax,
    y >= st_bbox(areas_proj)$ymin,
    y <= st_bbox(areas_proj)$ymax
  )

p1 <- areas_proj |>
  ggplot() +
  geom_sf(fill = NA, colour = "black") +
  geom_tile(
    data = elev_df_bin_crop,
    aes(x = x, y = y, fill = elev_bin),
    width = diff(range(elev_df_bin_crop$x)) / length(unique(elev_df_bin_crop$x)),
    height = diff(range(elev_df_bin_crop$y)) / length(unique(elev_df_bin_crop$y)),
  ) +
  cecl_theme(nejm_pal = FALSE, legend.position = "right") +
  labs(x = "", y = "") +
  # scico::scale_fill_scico_d(
  #   palette = "bilbao",
  #   name = "Elevation\n(m)",
  #   direction = -1
  # ) +
  scale_fill_manual(
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
  ) +
  NULL
# p1

# ggsave(plot = p1, "spain_elevation_binned.png", width = 10, height = 8)
ggsave(plot = p1, "spain_elevation_binned_reduced.png", width = 10, height = 8)

# save both objects
# readr::write_csv(elev_df_crop, "data/spain_elev.csv")
# readr::write_csv(elev_df_bin_crop, "data/spain_elev_binned.csv")
readr::write_csv(elev_df_crop, "data/spain_elev_reduced.csv")
readr::write_csv(elev_df_bin_crop, "data/spain_elev_binned_reduced.csv")

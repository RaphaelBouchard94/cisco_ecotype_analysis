library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
sf::sf_use_s2(FALSE)

# ---------------------------
# User inputs
# ---------------------------
bbox <- c(-79.5, 51, -78.35, 52)  # lon_min, lat_min, lon_max, lat_max
rup <- st_read("data/map/rupert/rupert.shp", quiet = TRUE) |> st_transform(4326)
not <- st_read("data/map/nottaway/nottaway.shp", quiet = TRUE) |> st_transform(4326)

was_pt <- st_as_sf(
  data.frame(name = "Waskaganish", lon = -78.76312, lat = 51.47979),
  coords = c("lon", "lat"),
  crs = 4326
)

smokey_pt <- st_as_sf(
  data.frame(name = "Waskaganish", lon = -78.513399, lat = 51.401414),
  coords = c("lon", "lat"),
  crs = 4326
)

kaa_pt <- st_as_sf(
  data.frame(name = "Waskaganish", lon = -78.808582, lat = 51.106724),
  coords = c("lon", "lat"),
  crs = 4326
)

# ---------------------------
# Natural Earth layers
# ---------------------------
bbox_sf <- st_as_sfc(st_bbox(c(xmin = bbox[1], ymin = bbox[2],
                               xmax = bbox[3], ymax = bbox[4]), crs = 4326))

# Land
land <- ne_countries(scale = "large", returnclass = "sf") |>
  st_crop(bbox_sf)

# Lakes (includes James Bay / Rupert Bay at large scale)
lakes <- ne_download(scale = "large", type = "lakes", category = "physical",
                     returnclass = "sf") |>
  st_crop(bbox_sf)

# Ocean / sea background
ocean <- ne_download(scale = "large", type = "ocean", category = "physical",
                     returnclass = "sf") |>
  st_crop(bbox_sf)

# Coastline (optional, for crisp edges)
coast <- ne_coastline(scale = "large", returnclass = "sf") |>
  st_crop(bbox_sf)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot() +
  # Ocean background
  geom_sf(data = ocean, fill = "white", color = NA) +
  # Land
  geom_sf(data = land, fill = "gray80", color = "gray80", linewidth = 0.3) +
  # Lakes (James Bay appears here)
  # geom_sf(data = lakes, fill = "#a8d8ea", linewidth = 0.3) +
  # Coastline
  # geom_sf(data = coast, color = "#5b9ab5", linewidth = 0.4) +
  # Rivers
  geom_sf(data = rup, color = "white", , fill = "white", linewidth = 0.8, lineend = "round") +
  geom_sf(data = not, color = "white", fill = "white",linewidth = 0.8, lineend = "round") +
  # # #Waskaganish
  # geom_sf(data = was_pt, shape = 21, fill = "black",
  #         size = 3.5, stroke = 0.8) +
  # geom_sf(data = smokey_pt, shape = 21, fill = "black",
  #         size = 3.5, stroke = 0.8) +
  # geom_sf(data = kaa_pt, shape = 21, fill = "black",
  #         size = 3.5, stroke = 0.8)+
  # Extent
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4]), expand = FALSE) +
  # Scale bar & north arrow
  annotation_scale(location = "bl", width_hint = 0.3,
                   style = "ticks", unit_category = "metric") +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering(),
                         height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
  labs(
    x = "Longitude", y = "Latitude"
  ) +
  scale_x_continuous(
    breaks = seq(-79.4, -78.4, by = 0.4)
  )+
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    plot.title = element_text(face = "bold"),
  )

print(p)

ggsave("figures/00_study_area_map.png", p, unit = "cm", width = 10, height = 15, dpi = 600)

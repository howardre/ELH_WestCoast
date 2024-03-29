# Libraries
library(ggOceanMaps)
library(usmap)
library(here)

# Create palettes and data frames
bathy_palette <- colorRampPalette(c("lightsteelblue1", "lightsteelblue4"))(9)
CA_bathy <- data.frame(lon = c(-126, -116),
                       lat = c(49, 32))
usa <- us_map(exclude = c("AK", "HI"))
usa_sf <- sf::st_as_sf(usa, 
                       coords = c("lon", "lat"),
                       crs = sf::st_crs(4326))

# Create labels for the states
text_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                          lat = c(47.5, 43.0, 37.0),
                          lon = c(-120.0, -121.0, -119.5))
text_sf <- sf::st_as_sf(text_labels, 
                        coords = c("lon", "lat"),
                        crs = sf::st_crs(4326))


# Make map
CCE_map <- basemap(data = CA_bathy,
                   bathymetry = TRUE,
                   rotate = TRUE,
                   legends = FALSE,
                   land.col = "wheat4",
                   grid.col = NA,
                   lon.interval = 4,
                   lat.interval = 3) +
  geom_polygon(data = transform_coord(CA_bathy),
               aes(x = lon, y = lat),
               fill = NA) +
  scale_fill_manual(values = bathy_palette) +
  labs(x = "Longitude",
       y = "Latitude")  +
  ggspatial::annotation_north_arrow(location = "tr",
                                    which_north = "true",
                                    style = ggspatial::north_arrow_nautical(text_family = "serif"),
                                    height = unit(3, "cm"),
                                    width = unit(3, "cm")) +
  # ggspatial::annotation_scale(location = "bl",
  #                             text_family = "serif",
  #                             style = "ticks",
  #                             height = unit(0.7, "cm"),
  #                             text_cex = 1.6,
  #                             line_width = 1.8) +
  geom_sf_text(data = text_sf,
               aes(label = name),
               size = 7,
               family = "serif") +
  theme(axis.text = element_text(family = "serif", size = 20),
        axis.title = element_text(family = "serif", size = 23),
        strip.text = element_text(family = "serif", size = 21))

CCE_map

# Save map
dev.copy(jpeg,
         here('manuscripts/Figures',
              'Figure1.jpg'),
         height = 18,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

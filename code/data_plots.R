### Title: Raw Data Plots
### Author: Rebecca Howard
### Date: 04/11/2024

library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB)
source(here('code/functions', 'read_data.R'))

# Load data
yoy_hake <- filter(read_data('yoy_hake.Rdata'), year > 2002) 
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42, year > 2013)
yoy_widow <- filter(read_data('yoy_widw.Rdata'), year > 2000)
yoy_shortbelly <- filter(read_data('yoy_sbly.Rdata'), year > 2000)
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), year > 2012) # doesn't work otherwise
yoy_squid <- read_data('yoy_squid.Rdata')

west_coast <- map_data("usa")

# Map function
make_maps <- function(data, size){
  ggplot() +  
    geom_polygon(aes(long, lat, group = group), data = west_coast,
                 fill = "lightyellow4", 
                 colour = "black") +
    geom_point(data = data, 
               aes(longitude, latitude,
                   size = size),
               color = "salmon") +
    coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 14),
          axis.title = element_text(family = "serif", size = 18),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          strip.text = element_text(family = "serif", size = 18),
          legend.title = element_text(family = "serif", size = 16),
          legend.text = element_text(family = "serif", size = 14)) +
    facet_wrap(~ year, nrow = 2) 
}

# Faceted maps
# Pacific Hake
make_maps(yoy_hake, yoy_hake$small)
make_maps(yoy_hake, yoy_hake$large)

# Maps for SVC terms
# Small sizes
# Low values of mean CU velocity
yoy_hake %>% filter(year == 2003 & 2004 & 2007 & 2010 & 2015 & 2016) %>% ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = west_coast,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(aes(longitude, latitude,
                 size = small),
             color = "salmon") +
  coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14))

# High values of mean CU velocity
yoy_hake %>% filter(year %in% 2006 & 2009 & 2011 & 2014 & 2018) %>% ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = west_coast,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(aes(longitude, latitude,
                 size = small),
             color = "salmon") +
  coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14))

# Large sizes
# Low values of 26 isopycnal depth
yoy_hake %>% filter(year == 2013 & 2007 & 2008 & 2009) %>% ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = west_coast,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(aes(longitude, latitude,
                 size = large),
             color = "salmon") +
  coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14))

# High values of depth of 26 isopycnal
yoy_hake %>% filter(year %in% 2010 & 2015 & 2016) %>% ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = west_coast,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(aes(longitude, latitude,
                 size = large),
             color = "salmon") +
  coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14))

# Shortbelly Rockfish
make_maps(yoy_shortbelly, yoy_shortbelly$small)
make_maps(yoy_shortbelly, yoy_shortbelly$large)

# Widow Rockfish
make_maps(yoy_widow, yoy_widow$small)
make_maps(yoy_widow, yoy_widow$large)

# Pacific Sanddab
make_maps(yoy_sdab, yoy_sdab$small)
make_maps(yoy_sdab, yoy_sdab$large)

# Market Squid
make_maps(yoy_squid, yoy_squid$small)
make_maps(yoy_squid, yoy_squid$large)

# Northern Anchovy
make_maps(yoy_anchovy, yoy_anchovy$small)
make_maps(yoy_anchovy, yoy_anchovy$large)
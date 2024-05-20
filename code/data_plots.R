### Title: Raw Data Plots
### Author: Rebecca Howard
### Date: 04/11/2024

library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB)
library(ggpubr)
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

size_maps <- function(data1, data2, length, plot_title, legend_title, max){
  ggplot(data = NULL) +
    geom_polygon(aes(long, lat, group = group), 
                 data = data1,
                 fill = "lightyellow4", 
                 colour = "black") +
    geom_point(data = data2,
               aes(longitude, 
                   latitude,
                   size = {{length}}),
               color = "salmon",
               alpha = 0.5) +
    coord_quickmap(xlim = c(-126, -117), ylim = c(32, 48)) +
    scale_size(limits = c(0, max)) +
    theme_classic() +
    labs(y = "Latitude \u00B0N",
         x = "Longitude \u00B0W",
         title = plot_title,
         size = legend_title) +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 14, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 14),
          axis.title = element_text(family = "serif", size = 14),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          legend.title = element_text(family = "serif", size = 13),
          legend.text = element_text(family = "serif", size = 11))
}

# Faceted maps
# Pacific Hake
make_maps(yoy_hake, yoy_hake$small)
make_maps(yoy_hake, yoy_hake$large)

# Maps for SVC terms
# Small sizes
yoy_hake_low <-  yoy_hake %>%  # Low values of mean CU velocity
  filter(year == 2003 & 2004 & 2007 & 2010 & 2015 & 2016)

yoy_hake_high <- yoy_hake %>%  # High values of mean CU velocity
  filter(year == 2006 & 2009 & 2011 & 2014 & 2018)

hake_low <- size_maps(west_coast, 
                      yoy_hake_low, 
                      small, 
                      "Low Mean California \nUndercurrent Velocity", 
                      "small size \nabundance",
                      max(yoy_hake_low$small))
  
hake_high <- size_maps(west_coast, 
                       yoy_hake_high, 
                       small, 
                       "High Mean California \nUndercurrent Velocity", 
                       "small size \nabundance",
                       max(yoy_hake_low$small)) 

ggarrange(hake_low, hake_high, nrow = 1, common.legend = TRUE, legend = "right")
  

# Large sizes
yoy_hake_low1 <-  yoy_hake %>%  # Low values of 26 isopycnal depth
  filter(year == 2013 & 2007 & 2008 & 2009)

yoy_hake_high1 <- yoy_hake %>%  # High values of 26 isopycnal depth
  filter(year == 2010 & 2015 & 2016)

hake_low1 <- size_maps(west_coast, 
                       yoy_hake_low1, 
                       large, 
                       "Shallower Depth of \nthe 26-isopycnal", 
                       "large size \nabundance",
                       max(yoy_hake_low1$large))

hake_high1 <- size_maps(west_coast, 
                        yoy_hake_high1, 
                        small, 
                        "Greater Depth of \nthe 26-isopycnal",
                        "large size \nabundance",
                        max(yoy_hake_low1$large)) 

ggarrange(hake_low1, hake_high1, nrow = 1, common.legend = TRUE, legend = "right")


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
### Title: sdmTMB Projecting
### Author: Rebecca Howard
### Date: 03/12/2024

# Load libraries ----
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB) 
library(visreg)
library(raster)
library(sf)
library(colorspace)
library(RColorBrewer)
library(mapdata)
library(fields)
library(ggpubr)
library(future)
library(scales)
library(tidync)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'distance_function.R'))

# Data
roms_ss <- readRDS(here('data', 'nep_ipsl.Rdata'))
roms_means <- readRDS(here('data', 'nep_ipsl_means.Rdata'))

yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42)
yoy_widow <- read_data('yoy_widw.Rdata') 
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36)
yoy_squid <- read_data('yoy_squid.Rdata')

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

extra_years <- c(2020:2100)

# Sandbox ----
data <- read_data('yoy_hake.Rdata') 

# Create spatiotemporal prediction grid
nlat = 40
nlon = 60
latd = seq(min(data$latitude), max(data$latitude), length.out = nlat)
lond = seq(min(data$longitude), max(data$longitude), length.out = nlon)

grid_extent <- expand.grid(lond, latd)
names(grid_extent) <- c('lon', 'lat')
grid_extent$dist <- NA # calculate distance from nearest station
for (k in 1:nrow(grid_extent)) {
  dist <-  distance_function(grid_extent$lat[k],
                             grid_extent$lon[k],
                             data$latitude,
                             data$longitude)
  grid_extent$dist[k] <- min(dist)
}

grid_extent$year <- 2025
grid_extent$week <- 20
grid_extent$year_week <- paste(grid_extent$year, grid_extent$week, sep = "-")
grid_extent$jday_scaled <- median(data$jday_scaled, na.rm = TRUE)
grid_extent <- add_utm_columns(grid_extent, c("lon", "lat"), 
                               utm_crs = 32610)
# Create fish and roms sf objects
grid_sf <- grid_extent %>%
  st_as_sf(coords = c("lon", "lat"), 
           remove = FALSE) %>%
  st_set_crs(32610)

ssroms_sf <- roms_ss %>%
  st_as_sf(coords = c("lon", "lat"),
           remove = FALSE) %>%
  st_set_crs(32610)

# Match SST and SSS roms to prediction grid
grid_combined <- do.call("rbind",
                        lapply(split(grid_sf, 1:nrow(grid_sf)), function(x) {
                          st_join(x, ssroms_sf[ssroms_sf$year_week == unique(x$year_week),],
                                  join = st_nearest_feature)
                        }))

grid_df <- as.data.frame(grid_combined)
grid_final <- grid_df %>%
  dplyr::select(-year_week.x, -year_week.y, -year.y, -lon.y, -lat.y, -geometry) %>%
  rename(year = year.x,
         longitude = lon.x,
         latitude = lat.x)

# Create latitude bins for different variables
roms_means$large_grp <- findInterval(roms_means$latitude,
                               c(32, 35, 38, 41, 44, 47))
roms_means$small_grp <- findInterval(roms_means$latitude,
                               c(32, 34, 36, 38, 40, 42, 44, 48))
grid_final$large_grp <- findInterval(grid_final$latitude,
                                  c(32, 35, 38, 41, 44, 47))
grid_final$small_grp <- findInterval(grid_final$latitude,
                                  c(32, 34, 36, 38, 40, 42, 44, 48))
nep_large <- roms_means %>% 
  group_by(years, large_grp) %>%
  summarize(across(c(vgeo, v_cu, vmax_cu), mean)) 
nep_small <- roms_means %>%
  group_by(years, small_grp) %>%
  summarize(across(c(u_vint_50m, u_vint_100m, depth_iso26, spice_iso26), mean))

large <- merge(grid_final,
               nep_large,
               by.x = c("year", "large_grp"),
               by.y = c("years", "large_grp"),
               all.x = TRUE)
small <- merge(large,
               nep_small,
               by.x = c("year", "small_grp"),
               by.y = c("years", "small_grp"),
               all.x = TRUE)

final <- dplyr::select(small, -small_grp, -large_grp)


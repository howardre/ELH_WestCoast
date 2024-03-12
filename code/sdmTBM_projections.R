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
grid_extent$jday_scaled <- median(data$jday_scaled, na.rm = TRUE)
grid_extent <- add_utm_columns(grid_extent, c("lon", "lat"), 
                               utm_crs = 32610)
grid_sf <- grid_extent %>%
  st_as_sf(coords = c("lon", "lat"), 
           remove = FALSE) %>%
  st_set_crs(32610)

ssroms_sf <- 

grid_extent$sst_scaled <- median(ss_roms$sst_scaled, na.rm = TRUE)
grid_extent$sss_scaled <- median(ss_roms$sss_scaled, na.rm = TRUE)
grid_extent$vgeo <- median(mean_roms$vgeo, na.rm = TRUE)
grid_extent$u_vint_50m <- median(mean_roms$u_vint_50m, na.rm = TRUE)
grid_extent$u_vint_100m <- median(mean_roms$u_vint_100m, na.rm = TRUE)
grid_extent$depth_iso26 <- median(mean_roms$depth_iso26, na.rm = TRUE)
grid_extent$spice_iso26 <- median(mean_roms$spice_iso26, na.rm = TRUE)
grid_extent$v_cu <- median(mean_roms$v_cu, na.rm = TRUE)
grid_extent$vmax_cu <- median(mean_roms$vmax_cu, na.rm = TRUE)



add_roms <- function(tows, roms){
  tows$week <- week(tows$date)
  tows$year_week <- paste(tows$year, tows$week, sep = "-")
  tow_sf <- tows %>%
    st_as_sf(coords = c("lon", "lat"), 
             remove = FALSE) %>%
    st_set_crs(32610)
  
  roms_sf <- roms %>%
    st_as_sf(coords = c("lon", "lat"),
             remove = FALSE) %>%
    st_set_crs(32610)
  tow_combined <- do.call("rbind",
                          lapply(split(tow_sf, 1:nrow(tow_sf)), function(x) {
                            st_join(x, roms_sf[roms_sf$year_week == unique(x$year_week),],
                                    join = st_nearest_feature)
                          }))
  
  tow_df <- as.data.frame(tow_combined)
  tow_final <- tow_df %>%
    dplyr::select(-year_week.x, -year_week.y, -year.y, -lon.y, -lat.y, -geometry) %>%
    rename(year = year.x,
           longitude = lon.x,
           latitude = lat.x)
  return(tow_final)
}

match_avgs <- function(species, roms){
  # Create latitude bins for different variables
  roms$large_grp <- findInterval(roms$latitude,
                                 c(32, 35, 38, 41, 44, 47))
  roms$small_grp <- findInterval(roms$latitude,
                                 c(32, 34, 36, 38, 40, 42, 44, 48))
  species$large_grp <- findInterval(species$latitude,
                                    c(32, 35, 38, 41, 44, 47))
  species$small_grp <- findInterval(species$latitude,
                                    c(32, 34, 36, 38, 40, 42, 44, 48))
  nep_large <- roms %>% 
    group_by(years, large_grp) %>%
    summarize(across(c(vgeo, v_cu, vmax_cu), mean)) 
  nep_small <- roms %>%
    group_by(years, small_grp) %>%
    summarize(across(c(u_vint_50m, u_vint_100m, depth_iso26, spice_iso26), mean))
  
  large <- merge(species,
                 nep_large,
                 by.x = c("year", "large_grp"),
                 by.y = c("years", "large_grp"),
                 all.x = TRUE)
  small <- merge(large,
                 nep_small,
                 by.x = c("year", "small_grp"),
                 by.y = c("years", "small_grp"),
                 all.x = TRUE)
  final <- dplyr::select(small, -small_grp, -large_grp, -haul_no, -cruise, -ctd_index,
                         -station_bottom_depth, -strata, -area,
                         -active, -numeric_date, -haul_date)
  return(final)
}
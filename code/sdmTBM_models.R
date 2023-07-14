### Title: sdmTMB Model Development
### Author: Rebecca Howard
### Date: 06/30/2023

# Load libraries ----
library(readxl)
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB)
library(visreg)
library(raster)
library(sf)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(roms_temperature, roms_salinity, roms_ssh, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year < 2020 &
             lat < 42) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10)
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  yoy_utm <- add_utm_columns(yoy, c("lon", "lat")) # add UTM coordinates
  return(yoy_utm)
}

# Data
yoy_hake <- read_data('yoy_hake.Rdata') 
# yoy_anchovy <- read_data('yoy_anch.Rdata') 
# yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
# yoy_widow <- read_data('yoy_widw.Rdata')
# yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
# yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
# yoy_sdab <- read_data('yoy_dab.Rdata') 

# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           cutoff = 10,
                           seed = 1993)
plot(yoy_hake_mesh)

# Fit model
# More info here: https://pbs-assess.github.io/sdmTMB/index.html
# Can add k-fold cross validation
hake_model <- sdmTMB(catch ~ 0 + as.factor(year) +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday) ,
                     # spatial_varying = ~ 0 + ssh_pos, # Not sure this is the right way to do this with SSH
                     data = yoy_hake,
                     mesh = yoy_hake_mesh,
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "rw")
sanity(hake_model)

tidy(hake_model)
tidy(hake_model,
     effect = "ran_pars", 
     conf.int = T)

# Plot covariates
visreg(hake_model, 
       xvar = "bottom_depth", 
       scale = "response")

# Predict and plot
# Need to make grid, may fix the varying coefficient issues below
# Prediction grid
nlat = 40
nlon = 60
latd = seq(min(yoy_hake$lat), max(yoy_hake$lat), length.out = nlat)
lond = seq(min(yoy_hake$lon), max(yoy_hake$lon), length.out = nlon)
spatial_grid <- expand.grid(lond, latd) # create grid
names(spatial_grid) <- c('lon', 'lat')
spatial_grid$dist <- NA # calculate distance from nearest station
for (k in 1:nrow(spatial_grid)) {
  dist <-  distance_function(spatial_grid$lat[k],
                             spatial_grid$lon[k],
                             yoy_hake$lat,
                             yoy_hake$lon)
  spatial_grid$dist[k] <- min(dist)
}
spatial_grid$year <- 2010
spatial_grid$bottom_depth <- median(yoy_hake$bottom_depth, na.rm = T)
spatial_grid$roms_temperature <- median(yoy_hake$roms_temperature, na.rm = T)
spatial_grid$roms_salinity <- median(yoy_hake$roms_salinity, na.rm = T)
spatial_grid$jday <- median(yoy_hake$jday, na.rm = T)
spatial_grid <- add_utm_columns(spatial_grid, c("lon", "lat"))

hake_pred <- predict(hake_model, newdata = spatial_grid, "response")
hake_pred$est[hake_pred$dist > 50000] <- NA # may want to find a way to mask with a polygon

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(hake_pred$est, na.rm = T),
                               max(hake_pred$est, na.rm = T))))
color_sequence = seq(max(hake_pred$est, na.rm = T), 
                     min(hake_pred$est, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(hake_pred$est, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

latd = seq(min(yoy_hake$Y), max(yoy_hake$Y), length.out = nlat)
lond = seq(min(yoy_hake$X), max(yoy_hake$X), length.out = nlon)

# Make map
windows(width = 7, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(yoy_hake$est,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(280, 1050),
      ylim = range(yoy_hake$Y, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(lond,
      latd,
      t(matrix(yoy_hake$est,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(280, 1050),
      ylim = range(yoy_hake$Y, na.rm = TRUE) + c(-.4, .5),
      cex.main = 2,
      cex.lab = 2,
      cex.axis = 1.8)
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T) # this does not work with UTM
# need to figure out how to get back to lat/lon
# Looks like there's something in the discussion on sdmTMB github 
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.18, .21, .17, .38),
           legend.cex = 1.3,
           axis.args = list(cex.axis = 1.6,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(hake_pred$est, na.rm = T), 
                    max(hake_pred$est, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.8,
                              line = 1.3,
                              family =  "serif"))

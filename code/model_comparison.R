### Title: ROMS vs. CTDs
### Author: Rebecca Howard
### Date: 12/12/2023

# Libraries ----
library(mgcv)
library(maps)
library(mapdata)
library(here)
library(purrr)
library(ggplot2)
library(Metrics)
library(dplyr)
library(fields)
library(colorspace)
library(ggthemes)

# Functions ----
source(here('code/functions', 'vis_gam_COLORS.R'))
read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(temperature, salinity, ssh_nc, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year > 1992 &
             year < 2015) %>%
    mutate(year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10)
  return(yoy)
}

plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam,
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = alpha('deepskyblue4', 0.4),
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       # seWithMean = T,
       scale = 0,
       cex.axis = 6.5,
       cex.lab = 6.5,
       family = "serif",
       lwd = 2.5)
}


# Load data ----
# See function for modifications made
# Currently only have hindcast up to 2018, so data is filtered to this
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata') 
yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
yoy_sdab <- read_data('yoy_dab.Rdata') 

# Hake ----
hake_small <- gam(small ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(temperature, k = 5) +
                    s(salinity, k = 5) + 
                    s(ssh_nc, k = 5),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake)
summary(hake_small)

hake_large <- gam(large ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(temperature, k = 5) +
                    s(salinity, k = 5) + 
                    s(ssh_nc, k = 5),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake)
summary(hake_large)

# Small model
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_small_insitu.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(hake_small,
              covariate = 2,
              bounds = c(-8, 3),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(hake_small,
              covariate = 4,
              bounds = c(-8, 3),
              "Temperature",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 5,
              bounds = c(-8, 3),
              "Salinity",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 6,
              bounds = c(-8, 3),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 3,
              bounds = c(-8, 3),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_large_insitu.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(hake_large,
              covariate = 2,
              bounds = c(-9, 8),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(hake_large,
              covariate = 4,
              bounds = c(-9, 8),
              "Temperature",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 5,
              bounds = c(-9, 8),
              "Salinity",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 6,
              bounds = c(-9, 8),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 3,
              bounds = c(-9, 8),
              "Day of Year",
              " ",
              "n")
dev.off()

# Widow ----
widow_small <- gam(small ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_widow)
summary(widow_small)

widow_large <- gam(large ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_widow)
summary(widow_large)

# Shortbelly ----
shortbelly_small <- gam(small ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_shortbelly)
summary(shortbelly_small)

shortbelly_large <- gam(large ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_shortbelly)
summary(shortbelly_large)

# Anchovy ----
anchovy_small <- gam(small ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_anchovy)
summary(anchovy_small)

anchovy_large <- gam(large ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(temperature, k = 5) +
                       s(salinity, k = 5) +
                       s(ssh_nc, k = 5),
                     family = tw(link = "log"),
                     method = "REML",
                     data = yoy_anchovy)
summary(anchovy_large)

# Sanddab ----
sdab_small <- gam(small ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(temperature, k = 5) +
                    s(salinity, k = 5) +
                    s(ssh_nc, k = 5),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_sdab)
summary(sdab_small)

sdab_large <- gam(large ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(temperature, k = 5) +
                    s(salinity, k = 5) +
                    s(ssh_nc, k = 5),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_sdab)
summary(sdab_large)

# Load libraries ----
library(maps)
library(mapdata)
library(ncdf4)
library(fields)
library(pracma)
library(viridis)
library(itsadug)
library(date)
library(here)
library(tidyverse)
library(lubridate)
library(magick)

# Load data and necessary functions ----
source(here('code/functions', 'distance_function.R'))

# NEP model output: avg bottom temp 2015-2019
nc <- nc_open(here('data','nep_srf_1995-2020.nc'))
lats <- ncvar_get(nc,"lat_rho")
lons <- ncvar_get(nc,"lon_rho")
lon1 <- ifelse(lons >= 180, lons -360, lons)
time <- ncvar_get(nc, varid = 'ocean_time')
time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
fillvalue_t <- ncatt_get(nc, 'temp', "_FillValue")
temp <- ncvar_get(nc, varid = 'temp')
nc_close(nc)

# Create figure of temperatures for 2015 - 2019 ----
# Get index for "month"
year <- unique(substr(time1, 1, 4))
month_year <- substr(time1, 1, 7)
month <- "02"

# Plot bottom temp in 'month'
windows(width = 9, height = 10)
par(mfrow = c(3, 2))
for (i in 1:length(year)) {
  starttime <- (1:length(month_year))[month_year == paste(year[i], month, sep = "-")][1]
  
  # Get temp
  tempplot <- temp[, , starttime]
  
  # Get date of image plot
  plotday1 <- as.character(time1[starttime])
  
  # show Temp
  image.plot(lon1,
             lats,
             tempplot,
             ylab = "Latitude (north)",
             col = viridis(100),
             xlab = "Longitude (east)",
             main = plotday1,
             xlim = c(-130, -116),
             ylim = c(30, 50))
  maps::map("world2",
            fill = T,
            col = "grey",
            add = T)
}

dev.copy(jpeg,
         'results/Surface_temp_2015_2019.jpg',
         height = 10,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Function for presentation figure
annual_temps <- function(number, month){
  tempplot <- temp[, , number]
  par(mfrow = c(1, 1),
      family = 'Lato',
      mar = c(6.4, 7.2, 4, 0.6) + 0.1,
      mgp = c(5, 2, 0))
  image.plot(lon1,
             lats,
             tempplot,
             ylab = "Latitude",
             col = viridis(100),
             xlab = "Longitude",
             main = "Monthly Temperature Change",
             legend.args = list("Temperature",
                                side = 4, 
                                cex = 7.6,
                                line = 6.2),
             xlim = c(-130, -116),
             ylim = c(30, 50),
             zlim =  c(-5, 23),
             cex.main = 7.7,
             cex.lab = 7.4,
             cex.axis = 7.4,
             legend.cex = 7.4,
             axis.args = list(cex.axis = 7.1))
  maps::map("world",
            fill = T,
            col = "wheat4",
            add = T)
  # image.plot(legend.only = T,
  #            col = viridis(100),
  #            legend.shrink = 0.2,
  #            smallplot = c(.76, .78, .64, .82),
  #            legend.cex = 0.8,
  #            axis.args = list(cex.axis = 0.8),
  #            legend.width = 0.5,
  #            legend.mar = 6,
  #            zlim = c(-5, 23),
  #            legend.args = list("Temperature",
  #                               side = 2, cex = 1))
  dev.copy(jpeg,
           paste('results/ROMS_hindcast/Surface_temp_', month, '_2015.jpg', sep = ''),
           height = 32,
           width = 20,
           res = 200,
           units = 'in',
           family = "Lato")
  dev.off()
}

# Presentation figure
# 1 = Jan, 5 = Feb, 9 = Mar, 14 = Apr, 18 = May, 23 = June
# 27 = July, 31 = Aug, 36 = Sept, 40 = Oct, 44 = Nov, 48 = Dec
windows(height = 28, width = 20, family = "serif")
annual_temps(1, '01')
windows(height = 28, width = 20, family = "serif")
annual_temps(5, '02')
windows(height = 28, width = 20, family = "serif")
annual_temps(9, '03')
windows(height = 28, width = 20, family = "serif")
annual_temps(14, '04')
windows(height = 28, width = 20, family = "serif")
annual_temps(18, '05')
windows(height = 28, width = 20, family = "serif")
annual_temps(23, '06')
windows(height = 28, width = 20, family = "serif")
annual_temps(27, '07')
windows(height = 28, width = 20, family = "serif")
annual_temps(31, '08')
windows(height = 28, width = 20, family = "serif")
annual_temps(36, '09')
windows(height = 28, width = 20, family = "serif")
annual_temps(40, '10')
windows(height = 28, width = 20, family = "serif")
annual_temps(44, '11')
windows(height = 28, width = 20, family = "serif")
annual_temps(48, '12')

# Gif
# Requires installation of Ghostscript and ImageMagick on machine
# Make sure Ghostscript is up to date
base_dir <- getwd()
temp_dir_out <- file.path(base_dir, 'results', 'ROMS', 'ROMS_hindcast')
temp_imgs <- list.files(temp_dir_out, full.names = TRUE)
temp_img_list <- lapply(temp_imgs, image_read)
temp_img_joined <- image_join(temp_img_list)
temp_img_animated <- image_animate(temp_img_joined, fps = 1)
image_write(image = temp_img_animated,
            path = here('results', 'ROMS', "temp_avgs.gif"))
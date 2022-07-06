# Libraries
library(marmap)

# Load fish data
yoy_hake_catch <- na.omit(readRDS(file = "../data/yoy_hake_catch.Rdata"))
yoy_hake_catch$presence <- 1 * (yoy_hake_catch$catch > 0)

# Load bathymetry
WC_bathy <- getNOAA.bathy(lon1 = -125.7, lon2 = -116.5, 
                          lat1 = 48.6, lat2 = 32.3, 
                          resolution = 1)
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))



windows(width = 20,
        height = 32,)
par(mfrow = c(1, 1),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
plot.bathy(
  WC_bathy,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  asp = NA,
  bpal = list(c(0,
                max(WC_bathy),
                greys),
              c(min(WC_bathy),
                0,
                blues)),
  xlim = c(-125.7,-116.5),
  ylim = range(yoy_hake_catch$lat, na.rm = TRUE) + c(-.4, .5),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  main = "",
  cex.lab = 1.2,
  cex.main = 1.7,
  cex.axis = 1)
points(yoy_hake_catch$lon[yoy_hake_catch$presence == 0 & yoy_hake_catch$year == 2015],
       yoy_hake_catch$lat[yoy_hake_catch$presence == 0 & yoy_hake_catch$year == 2015],
       pch = 4,
       col = 'red',
       cex = .9)
points(yoy_hake_catch$lon[yoy_hake_catch$presence == 1 & yoy_hake_catch$year == 2015],
       yoy_hake_catch$lat[yoy_hake_catch$presence == 1 & yoy_hake_catch$year == 2015],
       pch = 18,
       col = 'black',
       cex = .8)

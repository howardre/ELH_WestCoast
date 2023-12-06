sdmTMB_grid <- function(df, model){
  nlat = 40
  nlon = 60
  latd = seq(min(df$lat), max(df$lat), length.out = nlat)
  lond = seq(min(df$lon), max(df$lon), length.out = nlon)
  spatial_grid <- expand.grid(lond, latd) # create grid
  names(spatial_grid) <- c('lon', 'lat')
  spatial_grid$dist <- NA # calculate distance from nearest station
  for (k in 1:nrow(spatial_grid)) {
    dist <-  distance_function(spatial_grid$lat[k],
                               spatial_grid$lon[k],
                               df$lat,
                               df$lon)
    spatial_grid$dist[k] <- min(dist)
  }
  spatial_grid$year <- 2014
  spatial_grid$bottom_depth <- median(df$bottom_depth, na.rm = TRUE)
  spatial_grid$roms_temperature <- median(df$roms_temperature, na.rm = TRUE)
  spatial_grid$roms_salinity <- median(df$roms_salinity, na.rm = TRUE)
  spatial_grid$ssh_anom <- median(df$ssh_anom, na.rm = TRUE)
  spatial_grid$jday <- median(df$jday, na.rm = TRUE)
  spatial_grid$ssh_annual_scaled <- median(df$ssh_annual_scaled, na.rm = TRUE)
  spatial_grid <- add_utm_columns(spatial_grid, c("lon", "lat"), 
                                  utm_crs = 32610,)
  
  preds <- predict(model, 
                   newdata = spatial_grid, 
                   "link")
  preds$est[preds$dist > 50000] <- NA # may want to find a way to mask with a polygon
  return(preds)
}

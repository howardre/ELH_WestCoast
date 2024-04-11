sdmTMB_grid <- function(df, model, roms_large, roms_small, year){
  nlat = 40
  nlon = 60
  latd = seq(min(df$latitude), max(df$latitude), length.out = nlat)
  lond = seq(min(df$longitude), max(df$longitude), length.out = nlon)
  spatial_grid <- expand.grid(lond, latd) # create grid
  names(spatial_grid) <- c('longitude', 'latitude')
  spatial_grid$dist <- NA # calculate distance from nearest station
  for (k in 1:nrow(spatial_grid)) {
    dist <-  distance_function(spatial_grid$latitude[k],
                               spatial_grid$longitude[k],
                               df$latitude,
                               df$longitude)
    spatial_grid$dist[k] <- min(dist)
  }
  spatial_grid$year <- year
  spatial_grid$depth_scaled <- median(df$depth_scaled, na.rm = TRUE)
  spatial_grid$sst_scaled <- median(df$sst_scaled, na.rm = TRUE)
  spatial_grid$sss_scaled <- median(df$sss_scaled, na.rm = TRUE)
  spatial_grid$jday_scaled <- median(df$jday_scaled, na.rm = TRUE)
  spatial_grid$large_grp <- findInterval(spatial_grid$latitude,
                                         c(32, 36, 40, 44, 48))
  spatial_grid$small_grp <- findInterval(spatial_grid$latitude,
                                         c(32, 34, 36, 38, 40, 42, 44, 48))
  
  large <- merge(spatial_grid,
                 roms_large,
                 by.x = c("year", "large_grp"),
                 by.y = c("years", "large_grp"),
                 all.x = TRUE)
  small <- merge(large,
                 roms_small,
                 by.x = c("year", "small_grp"),
                 by.y = c("years", "small_grp"),
                 all.x = TRUE)
  final <- dplyr::select(small, -small_grp, -large_grp)
  
  final <- add_utm_columns(final, c("longitude", "latitude"), 
                                  utm_crs = 32610)
  
  preds <- predict(model, 
                   newdata = final, 
                   "link")
  preds$est[preds$dist > 50000] <- NA # may want to find a way to mask with a polygon
  preds$preds_scaled <- rescale(exp(preds$est))
  preds$preds_scaled_log <- rescale(preds$est)
  return(preds)
}

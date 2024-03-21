sdmTMB_grid <- function(df, model){
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
  spatial_grid$year <- 2020
  spatial_grid$depth_scaled <- median(df$depth_scaled, na.rm = TRUE)
  spatial_grid$sst_scaled <- median(df$sst_scaled, na.rm = TRUE)
  spatial_grid$sss_scaled <- median(df$sss_scaled, na.rm = TRUE)
  spatial_grid$jday_scaled <- median(df$jday_scaled, na.rm = TRUE)
  spatial_grid$vgeo <- median(df$vgeo, na.rm = TRUE)
  spatial_grid$u_vint_50m <- median(df$u_vint_50m, na.rm = TRUE)
  spatial_grid$u_vint_100m <- median(df$u_vint_100m, na.rm = TRUE)
  spatial_grid$depth_iso26 <- median(df$depth_iso26, na.rm = TRUE)
  spatial_grid$spice_iso26 <- median(df$spice_iso26, na.rm = TRUE)
  spatial_grid$v_cu <- median(df$v_cu, na.rm = TRUE)
  spatial_grid$vmax_cu <- median(df$vmax_cu, na.rm = TRUE)
  spatial_grid <- add_utm_columns(spatial_grid, c("longitude", "latitude"), 
                                  utm_crs = 32610)
  
  preds <- predict(model, 
                   newdata = spatial_grid, 
                   "link")
  preds$est[preds$dist > 50000] <- NA # may want to find a way to mask with a polygon
  preds$preds_scaled <- rescale(exp(preds$est))
  preds$preds_scaled_log <- rescale(preds$est)
  return(preds)
}

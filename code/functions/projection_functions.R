# Functions to forecast distributions

# Project for a specific year
sdm_project <- function(data, the_year, formula_small, formula_large,
                        roms_means, roms_ss){
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
  
  grid_extent$year <- the_year
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
  
  ssroms_sf <- filter(ssroms_sf, year == the_year)
  
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
  
  final <- dplyr::select(small, -small_grp, -large_grp) %>%
    mutate(sss_scaled = scale(sss)[, 1],
           sst_scaled = scale(sst)[, 1],
           vgeo = scale(vgeo)[, 1], 
           v_cu = scale(v_cu)[, 1], 
           vmax_cu = scale(vmax_cu)[, 1],
           u_vint_50m = scale(u_vint_50m)[, 1], 
           u_vint_100m = scale(u_vint_100m)[, 1], 
           depth_iso26 = scale(depth_iso26)[, 1], 
           spice_iso26 = scale(spice_iso26)[, 1])
  
  # Predict on forecasted output
  pred_small <- predict(formula_small,
                        newdata = final)
  pred_large <- predict(formula_large,
                        newdata = final)
  
  # Remove values outside of sample stations
  pred_small$small_pred <- rescale(exp(pred_small$est))
  pred_small$est[pred_small$dist > 60000] <- NA
  pred_small$pred <- rescale(pred_small$est)
  pred_small$zeta_s_vgeo[pred_small$dist > 60000] <- NA
  pred_small$omega_s[pred_small$dist > 60000] <- NA
  pred_small$est_non_rf[pred_small$dist > 60000] <- NA
  pred_small$est_rf[pred_small$dist > 60000] <- NA
  
  pred_large$est[pred_large$dist > 60000] <- NA
  pred_large$pred <- rescale(pred_large$est)
  pred_large$zeta_s_vgeo[pred_large$dist > 60000] <- NA
  pred_large$omega_s[pred_large$dist > 60000] <- NA
  pred_large$est_non_rf[pred_large$dist > 60000] <- NA
  pred_large$est_rf[pred_large$dist > 60000] <- NA
  
  # Calculate niche overlap
  pred_small$large_pred <- rescale(exp(pred_large$est))
  pred_small$p_small <- pred_small$pred / sum(pred_small$pred, na.rm = T)
  pred_small$p_large <- pred_small$large_pred / sum(pred_small$large_pred, na.rm = T)
  schoener <- 1 - 0.5 * sum(abs(pred_small$p_small - pred_small$p_large), na.rm = TRUE)
  
  final_list <- list(pred_small, pred_large, schoener)
  
  return(final_list)
}


## Loop through all the years
sdm_loop <- function(data, formula_small, 
                     formula_large, roms_means, 
                     roms_ss, range){
  grids <- list()
  for(j in range) {
    grid <- sdm_project(data, j, formula_small, 
                        formula_large, roms_means, roms_ss)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}


## Predict averages for each cell
sdm_cells <- function(data, formula_small, 
                      formula_large, roms_means, 
                      roms_ss, range){
  preds <- sdm_loop(data, formula_small, 
                    formula_large, roms_means, 
                    roms_ss, range)
  small_list <- lapply(preds, "[[", 1)
  large_list <- lapply(preds, "[[", 2)
  schoener <- lapply(preds, "[[", 3)
  small_df <- data.frame(lat = small_list[[1]]$latitude,
                         lon = small_list[[1]]$longitude,
                         avg_pred = rowMeans(do.call(cbind, lapply(small_list, "[", "pred"))))
  large_df <- data.frame(lat = large_list[[1]]$latitude,
                         lon = large_list[[1]]$longitude,
                         avg_pred = rowMeans(do.call(cbind, lapply(large_list, "[", "pred"))))
  small_df$pred_scaled <- rescale(small_df$avg_pred)
  large_df$pred_scaled <- rescale(large_df$avg_pred)
  final_list <- list(small_df, large_df, schoener)
  return(final_list)
}
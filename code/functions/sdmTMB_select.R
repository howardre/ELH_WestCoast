sdmTMB_select_small <- function(df, fish_mesh){
  sdm_full <- sdmTMB(small ~ 0 +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday, k = 15),
                     spatial_varying = ~ 0 + ssh_annual_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "ar1",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_sss <- sdmTMB(small ~ 0 +
                      s(bottom_depth, k = 5) +
                      s(roms_temperature, k = 5) +
                      s(jday, k = 15),
                    spatial_varying = ~ 0 + ssh_annual_scaled,
                    data = df,
                    mesh = fish_mesh,
                    spatial = "on",
                    time = "year",
                    family = tweedie(link = "log"),
                    spatiotemporal = "ar1",
                    control = sdmTMBcontrol(newton_loops = 1,
                                            nlminb_loops = 2))
  sdm_sst <- sdmTMB(small ~ 0 +
                      s(bottom_depth, k = 5) +
                      s(jday, k = 15),
                    spatial_varying = ~ 0 + ssh_annual_scaled,
                    data = df,
                    mesh = fish_mesh,
                    spatial = "on",
                    time = "year",
                    family = tweedie(link = "log"),
                    spatiotemporal = "ar1",
                    control = sdmTMBcontrol(newton_loops = 1,
                                            nlminb_loops = 2))
  sdm_jday <- sdmTMB(small ~ 0 +
                       s(jday, k = 15),
                     spatial_varying = ~ 0 + ssh_annual_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "ar1",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_list <- list(sdm_full, sdm_sss, sdm_sst, sdm_jday)
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list),
                                         function(x) min(sdm_list[[x]]$aic)))]]
  return_list <- list(sdm_list, best_sdm)
}

sdmTMB_select_large <- function(df, fish_mesh){
  sdm_full <- sdmTMB(large ~ 0 +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday, k = 15),
                     spatial_varying = ~ 0 + ssh_annual_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "ar1",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_sss <- sdmTMB(large ~ 0 +
                      s(bottom_depth, k = 5) +
                      s(roms_temperature, k = 5) +
                      s(jday, k = 15),
                    spatial_varying = ~ 0 + ssh_annual_scaled,
                    data = df,
                    mesh = fish_mesh,
                    spatial = "on",
                    time = "year",
                    family = tweedie(link = "log"),
                    spatiotemporal = "ar1",
                    control = sdmTMBcontrol(newton_loops = 1,
                                            nlminb_loops = 2))
  sdm_sst <- sdmTMB(large ~ 0 +
                      s(bottom_depth, k = 5) +
                      s(jday, k = 15),
                    spatial_varying = ~ 0 + ssh_annual_scaled,
                    data = df,
                    mesh = fish_mesh,
                    spatial = "on",
                    time = "year",
                    family = tweedie(link = "log"),
                    spatiotemporal = "ar1",
                    control = sdmTMBcontrol(newton_loops = 1,
                                            nlminb_loops = 2))
  sdm_jday <- sdmTMB(large ~ 0 +
                       s(jday, k = 15),
                     spatial_varying = ~ 0 + ssh_annual_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "ar1",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_list <- list(sdm_full, sdm_sss, sdm_sst, sdm_jday)
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list),
                                         function(x) min(sdm_list[[x]]$aic)))]]
  return_list <- list(sdm_list, best_sdm)
}

sdmTMB_select_small <- function(df, fish_mesh){
  try({
  sdm_v_cu <- sdmTMB(small ~ 0 + v_cu +
                       sst_scaled +
                       sss_scaled,
                     spatial_varying = ~ 0 + v_cu,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vgeo <- sdmTMB(small ~ 0 + vgeo +
                       sst_scaled +
                       sss_scaled,
                     spatial_varying = ~ 0 + vgeo,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vmax_cu <- sdmTMB(small ~ 0 + vmax_cu +
                          sst_scaled +
                          sss_scaled,
                        spatial_varying = ~ 0 + vmax_cu,
                        data = df,
                        mesh = fish_mesh,
                        spatial = "off",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "iid",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
  sdm_uvint50m <- sdmTMB(small ~ 0 + u_vint_50m +
                           sst_scaled +
                           sss_scaled,
                         spatial_varying = ~ 0 + u_vint_50m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "off",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_uvint100m <- sdmTMB(small ~ 0 + u_vint_100m +
                            sst_scaled +
                            sss_scaled,
                          spatial_varying = ~ 0 + u_vint_100m,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "off",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "iid",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2))
  sdm_iso26 <- sdmTMB(small ~ 0 + depth_iso26 +
                        sst_scaled +
                        sss_scaled,
                      spatial_varying = ~ 0 + depth_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "off",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_spice <- sdmTMB(small ~ 0 + spice_iso26 +
                        sst_scaled +
                        sss_scaled,
                      spatial_varying = ~ 0 + spice_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "off",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_base <- sdmTMB(small ~ 0 + 
                       sst_scaled +
                       sss_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  }) # allows function to move past non-working models
  sdm_list <- mget(ls(pattern = "sdm_*"))
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list), 
                                         function(x) AIC(sdm_list[[x]])))]]
  return_list <- list(sdm_list, best_sdm)
}

sdmTMB_select_large <- function(df, fish_mesh){
  try({
  sdm_v_cu <- sdmTMB(large ~ 0 + v_cu +
                       sst_scaled +
                       sss_scaled,
                     spatial_varying = ~ 0 + v_cu,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vgeo <- sdmTMB(large ~ 0 + vgeo +
                       sst_scaled +
                       sss_scaled,
                     spatial_varying = ~ 0 + vgeo,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vmax_cu <- sdmTMB(large ~ 0 + vmax_cu +
                          sst_scaled +
                          sss_scaled,
                        spatial_varying = ~ 0 + vmax_cu,
                        data = df,
                        mesh = fish_mesh,
                        spatial = "off",
                        time = "year",
                        family = tweedie(link = "log"),
                        priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 75, sigma_lt = 5)),
                        spatiotemporal = "iid",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
  sdm_uvint50m <- sdmTMB(large ~ 0 + u_vint_50m +
                           sst_scaled +
                           sss_scaled,
                         spatial_varying = ~ 0 + u_vint_50m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "off",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_uvint100m <- sdmTMB(large ~ 0 + u_vint_100m +
                            sst_scaled +
                            sss_scaled,
                         spatial_varying = ~ 0 + u_vint_100m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "off",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_iso26 <- sdmTMB(large ~ 0 + depth_iso26 +
                        sst_scaled +
                        sss_scaled,
                      spatial_varying = ~ 0 + depth_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "off",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_spice <- sdmTMB(large ~ 0 + spice_iso26 +
                        sst_scaled +
                        sss_scaled,
                      spatial_varying = ~ 0 + spice_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "off",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_base <- sdmTMB(large ~ 
                       sst_scaled +
                       sss_scaled,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "off",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  })
  sdm_list <- mget(ls(pattern = "sdm_*"))
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list), 
                                         function(x) AIC(sdm_list[[x]])))]]
  return_list <- list(sdm_list, best_sdm)
}

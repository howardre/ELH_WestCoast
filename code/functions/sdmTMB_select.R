sdmTMB_select_small <- function(df, fish_mesh){
  sdm_v_cu <- sdmTMB(small ~ v_cu +
                       s(sst_scaled, k = 5) +
                       s(sss_scaled, k = 5),
                     spatial_varying = ~ 0 + v_cu,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vgeo <- sdmTMB(small ~ vgeo +
                       s(sst_scaled, k = 5) +
                       s(sss_scaled, k = 5),
                     spatial_varying = ~ 0 + vgeo,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vmax_cu <- sdmTMB(small ~ vmax_cu +
                          s(sst_scaled, k = 5) +
                          s(sss_scaled, k = 5),
                        spatial_varying = ~ 0 + vmax_cu,
                        data = df,
                        mesh = fish_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "iid",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
  sdm_uvint50m <- sdmTMB(small ~ u_vint_50m +
                           s(sst_scaled, k = 5) +
                           s(sss_scaled, k = 5),
                         spatial_varying = ~ 0 + u_vint_50m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_uvint100m <- sdmTMB(small ~ u_vint_100m +
                            s(sst_scaled, k = 5) +
                            s(sss_scaled, k = 5),
                          spatial_varying = ~ 0 + u_vint_100m,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "iid",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2))
  sdm_iso26 <- sdmTMB(small ~ depth_iso26 +
                        s(sst_scaled, k = 5) +
                        s(sss_scaled, k = 5),
                      spatial_varying = ~ 0 + depth_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_spice <- sdmTMB(small ~ spice_iso26 +
                        s(sst_scaled, k = 5) +
                        s(sss_scaled, k = 5),
                      spatial_varying = ~ 0 + spice_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_list <- list(sdm_v_cu, sdm_vgeo, sdm_vmax_cu, sdm_uvint50m, sdm_uvint100m, sdm_iso26, sdm_spice)
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list), 
                                         function(x) AIC(sdm_list[[x]])))]]
  return_list <- list(sdm_list, best_sdm)
}

sdmTMB_select_large <- function(df, fish_mesh){
  sdm_v_cu <- sdmTMB(large ~ v_cu +
                       s(sst_scaled, k = 5) +
                       s(sss_scaled, k = 5),
                     spatial_varying = ~ 0 + v_cu,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vgeo <- sdmTMB(large ~ vgeo +
                       s(sst_scaled, k = 5) +
                       s(sss_scaled, k = 5),
                     spatial_varying = ~ 0 + vgeo,
                     data = df,
                     mesh = fish_mesh,
                     spatial = "on",
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "iid",
                     control = sdmTMBcontrol(newton_loops = 1,
                                             nlminb_loops = 2))
  sdm_vmax_cu <- sdmTMB(large ~ vmax_cu +
                          s(sst_scaled, k = 5) +
                          s(sss_scaled, k = 5),
                        spatial_varying = ~ 0 + vmax_cu,
                        data = df,
                        mesh = fish_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 75, sigma_lt = 5)),
                        spatiotemporal = "iid",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
  sdm_uvint50m <- sdmTMB(large ~ u_vint_50m +
                           s(sst_scaled, k = 5) +
                           s(sss_scaled, k = 5),
                         spatial_varying = ~ 0 + u_vint_50m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_uvint100m <- sdmTMB(large ~ u_vint_100m +
                           s(sst_scaled, k = 5) +
                           s(sss_scaled, k = 5),
                         spatial_varying = ~ 0 + u_vint_100m,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "iid",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2))
  sdm_iso26 <- sdmTMB(large ~ depth_iso26 +
                        s(sst_scaled, k = 5) +
                        s(sss_scaled, k = 5),
                      spatial_varying = ~ 0 + depth_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_spice <- sdmTMB(large ~ spice_iso26 +
                        s(sst_scaled, k = 5) +
                        s(sss_scaled, k = 5),
                      spatial_varying = ~ 0 + spice_iso26,
                      data = df,
                      mesh = fish_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "iid",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
  sdm_list <- list(sdm_v_cu, sdm_vgeo, sdm_vmax_cu, sdm_uvint50m, sdm_uvint100m, sdm_iso26, sdm_spice)
  best_sdm <- sdm_list[[which.min(sapply(1:length(sdm_list), 
                                         function(x) AIC(sdm_list[[x]])))]]
  return_list <- list(sdm_list, best_sdm)
}

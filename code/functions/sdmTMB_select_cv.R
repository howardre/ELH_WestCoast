sdmTMB_cv_small <- function(df, fish_mesh){
  set.seed(1993)
  sdm_v_cu <- try(sdmTMB_cv(small ~ 0 + v_cu +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + v_cu,
                            data = df,
                            mesh = fish_mesh,
                            time = "year",
                            spatial = "on",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  sdm_vgeo <- try(sdmTMB_cv(small ~ 0 + vgeo +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + vgeo,
                            data = df,
                            mesh = fish_mesh,
                            time = "year",
                            spatial = "on",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  sdm_vmax_cu <- try(sdmTMB_cv(small ~ 0 + vmax_cu +
                                 s(jday_scaled, k = 3) +
                                 s(sst_scaled, k = 3) +
                                 s(sss_scaled, k = 3),
                               spatial_varying = ~ 0 + vmax_cu,
                               data = df,
                               mesh = fish_mesh,
                               spatial = "on",
                               time = "year",
                               family = tweedie(link = "log"),
                               spatiotemporal = "off",
                               k_folds = 5,
                               parallel = TRUE),
                     silent = TRUE)
  sdm_uvint50m <- try(sdmTMB_cv(small ~ 0 + u_vint_50m +
                                  s(jday_scaled, k = 3) +
                                  s(sst_scaled, k = 3) +
                                  s(sss_scaled, k = 3),
                                spatial_varying = ~ 0 + u_vint_50m,
                                data = df,
                                mesh = fish_mesh,
                                spatial = "on",
                                time = "year",
                                family = tweedie(link = "log"),
                                spatiotemporal = "off",
                                k_folds = 5,
                                parallel = TRUE),
                      silent = TRUE)
  sdm_uvint100m <- try(sdmTMB_cv(small ~ 0 + u_vint_100m +
                                   s(jday_scaled, k = 3) +
                                   s(sst_scaled, k = 3) +
                                   s(sss_scaled, k = 3),
                                 spatial_varying = ~ 0 + u_vint_100m,
                                 data = df,
                                 mesh = fish_mesh,
                                 spatial = "on",
                                 time = "year",
                                 family = tweedie(link = "log"),
                                 spatiotemporal = "off",
                                 k_folds = 5,
                                 parallel = TRUE),
                       silent = TRUE)
  sdm_iso26 <- try(sdmTMB_cv(small ~ 0 + depth_iso26 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + depth_iso26,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             k_folds = 5,
                             parallel = TRUE),
                   silent = TRUE)
  sdm_spice <- try(sdmTMB_cv(small ~ 0 + spice_iso26 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + spice_iso26,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             k_folds = 5,
                             parallel = TRUE),
                   silent = TRUE)
  sdm_base <- try(sdmTMB_cv(small ~ 0 +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            data = df,
                            mesh = fish_mesh,
                            spatial = "on",
                            time = "year",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  cv_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  best_model <- cv_list[[which.min(sapply(1:length(cv_list),
                                      function(x) cv_list[[x]]$sum_loglik))]]
  the_list <- list(cv_list, best_model)
  return(the_list)
  }

sdmTMB_cv_large <- function(df, fish_mesh){
  set.seed(1993)
  sdm_v_cu <- try(sdmTMB_cv(large ~ 0 + v_cu +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + v_cu,
                            data = df,
                            mesh = fish_mesh,
                            time = "year",
                            spatial = "on",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  sdm_vgeo <- try(sdmTMB_cv(large ~ 0 + vgeo +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + vgeo,
                            data = df,
                            mesh = fish_mesh,
                            time = "year",
                            spatial = "on",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  sdm_vmax_cu <- try(sdmTMB_cv(large ~ 0 + vmax_cu +
                                 s(jday_scaled, k = 3) +
                                 s(sst_scaled, k = 3) +
                                 s(sss_scaled, k = 3),
                               spatial_varying = ~ 0 + vmax_cu,
                               data = df,
                               mesh = fish_mesh,
                               spatial = "on",
                               time = "year",
                               family = tweedie(link = "log"),
                               spatiotemporal = "off",
                               k_folds = 5,
                               parallel = TRUE),
                     silent = TRUE)
  sdm_uvint50m <- try(sdmTMB_cv(large ~ 0 + u_vint_50m +
                                  s(jday_scaled, k = 3) +
                                  s(sst_scaled, k = 3) +
                                  s(sss_scaled, k = 3),
                                spatial_varying = ~ 0 + u_vint_50m,
                                data = df,
                                mesh = fish_mesh,
                                spatial = "on",
                                time = "year",
                                family = tweedie(link = "log"),
                                spatiotemporal = "off",
                                k_folds = 5,
                                parallel = TRUE),
                      silent = TRUE)
  sdm_uvint100m <- try(sdmTMB_cv(large ~ 0 + u_vint_100m +
                                   s(jday_scaled, k = 3) +
                                   s(sst_scaled, k = 3) +
                                   s(sss_scaled, k = 3),
                                 spatial_varying = ~ 0 + u_vint_100m,
                                 data = df,
                                 mesh = fish_mesh,
                                 spatial = "on",
                                 time = "year",
                                 family = tweedie(link = "log"),
                                 spatiotemporal = "off",
                                 k_folds = 5,
                                 parallel = TRUE),
                       silent = TRUE)
  sdm_iso26 <- try(sdmTMB_cv(large ~ 0 + depth_iso26 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + depth_iso26,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             k_folds = 5,
                             parallel = TRUE),
                   silent = TRUE)
  sdm_spice <- try(sdmTMB_cv(large ~ 0 + spice_iso26 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + spice_iso26,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             k_folds = 5,
                             parallel = TRUE),
                   silent = TRUE)
  sdm_base <- try(sdmTMB_cv(large ~ 0 +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            data = df,
                            mesh = fish_mesh,
                            spatial = "on",
                            time = "year",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            k_folds = 5,
                            parallel = TRUE),
                  silent = TRUE)
  cv_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  best_model <- cv_list[[which.min(sapply(1:length(cv_list),
                                          function(x) cv_list[[x]]$sum_loglik))]]
  the_list <- list(cv_list, best_model)
  return(the_list)
}

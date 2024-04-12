sdmTMB_select_small <- function(df, fish_mesh){
  set.seed(1993)
  sdm_v_cu <- try(sdmTMB(small ~ 0 + v_cu +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + v_cu,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vgeo <- try(sdmTMB(small ~ 0 + vgeo +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + vgeo,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vmax_cu <- try(sdmTMB(small ~ 0 + vmax_cu +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + vmax_cu,
                            extra_time = extra_years,
                            data = df,
                            mesh = fish_mesh,
                            spatial = "on",
                            time = "year",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            control = sdmTMBcontrol(newton_loops = 1,
                                                    nlminb_loops = 2)),
                     silent = TRUE)
  sdm_uvint50m <- try(sdmTMB(small ~ 0 + u_vint_50m +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + u_vint_50m,
                             extra_time = extra_years,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             control = sdmTMBcontrol(newton_loops = 1,
                                                     nlminb_loops = 2)),
                      silent = TRUE)
  sdm_uvint100m <- try(sdmTMB(small ~ 0 + u_vint_100m +
                                s(jday_scaled, k = 3) +
                                s(sst_scaled, k = 3) +
                                s(sss_scaled, k = 3),
                              spatial_varying = ~ 0 + u_vint_100m,
                              extra_time = extra_years,
                              data = df,
                              mesh = fish_mesh,
                              spatial = "on",
                              time = "year",
                              family = tweedie(link = "log"),
                              spatiotemporal = "off",
                              control = sdmTMBcontrol(newton_loops = 1,
                                                      nlminb_loops = 2)),
                       silent = TRUE)
  sdm_iso26 <- try(sdmTMB(small ~ 0 + depth_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + depth_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_spice <- try(sdmTMB(small ~ 0 + spice_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + spice_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_base <- try(sdmTMB(small ~ 0 +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  return(sdm_list)
}

sdmTMB_select_large <- function(df, fish_mesh){
  set.seed(1993)
  sdm_v_cu <- try(sdmTMB(large ~ 0 + v_cu +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + v_cu,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vgeo <- try(sdmTMB(large ~ 0 + vgeo +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + vgeo,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vmax_cu <- try(sdmTMB(large ~ 0 + vmax_cu +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + vmax_cu,
                            extra_time = extra_years,
                            data = df,
                            mesh = fish_mesh,
                            spatial = "on",
                            time = "year",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            control = sdmTMBcontrol(newton_loops = 1,
                                                    nlminb_loops = 2)),
                     silent = TRUE)
  sdm_uvint50m <- try(sdmTMB(large ~ 0 + u_vint_50m +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + u_vint_50m,
                             extra_time = extra_years,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             control = sdmTMBcontrol(newton_loops = 1,
                                                     nlminb_loops = 2)),
                      silent = TRUE)
  sdm_uvint100m <- try(sdmTMB(large ~ 0 + u_vint_100m +
                                s(jday_scaled, k = 3) +
                                s(sst_scaled, k = 3) +
                                s(sss_scaled, k = 3),
                              spatial_varying = ~ 0 + u_vint_100m,
                              extra_time = extra_years,
                              data = df,
                              mesh = fish_mesh,
                              spatial = "on",
                              time = "year",
                              family = tweedie(link = "log"),
                              spatiotemporal = "off",
                              control = sdmTMBcontrol(newton_loops = 1,
                                                      nlminb_loops = 2)),
                       silent = TRUE)
  sdm_iso26 <- try(sdmTMB(large ~ 0 + depth_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + depth_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_spice <- try(sdmTMB(large ~ 0 + spice_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + spice_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_base <- try(sdmTMB(large ~ 0 +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  return(sdm_list)
}

sdmTMB_select <- function(df, fish_mesh){
  set.seed(1993)
  sdm_v_cu <- try(sdmTMB(catch ~ 0 + v_cu +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + v_cu,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vgeo <- try(sdmTMB(catch ~ 0 + vgeo +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         spatial_varying = ~ 0 + vgeo,
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_vmax_cu <- try(sdmTMB(catch ~ 0 + vmax_cu +
                              s(jday_scaled, k = 3) +
                              s(sst_scaled, k = 3) +
                              s(sss_scaled, k = 3),
                            spatial_varying = ~ 0 + vmax_cu,
                            extra_time = extra_years,
                            data = df,
                            mesh = fish_mesh,
                            spatial = "on",
                            time = "year",
                            family = tweedie(link = "log"),
                            spatiotemporal = "off",
                            control = sdmTMBcontrol(newton_loops = 1,
                                                    nlminb_loops = 2)),
                     silent = TRUE)
  sdm_uvint50m <- try(sdmTMB(catch ~ 0 + u_vint_50m +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + u_vint_50m,
                             extra_time = extra_years,
                             data = df,
                             mesh = fish_mesh,
                             spatial = "on",
                             time = "year",
                             family = tweedie(link = "log"),
                             spatiotemporal = "off",
                             control = sdmTMBcontrol(newton_loops = 1,
                                                     nlminb_loops = 2)),
                      silent = TRUE)
  sdm_uvint100m <- try(sdmTMB(catch ~ 0 + u_vint_100m +
                                s(jday_scaled, k = 3) +
                                s(sst_scaled, k = 3) +
                                s(sss_scaled, k = 3),
                              spatial_varying = ~ 0 + u_vint_100m,
                              extra_time = extra_years,
                              data = df,
                              mesh = fish_mesh,
                              spatial = "on",
                              time = "year",
                              family = tweedie(link = "log"),
                              spatiotemporal = "off",
                              control = sdmTMBcontrol(newton_loops = 1,
                                                      nlminb_loops = 2)),
                       silent = TRUE)
  sdm_iso26 <- try(sdmTMB(catch ~ 0 + depth_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + depth_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_spice <- try(sdmTMB(catch ~ 0 + spice_iso26 +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + spice_iso26,
                          extra_time = extra_years,
                          data = df,
                          mesh = fish_mesh,
                          spatial = "on",
                          time = "year",
                          family = tweedie(link = "log"),
                          spatiotemporal = "off",
                          control = sdmTMBcontrol(newton_loops = 1,
                                                  nlminb_loops = 2)),
                   silent = TRUE)
  sdm_base <- try(sdmTMB(catch ~ 0 +
                           s(jday_scaled, k = 3) +
                           s(sst_scaled, k = 3) +
                           s(sss_scaled, k = 3),
                         extra_time = extra_years,
                         data = df,
                         mesh = fish_mesh,
                         spatial = "on",
                         time = "year",
                         family = tweedie(link = "log"),
                         spatiotemporal = "off",
                         control = sdmTMBcontrol(newton_loops = 1,
                                                 nlminb_loops = 2)),
                  silent = TRUE)
  sdm_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  return(sdm_list)
}

calc_stat_large <- function(model_list, fish_mesh, df){
  null_mod <- sdmTMB::sdmTMB(large ~ 1,
                         spatial = "off",
                         mesh = fish_mesh,
                         data = df)
  null_dev <- -2 * as.numeric(logLik(null_mod))
  
  log_lik <- data.frame(NA)
  for(i in 1:length(model_list)) {
    log_lik[i, ] <- as.numeric(logLik(model_list[[i]]))
  }
  
  rownames(log_lik) <- names(model_list)
  log_lik[, 2] <- -2 * log_lik
  colnames(log_lik) <- c("log_likelihood", "resid_dev")
  log_lik$deviance <- 100 * (null_dev - log_lik$resid_dev) / null_dev

  return(log_lik)
  }

calc_stat_small <- function(model_list, fish_mesh, df){
  null_mod <- sdmTMB::sdmTMB(small ~ 1,
                             spatial = "off",
                             mesh = fish_mesh,
                             data = df)
  null_dev <- -2 * as.numeric(logLik(null_mod))
  
  log_lik <- data.frame(NA)
  for(i in 1:length(model_list)) {
    log_lik[i, ] <- as.numeric(logLik(model_list[[i]]))
  }
  
  rownames(log_lik) <- names(model_list)
  log_lik[, 2] <- -2 * log_lik
  colnames(log_lik) <- c("log_likelihood", "resid_dev")
  log_lik$deviance <- 100 * (null_dev - log_lik$resid_dev) / null_dev
  
  return(log_lik)
}

calc_stat <- function(model_list, fish_mesh, df){
  null_mod <- sdmTMB::sdmTMB(catch ~ 1,
                             spatial = "off",
                             mesh = fish_mesh,
                             data = df)
  null_dev <- -2 * as.numeric(logLik(null_mod))
  
  log_lik <- data.frame(NA)
  for(i in 1:length(model_list)) {
    log_lik[i, ] <- as.numeric(logLik(model_list[[i]]))
  }
  
  rownames(log_lik) <- names(model_list)
  log_lik[, 2] <- -2 * log_lik
  colnames(log_lik) <- c("log_likelihood", "resid_dev")
  log_lik$deviance <- 100 * (null_dev - log_lik$resid_dev) / null_dev
  
  return(log_lik)
}

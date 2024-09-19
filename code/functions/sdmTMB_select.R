sdmTMB_select_small <- function(df, fish_mesh, size){
  set.seed(1993)
  if(size == "small"){
    sdm_formula = as.formula(small ~ 0 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3))}
  if(size == "large"){
    sdm_formula = as.formula(large ~ 0 +
                               s(jday_scaled, k = 3) +
                               s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3))}
  sdm_base <- try(sdmTMB(formula = sdm_formula,
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
  v_cu_formula <- update(sdm_formula, . ~ . + v_cu)
  sdm_v_cu <- try(sdmTMB(v_cu_formula,
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
  vgeo_formula <- update(sdm_formula, . ~ . + vgeo)
  sdm_vgeo <- try(sdmTMB(vgeo_formula,
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
  vmax_cu_formula <- update(sdm_formula, . ~ . + vmax_cu)
  sdm_vmax_cu <- try(sdmTMB(vmax_cu_formula,
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
  uvint50m_formula <- update(sdm_formula, . ~ . + u_vint_50m)
  sdm_uvint50m <- try(sdmTMB(uvint50m_formula,
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
  uvint100m_formula <- update(sdm_formula, . ~ . + u_vint_100m)
  sdm_uvint100m <- try(sdmTMB(uvint100m_formula,
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
  iso26_formula <- update(sdm_formula, . ~ . + depth_iso26)
  sdm_iso26 <- try(sdmTMB(iso26_formula,
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
  spice_formula <- update(sdm_formula, . ~ . + spice_iso26)
  sdm_spice <- try(sdmTMB(spice_formula,
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

  sdm_list <- Filter(is.list, mget(ls(pattern = "sdm_*")))
  return(sdm_list)
}



calc_stat <- function(model_list, fish_mesh, df, size){
  if(size == "small"){
    sdm_formula = small ~ 1
  }
  
  if(size == "large"){
    sdm_formula = large ~ 1
  }
  null_mod <- sdmTMB::sdmTMB(sdm_formula,
                         spatial = "off",
                         mesh = fish_mesh,
                         data = df)
  null_dev <- -2 * as.numeric(logLik(null_mod))
  
  log_lik <- data.frame(NA)
  for(i in 1:length(model_list)) {
    log_lik[i, ] <- as.numeric(logLik(model_list[[i]]))
  }
  
  aic <- data.frame(NA)
  for(i in 1:length(model_list)) {
    aic[i, ] <- as.numeric(AIC(model_list[[i]]))
  }
  
  colnames(aic) <- "AIC"
  
  rownames(log_lik) <- names(model_list)
  log_lik[, 2] <- -2 * log_lik
  colnames(log_lik) <- c("log_likelihood", "resid_dev")
  log_lik$deviance <- 100 * (null_dev - log_lik$resid_dev) / null_dev

  final_stat <- cbind(log_lik, aic)
  
  return(final_stat)
}

gam_select_small <- function(df){
  gam_v_cu <- gam(small ~ year_f +
                      s(longitude, latitude) +
                      s(jday) +
                      s(sst, k = 5) +
                      s(sss, k = 5) +
                      s(longitude, latitude, by = v_cu),
                    family = tw(link = "log"),
                    method = "REML",
                    data = df)
  gam_vgeo <- gam(small ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = vgeo),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_vmax_cu <- gam(small ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = vmax_cu),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_uvint50m <- gam(small ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = u_vint_50m),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_uvint100m <- gam(small ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = u_vint_100m),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_iso26 <- gam(small ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = depth_iso26),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_spice <- gam(small ~ year_f +
                     s(longitude, latitude) +
                     s(jday) +
                     s(sst, k = 5) +
                     s(sss, k = 5) +
                     s(longitude, latitude, by = spice_iso26),
                   family = tw(link = "log"),
                   method = "REML",
                   data = df)
  gam_list <- list(gam_v_cu, gam_vgeo, gam_vmax_cu, gam_uvint50m, gam_uvint100m, gam_iso26, gam_spice)
  best_sdm <- gam_list[[which.min(sapply(1:length(gam_list), 
                                         function(x) AIC(gam_list[[x]])))]]
  return_list <- list(gam_list, best_sdm)
}

gam_select_large <- function(df){
  gam_v_cu <- gam(large ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = v_cu),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_vgeo <- gam(large ~ year_f +
                    s(longitude, latitude) +
                    s(jday) +
                    s(sst, k = 5) +
                    s(sss, k = 5) +
                    s(longitude, latitude, by = vgeo),
                  family = tw(link = "log"),
                  method = "REML",
                  data = df)
  gam_vmax_cu <- gam(large ~ year_f +
                       s(longitude, latitude) +
                       s(jday) +
                       s(sst, k = 5) +
                       s(sss, k = 5) +
                       s(longitude, latitude, by = vmax_cu),
                     family = tw(link = "log"),
                     method = "REML",
                     data = df)
  gam_uvint50m <- gam(large ~ year_f +
                        s(longitude, latitude) +
                        s(jday) +
                        s(sst, k = 5) +
                        s(sss, k = 5) +
                        s(longitude, latitude, by = u_vint_50m),
                      family = tw(link = "log"),
                      method = "REML",
                      data = df)
  gam_uvint100m <- gam(large ~ year_f +
                         s(longitude, latitude) +
                         s(jday) +
                         s(sst, k = 5) +
                         s(sss, k = 5) +
                         s(longitude, latitude, by = u_vint_100m),
                       family = tw(link = "log"),
                       method = "REML",
                       data = df)
  gam_iso26 <- gam(large ~ year_f +
                     s(longitude, latitude) +
                     s(jday) +
                     s(sst, k = 5) +
                     s(sss, k = 5) +
                     s(longitude, latitude, by = depth_iso26),
                   family = tw(link = "log"),
                   method = "REML",
                   data = df)
  gam_spice <- gam(large ~ year_f +
                       s(longitude, latitude) +
                       s(jday) +
                       s(sst, k = 5) +
                       s(sss, k = 5) +
                       s(longitude, latitude, by = spice_iso26),
                     family = tw(link = "log"),
                     method = "REML",
                     data = df)
  gam_list <- list(gam_v_cu, gam_vgeo, gam_vmax_cu, gam_uvint50m, gam_uvint100m, gam_iso26, gam_spice)
  best_sdm <- gam_list[[which.min(sapply(1:length(gam_list), 
                                         function(x) AIC(gam_list[[x]])))]]
  return_list <- list(gam_list, best_sdm)
}

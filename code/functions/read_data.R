read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(sst, sss, bottom_depth, year, jday, latitude, longitude) %>%
    filter(catch < 2500 &
             year < 2020 & year > 2010) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_scaled = scale(year)[, 1],
           year_f = as.factor(year),
           depth_scaled = scale(bottom_depth)[, 1],
           sss_scaled = scale(sss)[, 1],
           sst_scaled = scale(sst)[, 1],
           jday_scaled = scale(jday)[, 1],
           vgeo = scale(vgeo)[, 1], 
           v_cu = scale(v_cu)[, 1], 
           vmax_cu = scale(vmax_cu)[, 1],
           u_vint_50m = scale(u_vint_50m)[, 1], 
           u_vint_100m = scale(u_vint_100m)[, 1], 
           depth_iso26 = scale(depth_iso26)[, 1], 
           spice_iso26 = scale(spice_iso26)[, 1])
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  yoy_utm <- add_utm_columns(yoy, 
                             utm_crs = 32610, # UTM 10
                             ll_names = c("longitude", "latitude")) # add UTM coordinates
  return(yoy_utm)
}

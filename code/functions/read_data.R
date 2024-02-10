read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(roms_temperature, roms_salinity, roms_ssh, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year < 2020) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10,
           ssh_annual_scaled = scale((year_ssh - mean(year_ssh)) / sd(year_ssh))[, 1],
           depth_scaled = scale(bottom_depth)[, 1],
           sss_scaled = scale(roms_salinity)[, 1],
           sst_scaled = scale(roms_temperature)[, 1],
           ssh_scaled = scale(ssh_anom)[, 1])
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  yoy_utm <- add_utm_columns(yoy, 
                             utm_crs = 32610, # UTM 10
                             ll_names = c("lon", "lat")) # add UTM coordinates
  return(yoy_utm)
}

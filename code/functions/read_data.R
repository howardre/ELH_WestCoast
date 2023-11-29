read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(roms_temperature, roms_salinity, roms_ssh, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year < 2020 &
             roms_salinity > 33 &
             bottom_depth < 2000) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10,
           ssh_annual_scaled = (year_ssh - mean(year_ssh)) / sd(year_ssh),
           depth_scaled = (bottom_depth - mean(bottom_depth)) / sd(bottom_depth),
           sss_scaled = (roms_salinity - mean(roms_salinity)) / sd(roms_salinity),
           sst_scaled = (roms_temperature - mean(roms_temperature)) / sd(roms_temperature),
           ssh_scaled = (ssh_anom - mean(ssh_anom)) / sd(ssh_anom))
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  yoy_utm <- add_utm_columns(yoy, 
                             utm_crs = 32610, # UTM 10
                             ll_names = c("lon", "lat")) # add UTM coordinates
  return(yoy_utm)
}

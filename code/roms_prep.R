### Title: Projection Prep
### Author: Rebecca Howard
### Date: 03/12/2024

# Load libraries ----
library(tibble)
library(here)
library(dplyr)
library(tidync)
library(ncdf4)
library(lubridate)

# Prepare ROMS
nep_ipsl <- tidync('F:/NEP_ROMS_ipsl/NEP_ipsl_sstsss_WeeklyMean_April2July_1995-2100.nc')

nc <- nc_open('F:/NEP_ROMS_ipsl/NEP_ipsl_sstsss_WeeklyMean_April2July_1995-2100.nc')
lats <- ncvar_get(nc, "latitude")
lons <- ncvar_get(nc, "longitude")
nc_close(nc)

nep_col <- tidync('F:/NEP_ROMS_ipsl/NEP_ipsl_sstsss_WeeklyMean_April2July_1995-2100.nc')
nep_ss <- nep_col %>%
  hyper_tibble()

nep_ss$longitude <- lons[cbind(nep_ss$lon, nep_ss$lat)]
nep_ss$latitude <- lats[cbind(nep_ss$lon, nep_ss$lat)]

nep_time <- nep_col %>%
  activate("D2") %>%
  hyper_tibble

nep_time$date <- ymd(paste(nep_time$year, nep_time$month, nep_time$day, sep = "-"))
nep_time$week <- week(nep_time$date)
nep_time$year_week <- paste(nep_time$year, nep_time$week, sep = "-")

nep_merge <- merge(nep_ss, nep_time, by = c("time"))
nep_final <- nep_merge %>%
  dplyr::select(-time, -lon, -lat, -day, -month, -date, -week) %>%
  rename(lat = latitude,
         lon = longitude)

saveRDS(nep_final, file = "../data/nep_ipsl.Rdata")

nep_ipsl_mean <- tidync('F:/NEP_ROMS_ipsl/NEP_ipsl_variables_AvgMarchMay_1995-2100.nc')

nep_avg_var <- nep_ipsl_mean %>%
  hyper_tibble()

nep_avg_lat <- nep_ipsl_mean %>%
  activate("latitude") %>%
  hyper_tibble()

nep_avg_year <- nep_ipsl_mean %>%
  activate("years") %>%
  hyper_tibble()

nep_lat <- merge(nep_avg_var, nep_avg_lat, by = c("lat"))
nep_avgs <- merge(nep_lat, nep_avg_year, by = c("time"))
nep_avgs <- dplyr::select(nep_avgs, -time, -lat)

saveRDS(nep_avgs, file = "../data/nep_ipsl_means.Rdata")
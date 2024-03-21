### Title: sdmTMB Projecting
### Author: Rebecca Howard
### Date: 03/12/2024

# Load libraries ----
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB) 
library(visreg)
library(raster)
library(sf)
library(colorspace)
library(RColorBrewer)
library(mapdata)
library(fields)
library(ggpubr)
library(future)
library(scales)
library(tidync)
library(magick)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'projection_functions.R'))
source(here('code/functions', 'map_project.R'))

# Data
roms_ss <- readRDS(here('data', 'nep_ipsl.Rdata'))
roms_means <- readRDS(here('data', 'nep_ipsl_means.Rdata'))
hindcast_ss <- readRDS(here('data', 'nep_final.Rdata'))
hindcast_means <- readRDS(here('data', 'nep_avgs.Rdata'))

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

extra_years <- c(2020:2100)
base_dir <- getwd()

### Pacific Hake ---------------------------------------------------------------------------------------------------------------------------------
yoy_hake <- read_data('yoy_hake.Rdata') 
hake_mesh <- make_mesh(yoy_hake, 
                       xy_cols = c("X", "Y"),
                       cutoff = 18)

set.seed(1993)
hake_small <- sdmTMB(small ~ 0 + v_cu +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + v_cu,
                        extra_time = extra_years,
                        data = yoy_hake,
                        mesh = hake_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
hake_large <- sdmTMB(large ~ 0 + depth_iso26 +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + depth_iso26,
                        extra_time = extra_years,
                        data = yoy_hake,
                        mesh = hake_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
hake_hindcast <- sdm_cells(yoy_hake, hake_small, hake_large,
                           hindcast_means, hindcast_ss, 1995:2019)
saveRDS(hake_hindcast, file = here("data", "hake_hindcast.rds"))

# hake_hindcast <- readRDS(here("data", "hake_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_hindcast[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_hindcast[[2]], "Large (36-81 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl1 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(hake_ipsl1, file = here("data", "hake_ipsl1.rds"))

# hake_ipsl1 <- readRDS(here("data", "hake_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl1[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl1[[2]], "Large (36-81 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl2 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(hake_ipsl2, file = here("data", "hake_ipsl2.rds"))

# hake_ipsl2 <- readRDS(here("data", "hake_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl2[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl2[[2]], "Large (36-81 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl3 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(hake_ipsl3, file = here("data", "hake_ipsl3.rds"))

# hake_ipsl3 <- readRDS(here("data", "hake_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl3[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl3[[2]], "Large (36-81 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
hake_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_hake')
hake_imgs <- list.files(hake_dir_out, full.names = T)
hake_img_list <- lapply(hake_imgs, image_read)
hake_img_joined <- image_join(hake_img_list)
hake_img_animated <- image_animate(hake_img_joined, fps = 1)
image_write(image = hake_img_animated,
            path = here('results', 'forecast_output', 'yoy_hake', "hake_gifs.gif"))

# Remove objects
rm(hake_hindcast, hake_ipsl1, hake_ipsl2,
   hake_ipsl3, yoy_hake, hake_mesh, hake_large,
   hake_small, hake_dir_out, hake_imgs, hake_img_list,
   hake_img_joined, hake_img_animated)


# Northern Anchovy --------------------------------------------------------------------------------------------------------------------------------
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42, year > 2013)
anchovy_mesh <- make_mesh(yoy_anchovy,
                          xy_cols = c("X", "Y"),
                          cutoff = 15)

set.seed(1993)
anchovy_small <- sdmTMB(small ~ 0 + u_vint_100m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_100m,
                        extra_time = extra_years,
                        data = yoy_anchovy,
                        mesh = anchovy_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
anchovy_large <- sdmTMB(large ~ 0 + u_vint_100m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_100m,
                        extra_time = extra_years,
                        data = yoy_anchovy,
                        mesh = anchovy_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
anchovy_hindcast <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                              hindcast_means, hindcast_ss, 2014:2019)
saveRDS(anchovy_hindcast, file = here("data", "anchovy_hindcast.rds"))

# anchovy_hindcast <- readRDS(here("data", "anchovy_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_hindcast[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_hindcast[[2]], "Large (36-85 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl1 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(anchovy_ipsl1, file = here("data", "anchovy_ipsl1.rds"))

# anchovy_ipsl1 <- readRDS(here("data", "anchovy_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl1[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl1[[2]], "Large (36-85 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl2 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(anchovy_ipsl2, file = here("data", "anchovy_ipsl2.rds"))

# anchovy_ipsl2 <- readRDS(here("data", "anchovy_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl2[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl2[[2]], "Large (36-85 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl3 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(anchovy_ipsl3, file = here("data", "anchovy_ipsl3.rds"))

# anchovy_ipsl3 <- readRDS(here("data", "anchovy_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl3[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl3[[2]], "Large (36-81 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
anchovy_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_anchovy')
anchovy_imgs <- list.files(anchovy_dir_out, full.names = T)
anchovy_img_list <- lapply(anchovy_imgs, image_read)
anchovy_img_joined <- image_join(anchovy_img_list)
anchovy_img_animated <- image_animate(anchovy_img_joined, fps = 1)
image_write(image = anchovy_img_animated,
            path = here('results', 'forecast_output', 'yoy_anchovy', "anchovy_gifs.gif"))

# Remove objects
rm(anchovy_hindcast, anchovy_ipsl1, anchovy_ipsl2,
   anchovy_ipsl3, yoy_anchovy, anchovy_mesh, anchovy_large,
   anchovy_small, anchovy_dir_out, anchovy_imgs, anchovy_img_list,
   anchovy_img_joined, anchovy_img_animated)


# Pacific Sanddab ---------------------------------------------------------------------------------------------------------------------------------
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36, year > 2012)
sdab_mesh <- make_mesh(yoy_sdab,
                       xy_cols = c("X", "Y"),
                       cutoff = 15)

set.seed(1993)
sdab_small <- sdmTMB(small ~ 0 + u_vint_50m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_50m,
                        extra_time = extra_years,
                        data = yoy_sdab,
                        mesh = sdab_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
sdab_large <- sdmTMB(large ~ 0 + u_vint_50m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_50m,
                        extra_time = extra_years,
                        data = yoy_sdab,
                        mesh = sdab_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
sdab_hindcast <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                              hindcast_means, hindcast_ss, 2013:2019)
saveRDS(sdab_hindcast, file = here("data", "sdab_hindcast.rds"))

# sdab_hindcast <- readRDS(here("data", "sdab_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdab_map_project(sdab_hindcast[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
sdab_map_project(sdab_hindcast[[2]], "Large (26-55 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl1 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2020:2040)
saveRDS(sdab_ipsl1, file = here("data", "sdab_ipsl1.rds"))

# sdab_ipsl1 <- readRDS(here("data", "sdab_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdab_map_project(sdab_ipsl1[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
sdab_map_project(sdab_ipsl1[[2]], "Large (26-55 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl2 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2050:2070)
saveRDS(sdab_ipsl2, file = here("data", "sdab_ipsl2.rds"))

# sdab_ipsl2 <- readRDS(here("data", "sdab_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdab_map_project(sdab_ipsl2[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
sdab_map_project(sdab_ipsl2[[2]], "Large (26-55 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl3 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2080:2100)
saveRDS(sdab_ipsl3, file = here("data", "sdab_ipsl3.rds"))

# sdab_ipsl3 <- readRDS(here("data", "sdab_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdab_map_project(sdab_ipsl3[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
sdab_map_project(sdab_ipsl3[[2]], "Large (26-55 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
sdab_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_sanddab')
sdab_imgs <- list.files(sdab_dir_out, full.names = T)
sdab_img_list <- lapply(sdab_imgs, image_read)
sdab_img_joined <- image_join(sdab_img_list)
sdab_img_animated <- image_animate(sdab_img_joined, fps = 1)
image_write(image = sdab_img_animated,
            path = here('results', 'forecast_output', 'yoy_sanddab', "sdab_gifs.gif"))

# Remove objects
rm(sdab_hindcast, sdab_ipsl1, sdab_ipsl2,
   sdab_ipsl3, yoy_sdab, sdab_mesh, sdab_large,
   sdab_small, sdab_dir_out, sdab_imgs, sdab_img_list,
   sdab_img_joined, sdab_img_animated)


# Shortbelly Rockfish ---------------------------------------------------------------------------------------------------------------------------------
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
shortbelly_mesh <- make_mesh(yoy_shortbelly,
                             xy_cols = c("X", "Y"),
                             cutoff = 18)

set.seed(1993)
shortbelly_small <- sdmTMB(small ~ 0 + depth_iso26 +
                             s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             s(sss_scaled, k = 3),
                           spatial_varying = ~ 0 + depth_iso26,
                           extra_time = extra_years,
                           data = yoy_shortbelly,
                           mesh = shortbelly_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
set.seed(1993)
shortbelly_large <- sdmTMB(large ~ 0 + vgeo +
                             s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             s(sss_scaled, k = 3),
                           spatial_varying = ~ 0 + vgeo,
                           extra_time = extra_years,
                           data = yoy_shortbelly,
                           mesh = shortbelly_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
shortbelly_hindcast <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                                 hindcast_means, hindcast_ss, 1995:2019)
saveRDS(shortbelly_hindcast, file = here("data", "shortbelly_hindcast.rds"))

# shortbelly_hindcast <- readRDS(here("data", "shortbelly_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_hindcast[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_hindcast[[2]], "Large (36-78 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl1 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2020:2040)
saveRDS(shortbelly_ipsl1, file = here("data", "shortbelly_ipsl1.rds"))

# shortbelly_ipsl1 <- readRDS(here("data", "shortbelly_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl1[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl1[[2]], "Large (36-78 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl2 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2050:2070)
saveRDS(shortbelly_ipsl2, file = here("data", "shortbelly_ipsl2.rds"))

# shortbelly_ipsl2 <- readRDS(here("data", "shortbelly_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl2[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl2[[2]], "Large (36-78 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl3 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2080:2100)
saveRDS(shortbelly_ipsl3, file = here("data", "shortbelly_ipsl3.rds"))

# shortbelly_ipsl3 <- readRDS(here("data", "shortbelly_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl3[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl3[[2]], "Large (36-78 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
shortbelly_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_shortbelly')
shortbelly_imgs <- list.files(shortbelly_dir_out, full.names = T)
shortbelly_img_list <- lapply(shortbelly_imgs, image_read)
shortbelly_img_joined <- image_join(shortbelly_img_list)
shortbelly_img_animated <- image_animate(shortbelly_img_joined, fps = 1)
image_write(image = shortbelly_img_animated,
            path = here('results', 'forecast_output', 'yoy_shortbelly', "shortbelly_gifs.gif"))

# Remove objects
rm(shortbelly_hindcast, shortbelly_ipsl1, shortbelly_ipsl2,
   shortbelly_ipsl3, yoy_shortbelly, shortbelly_mesh, shortbelly_large,
   shortbelly_small, shortbelly_dir_out, shortbelly_imgs, shortbelly_img_list,
   shortbelly_img_joined, shortbelly_img_animated)


# Widow Rockfish ---------------------------------------------------------------------------------------------------------------------------------
yoy_widow <- read_data('yoy_widw.Rdata') 
widow_mesh <- make_mesh(yoy_widow,
                        xy_cols = c("X", "Y"),
                        cutoff = 18)

set.seed(1993)
widow_small <- sdmTMB(small ~ 0 + vmax_cu +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + vmax_cu,
                        extra_time = extra_years,
                        data = yoy_widow,
                        mesh = widow_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
widow_large <- sdmTMB(large ~ 0 + spice_iso26 +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + spice_iso26,
                        extra_time = extra_years,
                        data = yoy_widow,
                        mesh = widow_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
widow_hindcast <- sdm_cells(yoy_widow, widow_small, widow_large,
                           hindcast_means, hindcast_ss, 1995:2019)
saveRDS(widow_hindcast, file = here("data", "widow_hindcast.rds"))

# widow_hindcast <- readRDS(here("data", "widow_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_hindcast[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_hindcast[[2]], "Large (33-64 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl1 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(widow_ipsl1, file = here("data", "widow_ipsl1.rds"))

# widow_ipsl1 <- readRDS(here("data", "widow_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl1[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl1[[2]], "Large (33-64 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl2 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(widow_ipsl2, file = here("data", "widow_ipsl2.rds"))

# widow_ipsl2 <- readRDS(here("data", "widow_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl2[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl2[[2]], "Large (33-64 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl3 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(widow_ipsl3, file = here("data", "widow_ipsl3.rds"))

# widow_ipsl3 <- readRDS(here("data", "widow_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl3[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl3[[2]], "Large (33-64 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
widow_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_widow')
widow_imgs <- list.files(widow_dir_out, full.names = T)
widow_img_list <- lapply(widow_imgs, image_read)
widow_img_joined <- image_join(widow_img_list)
widow_img_animated <- image_animate(widow_img_joined, fps = 1)
image_write(image = widow_img_animated,
            path = here('results', 'forecast_output', 'yoy_widow', "widow_gifs.gif"))

# Remove objects
rm(widow_hindcast, widow_ipsl1, widow_ipsl2,
   widow_ipsl3, yoy_widow, widow_mesh, widow_large,
   widow_small, widow_dir_out, widow_imgs, widow_img_list,
   widow_img_joined, widow_img_animated)


# Market Squid ------------------------------------------------------------------------------------------------------------------------------------
yoy_squid <- read_data('yoy_squid.Rdata')
squid_mesh <- make_mesh(yoy_squid,
                        xy_cols = c("X", "Y"),
                        cutoff = 15)

set.seed(1993)
squid_small <- sdmTMB(small ~ 0 + depth_iso26 +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + depth_iso26,
                      extra_time = extra_years,
                      data = yoy_squid,
                      mesh = squid_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
set.seed(1993)
squid_large <- sdmTMB(large ~ 0 + spice_iso26 +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + spice_iso26,
                      extra_time = extra_years,
                      data = yoy_squid,
                      mesh = squid_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
squid_hindcast <- sdm_cells(yoy_squid, squid_small, squid_large,
                                 hindcast_means, hindcast_ss, 2004:2019)
saveRDS(squid_hindcast, file = here("data", "squid_hindcast.rds"))

# squid_hindcast <- readRDS(here("data", "squid_hindcast.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_hindcast[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_hindcast[[2]], "Large (71-139 mm)", "")
mtext(substitute(paste(bold("Hindcast"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl1 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2020:2040)
saveRDS(squid_ipsl1, file = here("data", "squid_ipsl1.rds"))

# squid_ipsl1 <- readRDS(here("data", "squid_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl1[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl1[[2]], "Large (71-139 mm)", "")
mtext(substitute(paste(bold("2020-2040"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl2 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2050:2070)
saveRDS(squid_ipsl2, file = here("data", "squid_ipsl2.rds"))

# squid_ipsl2 <- readRDS(here("data", "squid_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl2[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl2[[2]], "Large (71-139 mm)", "")
mtext(substitute(paste(bold("2050-2070"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl3 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2080:2100)
saveRDS(squid_ipsl3, file = here("data", "squid_ipsl3.rds"))

# squid_ipsl3 <- readRDS(here("data", "squid_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 3.5, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl3[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl3[[2]], "Large (71-139 mm)", "")
mtext(substitute(paste(bold("2080-2100"))), 
      line = 0, side = 3, outer = TRUE, cex = 4, family = "serif")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl3[[3]]))

##### GIFs -------------------------------------------------------------------------------------------------------------------------
squid_dir_out <- file.path(base_dir, 'results', 'forecast_output', 'yoy_squid')
squid_imgs <- list.files(squid_dir_out, full.names = T)
squid_img_list <- lapply(squid_imgs, image_read)
squid_img_joined <- image_join(squid_img_list)
squid_img_animated <- image_animate(squid_img_joined, fps = 1)
image_write(image = squid_img_animated,
            path = here('results', 'forecast_output', 'yoy_squid', "squid_gifs.gif"))

# Remove objects
rm(squid_hindcast, squid_ipsl1, squid_ipsl2,
   squid_ipsl3, yoy_squid, squid_mesh, squid_large,
   squid_small, squid_dir_out, squid_imgs, squid_img_list,
   squid_img_joined, squid_img_animated)

### Title: SVC term plots
### Author: Rebecca Howard
### Date: 03/27/2024

# Load libraries ----
library(here)
library(ggplot2)
library(dplyr)
library(colorspace)
library(RColorBrewer)
library(viridis)


# Data ----
roms_means <- readRDS(here('data', 'nep_ipsl_means.Rdata'))
hindcast_means <- readRDS(here('data', 'nep_avgs.Rdata'))

# Functions ----
plot_term1 <- function(roms, variable, ylab){
  ggplot(roms, aes(x = latitude, 
                   y = variable)) +
    geom_line(aes(group = years,
                  color = years),
              linewidth = 1.3) +
    scale_color_viridis(discrete = FALSE,
                        option = "B",
                        end = .9)  +
    labs(x = "Latitude \u00B0N",
         y = ylab) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 14),
          axis.title = element_text(family = "serif", size = 18))
}

plot_term2 <- function(roms, variable, ylab){
  ggplot(roms, aes(x = years, 
                   y = variable)) +
    geom_line(aes(group = latitude,
                  color = latitude),
              linewidth = 1.3) +
    scale_color_viridis(discrete = FALSE,
                        option = "B",
                        end = .9)  +
    labs(x = "Years",
         y = ylab) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 14),
          axis.title = element_text(family = "serif", size = 18))
}

# Plot SVC variables ----
# Latitude y-axis
plot_term1(hindcast_means, hindcast_means$spice_iso26, "Spiciness")
plot_term1(hindcast_means, hindcast_means$v_cu, "CU Mean Velocity")
plot_term1(hindcast_means, hindcast_means$vmax_cu, "CU Maximum Velocity")
plot_term1(hindcast_means, hindcast_means$depth_iso26, "26 kg/m\u00B3 Isopycnal Depth")
plot_term1(hindcast_means, hindcast_means$u_vint_100m, "Eastward u vertically integrated 0-100m")

# Year y-axis
plot_term2(hindcast_means, hindcast_means$spice_iso26, "Spiciness")
plot_term2(hindcast_means, hindcast_means$v_cu, "CU Mean Velocity")
plot_term2(hindcast_means, hindcast_means$vmax_cu, "CU Maximum Velocity")
plot_term2(hindcast_means, hindcast_means$depth_iso26, "26 kg/m\u00B3 Isopycnal Depth")
plot_term2(hindcast_means, hindcast_means$u_vint_100m, "Eastward u vertically integrated 0-100m")

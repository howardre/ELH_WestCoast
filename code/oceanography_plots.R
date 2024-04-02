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
library(ggpubr)

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
                        end = .9,
                        name = "Year")  +
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
                        end = .9,
                        name = "Latitude \u00B0N")  +
    labs(x = "Years",
         y = ylab) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(6, 6, 6, 6), "mm"),
          plot.title = element_text(size = 30, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 22),
          axis.title = element_text(family = "serif", size = 26),
          legend.text = element_text(family = "serif", size = 23),
          legend.title = element_text(family = "serif", size = 26))
}

# Plot SVC variables ----
# Latitude y-axis
plot_term1(hindcast_means, hindcast_means$spice_iso26, "Spiciness")
plot_term1(hindcast_means, hindcast_means$v_cu, "CU Mean Velocity")
plot_term1(hindcast_means, hindcast_means$vmax_cu, "CU Maximum Velocity")
plot_term1(hindcast_means, hindcast_means$depth_iso26, "26 kg/m\u00B3 Isopycnal Depth")
plot_term1(hindcast_means, hindcast_means$u_vint_100m, "Eastward u vertically integrated 0-100m")

# Year y-axis
p1 <- plot_term2(hindcast_means, hindcast_means$spice_iso26, "Spiciness")
p2 <- plot_term2(hindcast_means, hindcast_means$v_cu, "CU Mean Velocity")
p3 <- plot_term2(hindcast_means, hindcast_means$vmax_cu, "CU Maximum Velocity")
p4 <- plot_term2(hindcast_means, hindcast_means$depth_iso26, "26 kg/m\u00B3 Isopycnal Depth")
p5 <- plot_term2(hindcast_means, hindcast_means$u_vint_100m, "Eastward u vertically integrated 0-100m")
p6 <- plot_term2(hindcast_means, hindcast_means$vgeo, "Geostrophic Current")

windows(height = 15, 
        width = 22)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.copy(jpeg, here('results', 
                    'variable_plots.jpg'), 
         height = 15, 
         width = 22, 
         units = 'in', 
         res = 200)
dev.off()

# Future values
p7 <- plot_term2(roms_means, roms_means$spice_iso26, "Spiciness")
p8 <- plot_term2(roms_means, roms_means$v_cu, "CU Mean Velocity")
p9 <- plot_term2(roms_means, roms_means$vmax_cu, "CU Maximum Velocity")
p10 <- plot_term2(roms_means, roms_means$depth_iso26, "26 kg/m\u00B3 Isopycnal Depth")
p11 <- plot_term2(roms_means, roms_means$u_vint_100m, "Eastward u vertically integrated 0-100m")
p12 <- plot_term2(roms_means, roms_means$vgeo, "Geostrophic Current")

windows(height = 15, 
        width = 22)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
ggarrange(p7, p8, p9, p10, p11, p12,
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.copy(jpeg, here('results', 
                    'forecast_plots.jpg'), 
         height = 15, 
         width = 22, 
         units = 'in', 
         res = 200)
dev.off()

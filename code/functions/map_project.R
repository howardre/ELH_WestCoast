# Function to make projection maps
map_project <- function(grid, title, latitude_label){
  nlat = 80
  nlon = 120
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(rev(c("#FFFFCC", "#FBF2A8", "#F9E585",
                                    "#F5D363", "#EFBA55", "#EAA352",
                                    "#E68C51", "#E0754F", "#D75C4D",
                                    "#BB4A48", "#994240", "#763931", 
                                    "#542D20", "#352311", "#191900")))
  image(lond,
        latd,
        t(matrix(grid$pred_scaled,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-126, -116),
        ylim = range(grid$lat, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(grid$pred_scaled,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(100), 
        ylab = latitude_label,
        xlab = "Longitude \u00B0W",
        xlim = c(-126, -116),
        ylim = range(grid$lat, na.rm = TRUE) + c(-.4, .5),
        zlim = c(min(grid$pred_scaled, na.rm = T), 
                 max(grid$pred_scaled, na.rm = T)),
        main = title,
        cex.lab = 2.9,
        cex.axis = 2.3,
        cex.main = 3)
  maps::map("state",
            boundary = FALSE,
            fill = TRUE,
            col = "wheat4",
            add = TRUE)
  text(x = state_labels$lon, 
       y = state_labels$lat,
       state_labels$name, 
       pos = 1,
       col = "black",
       cex = 2.6,
       family = "serif")
  image.plot(legend.only = T,
             col = my_color(100),
             legend.shrink = 0.2,
             smallplot = c(.33, .38, .11, .24),
             legend.cex = 1.5,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(grid$pred_scaled, na.rm = T), 
                      max(grid$pred_scaled, na.rm = T)),
             legend.args = list("Scaled \n Abundance",
                                side = 2, 
                                cex = 2.2,
                                family = "serif",
                                line = 1))
}


svc_hindcast <- function(grid, title, latitude_label, legend){
  nlat = 80
  nlon = 120
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))
  image(lond,
        latd,
        t(matrix(grid$avg_zeta,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-126, -116),
        ylim = range(grid$lat, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(grid$avg_zeta,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(100), 
        ylab = latitude_label,
        xlab = "Longitude \u00B0W",
        xlim = c(-126, -116),
        ylim = range(grid$lat, na.rm = TRUE) + c(-.4, .5),
        zlim = c(min(grid$avg_zeta, na.rm = T), 
                 max(grid$avg_zeta, na.rm = T)),
        main = title,
        cex.lab = 2.9,
        cex.axis = 2.3,
        cex.main = 3)
  maps::map("state",
            boundary = FALSE,
            fill = TRUE,
            col = "wheat4",
            add = TRUE)
  text(x = state_labels$lon, 
       y = state_labels$lat,
       state_labels$name, 
       pos = 1,
       col = "black",
       cex = 2.6,
       family = "serif")
  image.plot(legend.only = T,
             col = my_color(100),
             legend.shrink = 0.2,
             smallplot = c(.33, .38, .11, .24),
             legend.cex = 1.5,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(grid$avg_zeta, na.rm = T), 
                      max(grid$avg_zeta, na.rm = T)),
             legend.args = list(legend,
                                side = 2, 
                                cex = 2.2,
                                family = "serif",
                                line = 1))
}


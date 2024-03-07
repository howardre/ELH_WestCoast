sdmTMB_map <- function(df, preds, title, lat_label){
  
  my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  color_levels = 100
  max_absolute_value = max(abs(c(min(preds$est, na.rm = T),
                                 max(preds$est, na.rm = T))))
  color_sequence = seq(max(preds$est, na.rm = T), 
                       min(preds$est, na.rm = T),
                       length.out = color_levels + 1)
  n_in_class = hist(preds$est, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  
  # Make map
  image(lond,
        latd,
        t(matrix(exp(preds$est),
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-126, -116),
        ylim = range(df$lat, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(preds$est,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include],
        main = title,
        ylab = lat_label,
        xlab = "Longitude \u00B0W",
        xlim = c(-126, -116),
        ylim = range(df$lat, na.rm = TRUE) + c(-.4, .5),
        cex.lab = 3.1,
        cex.axis = 2.3,
        cex.main = 3.4)
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
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.3, .35, .11, .24),
             legend.cex = 1.5,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(preds$est, na.rm = T), 
                      max(preds$est, na.rm = T)),
             legend.args = list("log(CPUE + 1)",
                                side = 2, 
                                cex = 2.2,
                                family = "serif",
                                line = 1))
}

sdmTMB_SVC <- function(df, preds, title, lat_label, var, variable){

  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))
  color_levels = 100
  max_absolute_value = max(abs(c(min(var, na.rm = T),
                                 max(var, na.rm = T))))
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(var, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  
  # Make map
  image(lond,
        latd,
        t(matrix(var,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-126, -116),
        ylim = range(df$lat, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(var,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include],
        main = title,
        ylab = lat_label,
        xlab = "Longitude \u00B0W",
        xlim = c(-126, -116),
        ylim = range(df$lat, na.rm = TRUE) + c(-.4, .5),
        cex.lab = 3.1,
        cex.axis = 2.3,
        cex.main = 3.4)
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
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.3, .35, .11, .24),
             legend.cex = 1.5,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(var, na.rm = T), 
                      max(var, na.rm = T)),
             legend.args = list(variable,
                                side = 2, 
                                cex = 1.8,
                                family = "serif",
                                line = 1))
}

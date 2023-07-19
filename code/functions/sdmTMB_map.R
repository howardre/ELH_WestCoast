sdmTMB_map <- function(df, preds){
  
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
  windows(width = 7, height = 10)
  par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0),
      family = "serif")
  image(lond,
        latd,
        t(matrix(preds$est,
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
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-126, -116),
        ylim = range(df$lat, na.rm = TRUE) + c(-.4, .5),
        cex.main = 2,
        cex.lab = 2,
        cex.axis = 1.8)
  map("worldHires",
      fill = T,
      col = "wheat4",
      add = T)
  image.plot(legend.only = T,
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.18, .21, .17, .38),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.6,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(preds$est, na.rm = T), 
                      max(preds$est, na.rm = T)),
             legend.args = list("CPUE",
                                side = 2,
                                cex = 1.8,
                                line = 1.3,
                                family =  "serif"))
  
}

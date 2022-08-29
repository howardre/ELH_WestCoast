### Libraries ----
library(rerddap)

### Bering 10K model output ----
# Using avg surface temperatures & salinity
RREAS_CTD_info <- info('FED_Rockfish_CTD')

ctds <- tabledap(RREAS_CTD_info,
                 fields = c('temperature', 'ctd_depth', 'time'),
                 'ctd_depth<=10')

cols <- c('temperature')
ctds[cols] <- sapply(ctds[cols], as.numeric)

ctd_means <- ctds  %>%
  group_by(time) %>%
  summarise_at(vars('temperature'), mean) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(date = format(as.POSIXct(time, tz = "UTC", format = "%Y-%m-%d")),
         year = lubridate::year(date),
         month = lubridate::month(date)) %>%
  select(temperature, year, month) %>%
  group_by(year) %>%
  summarise_at(vars('temperature'), mean, na.rm = T)


saveRDS(ctd_means, file = here('data', 'RREAS_ctd_means.rds'))

# Check to make sure logical - use code in function to get the temp_filtered df
# temp_filtered %>%
#   slice(1:10000) %>%
# ggplot(aes(x = lon, y = lat, color = temp)) +
#   geom_point(size = 5, alpha = 0.5) +
#   scale_color_viridis_b()

# Plot the mean temps over year
ggplot(data = ctd_means) +
  geom_point(aes(x = year, y = temperature), color = "aquamarine4", size = 3) +
  labs(title = "Mean CTD Surface Temperature",
       y = "Temperature (C)",
       x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 24),
        text = element_text(family = "serif"))

dev.copy(jpeg, here('results', 'mean_RREAS_ctd_temps.jpeg'), 
         height = 10, width = 10, units = 'in', res = 200)
dev.off()

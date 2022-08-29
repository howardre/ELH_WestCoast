# Libraries 
library(rerddap)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(here)

# Using avg surface temperatures
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
         year = lubridate::year(date)) %>%
  select(temperature, year) %>%
  group_by(year) %>%
  summarise_at(vars('temperature'), mean, na.rm = T)

saveRDS(ctd_means, file = here('data', 'RREAS_ctd_means.rds'))

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

### Exploration of hake and sanddab data
### Date: 12/08/2021
### Author: Rebecca Howard

### Load libraries
library(readxl)
library(plyr)
library(tibble)
library(here)
library(ggplot2)
library(lubridate)

### Functions
# Get data into appropriate format
clean_data <- function(data_list, i){
  data_frame <- as.data.frame(data_list[[i]])
  data_frame <- mutate(data_frame,
                     year = lubridate::year(haul_date),
                     month = lubridate::month(haul_date),
                     day = lubridate::day(haul_date))
  data_frame$doy <- as.numeric(mdy.date(data_frame$month, data_frame$day, 1960))
  data_frame$lat <- data_frame$net_in_lat * 0.01
  data_frame$lon <- data_frame$net_in_lon * -0.01
  return(data_frame)
}

# Plot relevant variables
plot_variables <- function(data){
  plot(table(data$year[data$catch > 0]),
       ylab = 'Frequency',
       xlab = 'Year',
       main = '')
  plot(table(data$doy[data$catch > 0]),
       ylab = 'Frequency',
       xlab = 'Day of Year',
       main = '')
  plot(table(data$lat[data$catch > 0]),
       ylab = 'Frequency',
       xlab = 'Latitude',
       main = '')
  plot(table(data$lon[data$catch > 0]),
       ylab = 'Frequency',
       xlab = 'Longitude',
       main = '')
  plot(table(data$strata[data$catch > 0]),
       ylab = 'Frequency',
       xlab = 'Strata',
       main = '')
  hist(data$catch, 
       breaks = 100,
       main = '',
       xlab = 'Catch') 
}

### Load data
exl_data <- here('data', 'hake_dabs.xlsx')
excel_sheets(path = exl_data) # view the worksheets
tab_names <- excel_sheets(path = exl_data)

data_list <- lapply(tab_names, function(x) read_excel(path = exl_data, sheet = x))

### Clean data
yoy_hake <- clean_data(data_list, 1)
yoy_dab <- clean_data(data_list, 2)
adult_dab <- clean_data(data_list, 3)

### Basic information
# Hake
par(mfrow = c(2, 3))
plot_variables(yoy_hake)

# Sanddab
par(mfrow = c(2, 3))
plot_variables(yoy_dab)

par(mfrow = c(2, 3))
plot_variables(adult_dab)


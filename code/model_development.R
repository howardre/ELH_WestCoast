### Title: Model Development
### Author: Rebecca Howard
### Date: 09/13/2022

# Libraries
library(mgcv)
library(maps)
library(mapdata)
library(here)

# Load data
RREAS_list <- readRDS(here('data', 'RREAS_list.rdata'))
lengths <- readRDS(here('data', 'length_list.rdata'))

# Functions

merge_data <- function(data1, data2){
  merge(data1, data2,
        by.x = c("cruise", "haul_no"),
        by.y = c("cruise", "haul_no"),
        all.x = F, all.y = T)
}

length_split <- function(length_data, split_value, species){
  lower <- aggregate(std_length ~ haul_no + cruise,
                     data = length_data,
                     FUN = length,
                     subset = std_length <= split_value) # less than or equal for the smaller fish
  yoy <- merge_data(lower, species)
  names(yoy)[names(yoy) == "std_length"] <- "lower_num"
  
  upper <- aggregate(std_length ~ haul_no + cruise,
                     data = length_data,
                     FUN = length,
                     subset = std_length > split_value)
  yoy <- merge_data(upper, yoy)
  names(yoy)[names(yoy) == "std_length"] <- "upper_num"
  
  total <- aggregate(std_length ~ haul_no + cruise,
                     data = length_data,
                     FUN = length)
  yoy <- merge_data(total, yoy)
  names(yoy)[names(yoy) == "std_length"] <- "total_num"
  
  yoy$upper_cpue <- (yoy$upper_num / yoy$total_num) * yoy$catch
  yoy$lower_cpue <- (yoy$lower_num / yoy$total_num) * yoy$catch
  
  yoy$upper_cpue[is.na(yoy$upper_cpue)] <- 0
  yoy$lower_cpue[is.na(yoy$lower_cpue)] <- 0
  
  yoy_filtered <- yoy[yoy$year > min(length_data$year), ]
  yoy_filtered$lncpue <- log(yoy_filtered$catch + 1)
  
  return(yoy_filtered)
}

# Hake ----
# Split into different life stages (<= 31 larvae, >31 mm juveniles - this happens to be the median)
# Using Phillips et al. (2007) to determine length to split at
hake_data <- length_split(lengths[[1]], 31, RREAS_list[[1]])

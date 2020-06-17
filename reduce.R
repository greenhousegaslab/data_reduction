#-----------------------------------------------------------------------------#
# SCRIPT: reduce.R                                                            #
# PURPOSE: Reduce large atmospheric satellite datasets.                       #
# Xiaoling Liu, June 11, 2019                                                 #
#-----------------------------------------------------------------------------#

#-----------------------------------#
# Clear working environment         #
#-----------------------------------#
rm(list=ls())

#-----------------------------------#
# Initialization                    #
#-----------------------------------#
library(R.matlab)
library(geosphere)
library(dplyr)
source("CalculateRange.R")
# File that saves the observations sounding id, extracted from satellite 
# observation information file (information.mat).
# Reason: R.matlab package can't import souding id from mat file correctly.
# Details for data format:
#   sounding_id (data.frame)
#     V1 (num)
load("sounding_id.Rdata")

#-----------------------------------#
# Constant Variables                #
#-----------------------------------#
#! CHANGE THE PARAMETERS BELOW FOR YOUR PARTICULAR SETUP
# File path to read in satellite observation information
# Details for data format:
#   outmat (cell)
#     outmat{1,1} = sounding_id (uint64)
#     outmat{2,1} = lon (single)
#     outmat{3,1} = lat (single)
#     outmat{4,1} = date (int16)
#     outmat{5,1} = orbit number (int32)
#     outmat{6,1} = xco2 (single)
#     outmat{7,1} = quality flag (int8)
#     outmat{8,1} = land or ocean (int8)
#     outmat{9,1} = sounding path number (uint8)
kFilePathInfo <- "information.mat"
# File path to read in background, i.e. background CO2 levels
# Details for data format:
#   background (double)
kFilePathBackground <- "background.mat"
# Level of reduction: how much we want to reduce the data
# i.e. range = 500km, kScalingFactor = 0.05
# Meaning that we want to choose the next observation at the distance of 
# 25km (500 * 0.05 = 25km)
kScalingFactor <- 0.05
# Acceptable deviation of distance when choosing the next observation
# i.e. range = 500km, kScalingFactor = 0.05, kDevFactor = 0.2
# Meaning that we allow to choose the next observation at the distance of 
# 30km (500 * 0.05 * (1 + 0.2) = 30km), with 5km deviation (500 * 0.05 * 0.2 = 5km).
kDevFactor <- 0.2
# File path to save result
# Details for data format:
#   z_compress_unique (data.frame)
#     sounding_id (num)
#     lon (num)
#     lat (num)
#     year (num)
#     month (num)
#     day (num)
#     xco2 (num)
#     quality_flag (num)
#     land_ocean (num)
#     sounding_path_number (num)
#     bg (num)
#     co2 (num)
#     co2_est (num)
#     range (num)
#     sill (num)
#     est_error (num)
kFileSaveZCompress <- "z_compress.Rdata"
# File path to save variance result
# Details for data format:
#   var_result (data.frame)
#     var_result (num)
kFileSaveVarResult <- "var_result.Rdata"

#-----------------------------------#
# Data pre-process                  #
#-----------------------------------#
bg <- readMat(kFilePathBackground)
bg <- bg[["background"]]
oco2info <- readMat(kFilePathInfo)
lon <- oco2info[['outmat']][[2]][[1]]
lat <- oco2info[['outmat']][[3]][[1]]
date <- oco2info[['outmat']][[4]][[1]]
orbit_number <- oco2info[['outmat']][[5]][[1]]
xco2 <- oco2info[['outmat']][[6]][[1]]
quality_flag <- oco2info[['outmat']][[7]][[1]]
land_ocean <- oco2info[['outmat']][[8]][[1]]
sounding_path_number <- oco2info[['outmat']][[9]][[1]]
oco2info_matrix <- cbind(sounding_id, lon, lat, date,orbit_number, xco2, 
                         quality_flag, land_ocean, sounding_path_number, bg)
# extract information needed
oco2info_need <- oco2info_matrix[, c(1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 16)]
colnames(oco2info_need) <- c('sounding_id', 'lon', 'lat', 'year', 'month', 
                             'day', 'xco2', 'quality_flag', 'land_ocean', 
                             'sounding_path_number', 'bg')
oco2info_need <- data.frame(oco2info_need)
# remove duplicate observations and average the background
oco2info_need <- oco2info_need %>% group_by(sounding_id) 
                 %>% dplyr::summarize_all(mean, na.rm = TRUE)
# remove background
oco2info_need$co2 <- oco2info_need$xco2 - oco2info_need$bg
# remove ocean data
oco2info_need <- subset(oco2info_need, land_ocean == 1)
# remove bad quality data
oco2info_need <- subset(oco2info_need, quality_flag == 0)
oco2info_need <- data.frame(oco2info_need)

#-----------------------------------#
# Data reduction                  #
#-----------------------------------#
path <- unique(oco2info_need[['sounding_path_number']])
result = list()
var_result = list()
for (path_index in 1:length(path)) {
  # going through each path and compress its data
  path_number = path[path_index]
  # extract observations from the same path
  oco2info_group <- subset(oco2info_need, 
                           oco2info_need['sounding_path_number'] == path_number)
  # sort observations from north to south
  oco2info_group_sort <- oco2info_group[order(oco2info_group['lat'], decreasing = T), ]
  # choose the northest observation as base point
  base_point <- oco2info_group_sort[1, ]
  # calculate correlation length (also called range) for base point
  result_ok <- CalculateRange(base_point, oco2info_need)
  # (result_ok: range, co2_est, sill, est_error)
  range <- result_ok[[1]]
  # add the first point into result
  result[[length(result) + 1]] <- cbind(base_point, result_ok[[2]], result_ok[[1]], 
                                        result_ok[[3]], result_ok[[4]])
  # calculate range for the southest observation
  last_point <- oco2info_group_sort[nrow(oco2info_group_sort), ]
  result_ok <- CalculateRange(last_point, oco2info_need)
  # add the southest observation into result
  result[[length(result) + 1]] <- cbind(last_point, result_ok[[2]], result_ok[[1]], 
                                        result_ok[[3]], result_ok[[4]])
  # calculate the distance for selecting next observation
  base_dist <- kScalingFactor * range
  # set the acceptable deviation for distance
  dev_dist <- kDevFactor * base_dist
  for (point_index in 1:nrow(oco2info_group_sort)){
    # going through points on the path and decide whether to take it
    # current_dist: unit km, distm: unit m
    current_dist <- 0.001 * distm(c(base_point$lon, base_point$lat),
                                  c(oco2info_group_sort[point_index, ]$lon,
                                    oco2info_group_sort[point_index, ]$lat))
    if (current_dist >= base_dist) {
      if ((current_dist - base_dist < base_dist - pre_dist) & 
          (current_dist - base_dist <= dev_dist)) {
        # calculate variance for observations being removed
        points_removed <- subset(oco2info_group_sort, 
                                 oco2info_group_sort['lat'] < base_point$lat)
        points_removed <- subset(points_removed, 
                                 points_removed['lat'] > 
                                 oco2info_group_sort[point_index, ]$lat)
        # add variance into var_result
        var_result[[length(var_result) + 1]] <- var(points_removed$co2)
        # select a new base point
        base_point <- oco2info_group_sort[point_index, ]
        pre_dist <- 0
        # calculate range for new base point
        result_ok <- CalculateRange(base_point, oco2info_need)
        range <- result_ok[[1]]
        base_dist <- kScalingFactor * range
        # add the new base point into result
        result[[length(result) + 1]] <- cbind(base_point, result_ok[[2]], 
                                              result_ok[[1]], result_ok[[3]], 
                                              result_ok[[4]])
      } else {
        # calculate variance for observations being removed
        points_removed <- subset(oco2info_group_sort, 
                                 oco2info_group_sort['lat'] < base_point$lat)
        points_removed <- subset(points_removed, 
                                 points_removed['lat'] > 
                                 oco2info_group_sort[point_index-1, ]$lat)
        # add variance into var_result
        var_result[[length(var_result) + 1]] <- var(points_remo21ved$co2)
        # select a new base point
        base_point <- oco2info_group_sort[point_index-1, ]
        pre_dist <- 0.001 * distm(c(base_point$lon, base_point$lat),
                                  c(oco2info_group_sort[point_index, ]$lon,
                                    oco2info_group_sort[point_index, ]$lat))
        # calculate range for new base point
        result_ok <- CalculateRange(base_point, oco2info_need)
        range <- result_ok[[1]]
        base_dist <- kScalingFactor * range
        # add the new base point into result
        result[[length(result) + 1]] <- cbind(base_point, result_ok[[2]], 
                                              result_ok[[1]], result_ok[[3]], 
                                              result_ok[[4]])
      }
    } else {
      pre_dist <- current_dist
    }
  } ## end of point_index loop
} ## end of path_index loop

#-----------------------------------#
# Save results                      #
#-----------------------------------#
var_result <- do.call('rbind', var_result)
var_result <- data.frame(var_result)
save(var_result, file = kFileSaveVarResult)

z_compress <- do.call('rbind', result)
z_compress <- data.frame(z_compress)
colnames(z_compress)[13] <- 'co2_est'
colnames(z_compress)[14] <- 'range'
colnames(z_compress)[15] <- 'sill'
colnames(z_compress)[16] <- 'est_error'
# remove duplicated points
z_compress_unique <- unique(z_compress, by = 'sounding_id')
save(z_compress_unique, file = kFileSaveZCompress)

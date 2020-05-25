#-----------------------------------------------------------------------------#
# SCRIPT: CalculateRange.R                                                    #
# PURPOSE: Calculate correlation length and use local kriging                 #
#          to estimate co2.                                                   #
# Xiaoling Liu, June 11, 2019                                                #
#-----------------------------------------------------------------------------#

CalculateRange <- function(base_point, oco2info_need){
  # Compute sill, range, estimated CO2, estimated error
  #
  # Args:
  #   base_point: current observation being selected
  #   oco2info_need: all observations filtered by data pre-process
  #
  # Returns:
  #   range, estimated CO2, sill, estimated error

  #-----------------------------------#
  # Initialization                    #
  #-----------------------------------#
  library(geoR)
  library(rgdal)

  #-----------------------------------#
  # Constant Variables                #
  #-----------------------------------#
  # Choose observations within a circle with radius kDistance
  # Observations chosen are used for variogram
  # (unit: m)
  kDistance <- 2000000
  # Reference point at one of the edges of your domain 
  # Used for coordinates conversion
  kRefLon <- -180
  kRefLat <- 9
  # Parameters for plotting experimental variogram
  # (unit for kMaxDist: km)
  # Note: kMaxDist should be equal to or smaller than kDistance.
  kBin <- 6
  kMaxDist <- 2000
  # Initial parameters for fitting variogram
  kInitialValues <- expand.grid(seq(0, 100, by = 5), seq(100, 135, by = 5))
  kCovModel = 'exp'
  kNuggetLand <- (1.3)^2
  kWeights = 'npairs'
  # If range < 10km, re-set bin and maximum distance for experimental variogram
  # Reason: we found that most of the fitting variograms with range smaller than 
  #         10km had bad fit, and so we re-set bin and max distance and re-fit.
  kReFitThreshold <- 10
  # Re-set bin and maximum distance for experimental variogram
  kBinAgain <- 5
  kMaxDistAgain <- 200
  # Threshold for range
  kMaxRange <- 1000
  kMinRange <- 20
  # Parameters for kriging
  kKriging <- 'ok'
  kKrigingCovModel <- 'exponential'

  #-----------------------------------#
  # Data pre-process                  #
  #-----------------------------------#
  # Note for if statement: we want to choose observations within 2 days. 
  # The time period of data for this script is 7/1 ~ 8/8. If the base point
  # is observed on Jul 1, we want to choose observations which are observed on
  # Jul 1, Jul 2, Jul 3 because we don't have observations observed in June.
  # Same for the rest of the if statement.
  #! CHANGE THE PARAMETERS BELOW FOR YOUR PARTICULAR SETUP
  if (base_point$month == 7 & base_point$day == 1) {
    data_base <- oco2info_need[apply(oco2info_need['month'], 1, 
                               function(x) any(x == 7)), ]
    data_base <- data_base[apply(data_base['day'], 1, 
                                 function(x) any(x == 1 | x == 2 | x == 3)), ]
  } else if (base_point$month == 7 & base_point$day == 2) {
    data_base <- oco2info_need[apply(oco2info_need['month'], 1, 
                               function(x) any(x == 7)), ]
    data_base <- data_base[apply(data_base['day'], 1, 
                           function(x) any(x == 1 | x == 2 | x == 3 | x == 4)), ]
  } else if (base_point$month == 8 & base_point$day == 7) {
    data_base <- oco2info_need[apply(oco2info_need['month'], 1, 
                               function(x) any(x == 8)), ]
    data_base <- data_base[apply(data_base['day'], 1, 
                           function(x) any(x == 5 | x == 6 | x == 7 | x == 8)), ]
  } else if (base_point$month == 8 & base_point$day == 8) {
    data_base <- oco2info_need[apply(oco2info_need['month'], 1, 
                               function(x) any(x == 8)), ]
    data_base <- data_base[apply(data_base['day'], 1, 
                           function(x) any(x == 6 | x == 7 | x == 8)), ]
  } else {
    data_base <- oco2info_need[apply(oco2info_need['month'], 1, 
                               function(x) any(x == base_point$month)), ]
    data_base <- data_base[apply(data_base['day'], 1, 
                                 function(x) any(x == base_point$day | 
                                 x == base_point$day-2 | x == base_point$day-1 | 
                                 x == base_point$day+1 | x == base_point$day+2)), ]
  }
  # choose points within a circle with radius kDistance
  mylist <- list()
  for (k in 1:nrow(data_base)) {
    co2_cal <- data_base[k, ]
    distance <- distm(c(base_point$lon, base_point$lat), c(co2_cal$lon,co2_cal$lat))
    if (distance <= kDistance) {
      mylist[[k]] <- co2_cal
    }
  }
  co2_cal_in <- do.call('rbind', mylist)
  rawdata <- data.frame(co2_cal_in)

  #-----------------------------------#
  # Convert coordiantes               #
  #-----------------------------------#
  # Note: this system wouldn’t work if your domain 
  #       crossed the international dateline
  x_coord_list <- list()
  y_coord_list <- list()
  for (k in 1:nrow(rawdata)) {
    cal_point <- rawdata[k, ]
    y_coord <- distm(c(kRefLon, kRefLat), c(kRefLon, cal_point$lat))
    x_coord <- distm(c(kRefLon, kRefLat), c(cal_point$lon, kRefLat))
    x_coord_list[[k]] <- x_coord
    y_coord_list[[k]] <- y_coord
  }
  x_coord_list <- data.frame(do.call('rbind', x_coord_list))
  y_coord_list <- data.frame(do.call('rbind', y_coord_list))
  data_proc <- cbind(x_coord_list, y_coord_list, rawdata$co2, rawdata$land_ocean)
  colnames(data_proc) <- c('lon', 'lat', 'co2', 'land_ocean')
  # change the unit to km
  data_proc['lon'] <- 0.001 * data_proc['lon']
  data_proc['lat'] <- 0.001 * data_proc['lat']
  
  #-----------------------------------#
  # Calculate correlation length      #
  #-----------------------------------#
  if (nrow(data_proc) > 1){
    # plot the experimental variogram
    geodata <- as.geodata(data_proc, coords.col = c(1:2), data.col = 3)
    vario_exp <- variog(geodata, uvec = kBin, max.dist = kMaxDist, bin.cloud = T)
    # fit variogram
    vario_wls <- variofit(vario_exp, ini.cov.pars = kInitialValues, 
                          cov.model = kCovModel, fix.nugget = T, 
                          nugget = kNuggetLand, weights = kWeights)
    range <- 3 * vario_wls$cov.pars[2]
    sill <- vario_wls$cov.pars[1]

    # Re-fit
    if (range < kReFitThreshold) {
      vario_exp <- variog(geodata, uvec = kBinAgain, max.dist = kMaxDistAgain,
                          bin.cloud = T)
      vario_wls <- variofit(vario_exp, ini.cov.pars = kInitialValues, 
                            cov.model = kCovModel, fix.nugget = T, 
                            nugget = kNuggetLand, weights = kWeights)
      range <- 3 * vario_wls$cov.pars[2]
      sill <- vario_wls$cov.pars[1]
    }
  } else {
    sill <- 2
    range <- 0
  }
  # If range is bigger than 1000km, set to 1000km
  if (range > kMaxRange) {
    range <- kMaxRange
  }
  # If range is smaller than 20km, set to 20km
  if (kMinRange <= kMinRange){
      kMinRange <- kMinRange
  }
  
  #---------------------------------------#
  # Use ordinary kriging to estimate co2  #
  #---------------------------------------#
  # convert coordinates for base point
  base_point_y_coord <- 0.001 * distm(c(kRefLon, kRefLat), c(kRefLon, base_point$lat))
  base_point_x_coord <- 0.001 * distm(c(kRefLon, kRefLat), c(base_point$lon, kRefLat))
  base_point_proc <- data.frame(cbind(base_point_x_coord, base_point_y_coord))
  colnames(base_point_proc) <- c('lon', 'lat')
  # estimate co2
  if (nrow(data_proc) > 1){
    ok <- krige.conv(geodata, 
                     loc = matrix(c(base_point_proc$lon, base_point_proc$lat), 
                                  nrow = 1, ncol = 2),
                     krige = krige.control(type.krige = kKriging,
                                           cov.model = kKrigingCovModel,
                                           cov.pars=c(sill, 1/3 * range), 
                                           nugget = kNuggetLand, micro.scale = 0),
                     output = output.control(signal = TRUE))
    co2_est <- ok$predict
    est_error <- ok$krige.var
  } else {
    co2_est <- base_point$co2
    est_error <- 0
  }
  return(list(range, co2_est, sill, est_error))
} # end of CalculateRange function

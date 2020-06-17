#-----------------------------------------------------------------------------#
# SCRIPT: CalculateRange.R                                                    #
# PURPOSE: Calculate correlation length and use ordinary kriging              #
#          to estimate co2.                                                   #
# Xiaoling Liu, June 11, 2019, adapted by X.Liu and A.Weinbren in Jun 2020    #
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
  #! CHANGE THE PARAMETERS BELOW FOR YOUR PARTICULAR SETUP
  # kDistance is the cut-off distance for choosing observations around 
  # the base point to plot the empirical variogram.
  # (unit: m)
  kDistance <- 2000000
  # Latitude and longtitude of the reference point.
  # The reference point is an original point on one of the edges of
  # the domain. It is used for converting the observation coordinates
  # from lat/lon to kilometers.
  kRefLon <- -180
  kRefLat <- 9
  # Parameters for plotting empirical variogram.
  # kBin represents bin. kMaxDist represents cut-off distance.
  # Note: kBin recommends 3~6. kMaxDist recommends to be equal to or smaller than kDistance.
  # (unit for kMaxDist: km)
  kBin <- 6
  kMaxDist <- 2000
  # Initial parameters for fitting variogram.
  # expand.grid(sill, 1/3 * range)
  # sill represents variability in observations that are located at long distance
  # from one another.
  # range represents distance at which two measurements are practically uncorrelated.
  kInitialValues <- expand.grid(seq(0, 100, by = 5), seq(100, 135, by = 5))
  # Parameters for fitting variogram.
  # kCovModel represents model used to fit variogram.
  # kNugget represents variability that is uncorrelated from one observation
  # to the next.
  # kWeights represents type weights used in the loss function.
  kCovModel <- 'exp'
  kNugget <- (1.3)^2
  kWeights <- 'npairs'
  # kReFitThreshold is the threshold to determine whether to re-calculate range.
  # Reason for using 10: in our case project, most of the fitting variograms 
  # with range smaller than 10km had bad fit, thus we re-set to kBinAgain and kMaxDistAgain 
  # and re-fit.
  kReFitThreshold <- 10
  kBinAgain <- 5
  kMaxDistAgain <- 200
  # Maximum and minimum threshold for range.
  # We set range that are larger than 1000km to be 1000km 
  # and range that are smaller than 20km to be 20km.
  kMaxRange <- 1000
  kMinRange <- 20
  # Parameters for kriging
  # kKriging represents type of kriging to be performed.
  # kKrigingCovModel represents the name of the model for the correlation function. 
  kKriging <- 'ok'
  kKrigingCovModel <- 'exponential'

  #-----------------------------------#
  # Data pre-process                  #
  #-----------------------------------#
  # choose obs within the same day
  data_base <- subset(oco2info_need, month == base_point$month)
  data_base <- subset(data_base, day == base_point$day)
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
  # Convert coordinates               #
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
    # plot the empirical variogram
    geodata <- as.geodata(data_proc, coords.col = c(1:2), data.col = 3)
    vario_exp <- variog(geodata, uvec = kBin, max.dist = kMaxDist, bin.cloud = T)
    # fit variogram
    vario_wls <- variofit(vario_exp, ini.cov.pars = kInitialValues, 
                          cov.model = kCovModel, fix.nugget = T, 
                          nugget = kNugget, weights = kWeights)
    range <- 3 * vario_wls$cov.pars[2]
    sill <- vario_wls$cov.pars[1]

    # Re-fit
    if (range < kReFitThreshold) {
      vario_exp <- variog(geodata, uvec = kBinAgain, max.dist = kMaxDistAgain,
                          bin.cloud = T)
      vario_wls <- variofit(vario_exp, ini.cov.pars = kInitialValues, 
                            cov.model = kCovModel, fix.nugget = T, 
                            nugget = kNugget, weights = kWeights)
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
                                           cov.pars = c(sill, 1/3 * range), 
                                           nugget = kNugget, micro.scale = 0),
                     output = output.control(signal = TRUE))
    co2_est <- ok$predict
    est_error <- ok$krige.var
  } else {
    co2_est <- base_point$co2
    est_error <- 0
  }
  return(list(range, co2_est, sill, est_error))
} # end of CalculateRange function

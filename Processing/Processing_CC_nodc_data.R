#####
# Processing data from the NODC database ftp://ftp.nodc.noaa.gov/pub/data.nodc/nodc/archive/data/0129374/
#####

# needed packages
setwd("D:/Claudia_Marlon/Climate_change")

pcakages <- c("raster", "maps", "maptools", "rgeos", "rgdal", "ncdf4")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

#####
#Precipitation
years <- c(1950:1979, 2010:2013) # years
sums <- list() # list for results

for (i in 1:length(years)) {
  months <- list()
  for (j in 1:12) {
    if (j < 10) {
      jc <- paste(0, j, sep = "")
    }else {
      jc <- j
    }
    
    name <- paste("Data1/livneh_NAmerExt_15Oct2014.", # name of files in your computer
                  years[i], jc, ".mon.nc", sep = "")
    
    # before
    months[[j]] <- stack(name, varname = "Prec") # reading raster files
  }
  
  sum <- stack(months) # stackin files
  vals <- values(sum) # getting only values from raster
  
  sums[[i]] <- apply(vals, 1, sum) # sum precipitations per year
  
  cat("year", i, "of", length(years), "\n")
}

# preparing big tables of data for calculating means for each period
before <- do.call(cbind, sums[1:30]) 
after <- do.call(cbind, sums[31:34])

# mean calculation
meanb <- sum[[1]]
meanb[] <- apply(before, 1, mean)

meana <- sum[[1]]
meana[] <- apply(after, 1, mean)

# difference in values
difba <- meana - meanb

##saving files just in case
dir.create("Processed_rasters")

writeRaster(meanb, filename = "Processed_rasters/prec_before.tif", format = "GTiff")
writeRaster(meana, filename = "Processed_rasters/prec_after.tif", format = "GTiff")
writeRaster(difba, filename = "Processed_rasters/prec_change.tif", format = "GTiff")

#####
#Temperature
# same thing except for calculation of means to get mean temperature per year
years <- c(1950:1979, 2010:2013)
means <- list()

for (i in 1:length(years)) {
  months <- list()
  for (j in 1:12) {
    if (j < 10) {
      jc <- paste(0, j, sep = "")
    }else {
      jc <- j
    }
    
    name <- paste("Data1/livneh_NAmerExt_15Oct2014.", years[i], jc, ".mon.nc", sep = "")
    
    # before
    min <- stack(name, varname = "Tmin")
    max <- stack(name, varname = "Tmax")
    
    months[[j]] <- (max + min) / 2 # mean temperature = mean of Tmax and Tmin
  }
  
  mean <- stack(months)
  vals <- values(mean)
  
  means[[i]] <- apply(vals, 1, mean) # mean instead of sum 
  
  cat("year", i, "of", length(years), "\n")
}

before <- do.call(cbind, means[1:30])
after <- do.call(cbind, means[31:34])

meanb <- mean[[1]]
meanb[] <- apply(before, 1, mean)

meana <- mean[[1]]
meana[] <- apply(after, 1, mean)

difba <- meana - meanb

##saving files just in case
writeRaster(meanb, filename = "Processed_rasters/tmean_before.tif", format = "GTiff")
writeRaster(meana, filename = "Processed_rasters/tmean_after.tif", format = "GTiff")
writeRaster(difba, filename = "Processed_rasters/tmean_change.tif", format = "GTiff")

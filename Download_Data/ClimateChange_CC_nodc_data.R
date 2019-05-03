#####
# Download data from the NODC database ftp://ftp.nodc.noaa.gov/pub/data.nodc/nodc/archive/data/0129374/
#####

# working directory
setwd("D:/Claudia_Marlon/Climate_change")

dir.create("Data1") # where to save the data

years <- c(1950:1979, 2010:2013) # years to which the data will be downloaded

# download in loop
for (i in years) {
  for (j in 1:12) {
    if (j < 10) {
      jc <- paste(0, j, sep = "")
    }else {
      jc <- j
    }
    
    name <- paste("Data1/livneh_NAmerExt_15Oct2014.", # name of files in your computer
                  i, jc, ".mon.nc", sep = "") 
    
    url <- paste("https://data.nodc.noaa.gov/nodc/archive/data/0129374/monthly/", # url of files in database
                 "livneh_NAmerExt_15Oct2014.", i, jc, ".mon.nc", sep = "")
    
    down <- download.file(url = url, destfile = name, mode = "wb", quiet = TRUE) # donwload the file 
  }
  
  cat("year", i, "of", "1950-1979 and 2010-2013\n")
}

#####
# Download and unzip Landsta Forest Cover Change 
# data from Global Land Cover Facility
#####

# needed packages
pcakages <- c("raster", "maps", "maptools", "rgeos", "rgdal", "R.utils")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# directory
setwd("D:/Claudia_Marlon/WRS2_descending")

# bringing bird hotspots and landsat grid
occ <- read.csv("Z:/ENM_USACan/BeforeAfterComparisons/HotSpot_changes/HotSpot_changes_R.csv")
lsat <- readOGR(dsn = ".", layer = "WRS2_descending")

points <- occ[, c(1,7:8)] # only coordinates

WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") # projection

# coordinates to spatialpoints and reprojection to mantain areas
occ_sp <- SpatialPointsDataFrame(coords = points[, 2:3], data = points,
                                     proj4string = WGS84)

# projection for mantaining areas
centroid <- gCentroid(occ_sp, byid = FALSE)

AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", centroid@coords[1],
                      " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", sep = ""))

occ_pr <- spTransform(occ_sp, AEQD)

# buffer
buff_area <- gBuffer(occ_pr, width = 50000, byid = TRUE)

# reproject to wgs84
buff <- spTransform(buff_area, WGS84)

writeOGR(buff, ".", "hostpot_buffer", driver = "ESRI Shapefile")

lsatp <- spTransform(lsat, WGS84)

# intersection 
intersection <- lsatp[buff,]

# getting the paths and rows
places <- intersection@data[, 7:8]
write.csv(places, "patrow.csv", row.names = F)
places <- read.csv("patrow.csv")

# download in loop
setwd("D:/Claudia_Marlon") # general directory
dir.create("LCC_files") # new folder for saving downloads 90-00
dir.create("LCC_files1") # new folder for saving downloads 00-05

res <- list() #for seeing which ones wwere not found

for (i in 1:dim(places)[1]) {
  # for adjusting path and row names
  if(places[i, 1] < 10){
    path <- paste("00", places[i, 1], sep = "")
  }
  if(places[i, 1] >= 10 & places[i, 1] < 100){
    path <- paste("0", places[i, 1], sep = "")
  }
  if(places[i, 1] >= 100) {
    path <- places[i, 1]
  }
  if(places[i, 2] < 10){
    row <- paste("00", places[i, 2], sep = "")
  }
  if(places[i, 2] >= 10 & places[i, 2] < 100){
    row <- paste("0", places[i, 2], sep = "")
  }
  if(places[i, 2] >= 100) {
    row <- places[i, 2]
  }
  
  # name of files
  name <- paste("LCC_files", paste("path", path, "row", row, "90-00.gz", sep = "_"), sep = "/")
  nam <- paste("LCC_files", paste("path", path, "row", row, "90-00", sep = "_"), sep = "/")
  names <- paste("LCC_files", paste("path", path, "row", row, "90-00.tif", sep = "_"), sep = "/")
  name1 <- paste("LCC_files1", paste("path", path, "row", row, "00-05.gz", sep = "_"), sep = "/")
  nam1 <- paste("LCC_files1", paste("path", path, "row", row, "00-05", sep = "_"), sep = "/")
  names1 <- paste("LCC_files1", paste("path", path, "row", row, "00-05.tif", sep = "_"), sep = "/")
  
  # url of each image 90-00
  url <- paste("ftp://ftp.glcf.umd.edu/glcf/LandsatFCC/stow/GLCF.NPM.AA2-004.00.FCC1990_2000_v1/WRS2/",
               paste("p", path, "/r", row, sep = ""),
               paste("/p", path, "r", row, "_FCC_19902000/", sep = ""),
               paste("p", path, "r", row, "_FCC_19902000_CM.tif.gz", sep = ""), sep = "")
  
  down <- try(download.file(url = url, destfile = name, mode = "wb", quiet = FALSE), # donwload the zipped file
              silent = TRUE)
  uz <- try(gunzip(filename = name, remove = TRUE), silent = TRUE) # unzip
  rn <- try(file.rename(nam, names), silent = TRUE) # rename files because they don't have .tif extention
  result <- class(down) # for seeing which ones gave errors
  
  # url of each image 00-05
  url1 <- paste("ftp://ftp.glcf.umd.edu/glcf/LandsatFCC/stow/GLCF.NPM.AA2-003.00.FCC2000_2005_v1/WRS2/",
                paste("p", path, "/r", row, sep = ""),
                paste("/p", path, "r", row, "_FCC_20002005/", sep = ""),
                paste("p", path, "r", row, "_FCC_20002005_CM.tif.gz", sep = ""), sep = "")
  
  down1 <- try(download.file(url = url1, destfile = name1, mode = "wb", quiet = FALSE), # donwload the zipped file
               silent = TRUE)
  uz1 <- try(gunzip(filename = name1, remove = TRUE), silent = TRUE) # unzip
  rn1 <- try(file.rename(nam1, names1), silent = TRUE) # rename files because they don't have .tif extention
  result1 <- class(down1) # for seeing which ones gave errors
  
  res[[i]] <- c(result, result1) # adding to the list of tryings
  names(res[[i]]) <- c(name, name1) #naming each try as the image to be downloaded
}

results <- unlist(res) # unlist tries 
nodata <- results[results == "try-error"] # getting only the ones that were not found
nodattable <- cbind(files = names(nodata), nodata) # making it a table

write.csv(nodattable, "files_nodata.csv", row.names = FALSE) # saving the table

# seeing the ones that were not found in a plot
## to check which polygons are the ones with no data
nodat <- paste(c(13, 13, 45, 48, 42, 17),
               c(33, 34, 34, 28, 37, 42), sep = "_") 

intersection@data$check <- paste(intersection@data$PATH,
                                 intersection@data$ROW, sep = "_")

last_nodata <- intersection[intersection@data$check %in% nodat, ]

## world map
w_map <- map(database = "world", fill = TRUE, plot = FALSE) # map of the world

w_po <- sapply(strsplit(w_map$names, ":"), function(x) x[1]) # preparing data to create polygon
world <- map2SpatialPolygons(w_map, IDs = w_po, proj4string = WGS84) # map to polygon

## plotting
plot(intersection)
plot(last_nodata, col = "red", add = TRUE)
plot(world, border = "blue", add = TRUE)
plot(buff, border = "purple", add = TRUE)

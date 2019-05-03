#####
# Seeing the difference of information in a combined raster
#####

# needed packages
pcakages <- c("raster", "rgdal", "sp")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# directory
setwd("D:/Claudia_Marlon")

# data
combin <- raster("LCC_combine/lcc_comb_mask1.tif") #combine
buff <- readOGR(dsn = "WRS2_descending", layer = "hostpot_buffer") # hotspot buffers

IDs <- 1:16 # values of interest in combine

# assessing type of change and mean change
statistics <- list()

for (i in buff$ID) {
  # getting data per hotspot buffer
  buffp <- as(buff[buff$ID == i, ], "SpatialPolygons")
  cp <- crop(combin, buffp)
  mask <- rasterize(buffp, cp, mask = TRUE)
  
  # assessing type of change
  comvals <- vector()
  for (j in 1:length(IDs)) {
    comvals[j] <- sum(na.omit(values(mask)) == IDs[j])
  }
  
  names(comvals) <- IDs
  
  statistics[[i]] <- comvals
  cat("Hotspot", i, "completed\n")
}

# preparing results
names(statistics) <- buff$ID

statisticsres <- do.call(rbind, statistics)

results <- data.frame(buff$ID, statisticsres)
names(results) <- c("Hotspot_ID", paste("Value", 1:16, sep = "_"))

# writing results
dir.create("LCC_results")
write.csv(results, "LCC_results/lcc_combine_pix_count.csv", row.names = FALSE)

# calculating results summary and ratios of change
merged <- data.frame((statisticsres[, 1] + statisticsres[, 4] + statisticsres[, 9] + statisticsres[, 10]), 
                     (statisticsres[, 3] + statisticsres[, 8] + statisticsres[, 12] + statisticsres[, 16]),
                     (statisticsres[, 2] + statisticsres[, 5] + statisticsres[, 14] + statisticsres[, 15]),
                     (statisticsres[, 6] + statisticsres[, 7] + statisticsres[, 11] + statisticsres[, 13]))

gain_ratio <- merged[, 3] / apply(merged[, 3:4], 1, sum)
loss_ratio <- merged[, 2] / apply(merged[, 1:2], 1, sum)

results1 <- data.frame(buff$ID, merged, gain_ratio, loss_ratio)
names(results1) <- c("Hotspot_ID", "Persistent_forest", "Forest_loss", "Forest_gain",
                     "Persistent_non_forest", "Forest_gain_ratio_(G/GpPNF)", "Forest_loss_ratio_(L/LpPF)")

write.csv(results1, "LCC_results/lcc_gain_loss_ratios.csv", row.names = FALSE)

#####
# Counting type of climate change USA-CAN
#####

# needed packages
pcakages <- c("raster", "rgdal", "sp")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# directory
setwd("D:/Claudia_Marlon/Climate_change/")

# precipiation, temperature, and buffers of hotspots
prec <- raster("Processed_rasters/prec_change.tif") 
temp <- raster("Processed_rasters/tmean_change.tif") 
vars <- list(prec, temp)

buff <- readOGR(dsn = "D:/Claudia_Marlon/WRS2_descending", layer = "hostpot_buffer") # hotspot buffers

# assessing type of change and mean change
statistics <- list()

for (h in 1:length(vars)) {
  stats <- list()
  for (i in buff$ID) {
    # getting data per hotspot buffer
    buffp <- as(buff[buff$ID == i, ], "SpatialPolygons")
    cp <- crop(vars[[h]], buffp)
    mask <- rasterize(buffp, cp, mask = TRUE)
    
    # type of change and mean change calculation
    inc <- sum(na.omit(values(mask)) > 0) / sum(!is.na(values(mask)))
    dec <- sum(na.omit(values(mask)) < 0) / sum(!is.na(values(mask)))
    sta <- sum(na.omit(values(mask)) == 0) / sum(!is.na(values(mask)))
    
    change <- mean(na.omit(values(mask)))
    
    all <- c(inc, dec, sta, change)
    names(all) <- c("Prop_increase", "Prop_pix_decrease", "Prop_pix_stable",
                    "Mean_change")
    
    stats[[i]] <- all
    
    cat("   Hotspot", i, "completed\n")
  }
  
  statistics[[h]] <- stats
  cat("variable", h, "of", length(vars), "completed\n")
}

# preparing results
prec_stats <- cbind(buff$ID, do.call(rbind, statistics[[1]]))
prec_stats[, 5] <- prec_stats[, 5] * 10 # precipitation was in cm, changing to mm
colnames(prec_stats)[1] <- "Hotspot_ID"

temp_stats <- cbind(buff$ID, do.call(rbind, statistics[[2]]))
colnames(temp_stats)[1] <- "Hotspot_ID"


# writing results
dir.create("CC_results")
write.csv(prec_stats, "CC_results/precipitation_change.csv", row.names = FALSE)
write.csv(temp_stats, "CC_results/temperature_change.csv", row.names = FALSE)

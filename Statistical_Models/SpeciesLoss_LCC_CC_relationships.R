#####
# Relationship
#####

# needed packages
pcakages <- c("stats", "raster", "rgdal", "sp", "rgeos")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# directory
setwd("Z:/ENM_USACan/CC_LCC_USA_CAN/Relationships")

# data
incidences <- read.csv("HotSpot_changes_R.csv")

lcc <- read.csv("lcc_gain_loss_ratios.csv")
colnames(lcc)[1] <- "ID"

prec <- read.csv("precipitation_change.csv")
colnames(prec)[1] <- "ID"

temp <- read.csv("temperature_change.csv")
colnames(temp)[1] <- "ID"

# merging data
all_data <- merge(incidences, lcc, by = "ID")
all_data <- merge(all_data, prec, by = "ID")
all_data <- merge(all_data, temp, by = "ID")

#####
# searching for potential relationships
## losses
x11()
par(mfrow = c(4, 2), mar = c(4.5, 4.5, 0.5, 0.5), cex = 0.7)
plot(all_data$Forest_gain_ratio_.G.GpPNF., all_data$loss_, xlab = "Proportion of forest cover gain",
     ylab = "Proportion of incidence loss")
plot(all_data$Forest_loss_ratio_.L.LpPF., all_data$loss_, xlab = "Proportion of forest cover loss",
     ylab = "")
plot(all_data$Mean_change.x, all_data$loss_, xlab = "Mean precipitation change",
     ylab = "Proportion of incidence loss")
plot(all_data$Mean_change.y, all_data$loss_, xlab = "Mean temperature change",
     ylab = "")
plot(all_data$Prop_pix_increase.x, all_data$loss_, xlab = "Proportion of area with increase in precipitation",
     ylab = "Proportion of incidence loss")
plot(all_data$Prop_pix_decrease.x, all_data$loss_, xlab = "Proportion of area with decrease in precipitation",
     ylab = "")
plot(all_data$Prop_pix_increase.y, all_data$loss_, xlab = "Proportion of area with increase in temperature",
     ylab = "Proportion of incidence loss")
plot(all_data$Prop_pix_decrease.y, all_data$loss_, xlab = "Proportion of area with decrease in temperature",
     ylab = "")

dev.off()

# testing for correlations
all_for_corr <- data.frame(all_data[, c("jac", "gain_", "stable_", "loss_", "Forest_gain_ratio_.G.GpPNF.",
                                        "Forest_loss_ratio_.L.LpPF.", "Mean_change.x", "Mean_change.y",
                                        "Prop_pix_increase.x", "Prop_pix_decrease.x", "Prop_pix_increase.y",
                                        "Prop_pix_decrease.y")])

correlations <- cor(all_for_corr)

# only for losses
loss_for_corr <- data.frame(all_data[, c("loss_", "Forest_gain_ratio_.G.GpPNF.",
                                        "Forest_loss_ratio_.L.LpPF.", "Mean_change.x", "Mean_change.y",
                                        "Prop_pix_increase.x", "Prop_pix_decrease.x", "Prop_pix_increase.y",
                                        "Prop_pix_decrease.y")])

write.csv(loss_for_corr, "loss_other_vars.csv", row.names = TRUE)

correlations_loss <- cor(loss_for_corr)
x11()
pairs(correlations_loss)

# correlation for loss and other variables and significance
cor_lossc <- c(loss_meanP = cor.test(all_data$Mean_change.x, all_data$loss_)$estimate, 
               loss_meanT = cor.test(all_data$Mean_change.y, all_data$loss_)$estimate, 
               loss_IP = cor.test(all_data$Prop_pix_increase.x, all_data$loss_)$estimate, 
               loss_DP = cor.test(all_data$Prop_pix_decrease.x, all_data$loss_)$estimate,
               loss_IT = cor.test(all_data$Prop_pix_increase.y, all_data$loss_)$estimate, 
               loss_DT = cor.test(all_data$Prop_pix_decrease.y, all_data$loss_)$estimate, 
               loss_FLR = cor.test(all_data$Forest_loss_ratio_.L.LpPF., all_data$loss_)$estimate, 
               loss_FGR = cor.test(all_data$Forest_gain_ratio_.G.GpPNF., all_data$loss_)$estimate)

cor_lossp <- c(loss_meanP = cor.test(all_data$Mean_change.x, all_data$loss_)$p.value, 
               loss_meanT = cor.test(all_data$Mean_change.y, all_data$loss_)$p.value, 
               loss_IP = cor.test(all_data$Prop_pix_increase.x, all_data$loss_)$p.value, 
               loss_DP = cor.test(all_data$Prop_pix_decrease.x, all_data$loss_)$p.value,
               loss_IT = cor.test(all_data$Prop_pix_increase.y, all_data$loss_)$p.value, 
               loss_DT = cor.test(all_data$Prop_pix_decrease.y, all_data$loss_)$p.value, 
               loss_FLR = cor.test(all_data$Forest_loss_ratio_.L.LpPF., all_data$loss_)$p.value, 
               loss_FGR = cor.test(all_data$Forest_gain_ratio_.G.GpPNF., all_data$loss_)$p.value)

correlaion_ress <- rbind(cor_lossc, cor_lossp)
rownames(correlaion_ress) <- c("Correlation", "p_value")

# write correlation results
write.csv(correlaion_ress, "correlation_loss.csv", row.names = TRUE)

#####
# checking for spatial autocorrelation
## matrix of distnaces
cords <- all_data[, c("Longitude", "Latitude")]

mat_dis <- dist(cords) # object type distance to be used after

md <- sapply(pointDistance(cords, lonlat = T), as.numeric) #distances in meters
md <- md[!is.na(md)] # erasing NAs
md <- md[md != 0] # erasing zeros

## truncation of the matrix distance of truncation (1000 km ~ ) truncation = D x 4
tdist <- 1000000 # truncation distance
md <- ifelse(md > tdist, md * 4, md) # tuncating values bigger than 1000 km

mat_dis[] <- md # replacing values in object type dist

## calculating eigenvectors via principal coordinate analysis
pcoa <- cmdscale(mat_dis, k = 107, eig = TRUE)

eigs <- pcoa$eig # all eigenvalues
pos_eigs <- 1:length(eigs[eigs > 0]) # position of positive eigenvalues

eigens <- pcoa$points[, pos_eigs] # getting only eigenvectors with positive eigenvalues

## selecting eigenvectors that should be considered in glms ********THREE EIGENVECTORS WERE SELECTED**********
### plots to select eigenvectors
jpeg(filename = "Eigenvector_selection.jpg", width = 80, height = 76,
     units = "mm", res = 600)
par(mar = c(4.5, 4.6, 0.5, 0.5), cex = 0.7)
plot(1:108, eigs / max(eigs), type = "l", col = "black", las = 1, xlab = "", ylab = "") # all eigenvalues
abline(h = 0, col = "grey75", lwd = 1, lty = 2)
title(xlab = "Rank of eigenvalues", ylab = "", cex.lab = 1.2, line = 3)
title(xlab = "", ylab = "Eigenvalues (normalized)", cex.lab = 1.2, line = 3.3)
abline(v = 3, col = "red", lwd = 1, lty = 2)
legend("topright", legend = "Third eigenvalue", lwd = 1, lty = 2, col = "red", bty = "n")
dev.off()

### fitting models for selecting vectors
#Incidence_loss <- all_data$loss_
#eigens_models <- list()

#for (i in 1:dim(eigens)[2]) {
#  loss_eigens <- data.frame(Incidence_loss, eigens[, 1:i])
#  eigens_models[[i]] <- glm(Incidence_loss ~ ., data = loss_eigens)
#}

#####
# regresions for predictions with all data
## making predictors independent and excluding correlated more than 0.8
Incidence_loss <- all_data$loss_
Precip_increase <- all_data$Prop_pix_increase.x
Temp_increase <- all_data$Prop_pix_increase.y
Precip_change <- all_data$Mean_change.x
Temp_change <- all_data$Mean_change.y
Forest_gain <- all_data$Forest_gain_ratio_.G.GpPNF.
Forest_loss <- all_data$Forest_loss_ratio_.L.LpPF.
E1 <- eigens[, 1]
E2 <- eigens[, 2]
E3 <- eigens[, 3]

## regressions with all combinations of predictors
variables <- c("Precip_increase", "Temp_increase", "Precip_change",
               "Temp_change", "Forest_gain", "Forest_loss")
var_comba <- list()
var_combi <- list()
objectsa <- list()
objectsi <- list()

for (i in 1:length(variables)) {
  comb <- combn(x = variables, m = i)
  comb_va <- vector()
  comb_vi <- vector()
  objecta <- vector()
  objecti <- vector()
  
  for (j in 1:dim(comb)[2]) {
    comb_va[j] <- paste("L", i, "V", j, "a", paste(paste(" <- glm(Incidence_loss", paste(comb[, j], collapse = " + "), sep = " ~ "),
                                                   " + E1 + E2 + E3, family = gaussian)\n", sep = ""), sep = "")
    comb_vi[j] <- paste("L", i, "V", j, "i", paste(paste(" <- glm(Incidence_loss", paste(comb[, j], collapse = " * "), sep = " ~ "), 
                                                   " + E1 + E2 + E3, family = gaussian)\n", sep = ""), sep = "")
    objecta[j] <- paste("L", i, "V", j, "a", sep = "")
    objecti[j] <- paste("L", i, "V", j, "i", sep = "")
    
  }
  
  var_comba[[i]] <- comb_va
  var_combi[[i]] <- comb_vi
  objectsa[[i]] <- paste(paste(objecta, collapse = ", "), "\n", sep = "")
  objectsi[[i]] <- paste(paste(objecti, collapse = ", "), "\n", sep = "")
}

var_combar <- do.call(c, var_comba)
var_combir <- do.call(c, var_combi)

objectsr <- paste("aic_models <- AIC(", paste(do.call(c, objectsa), collapse = ", "), ", \n",
                  paste(do.call(c, objectsi), collapse = ", "), ")", sep = "")

sink("model_cal1.R")
cat(var_combar, "\n")
cat(var_combir, "\n")
cat(objectsr, "\n")
sink()

source("model_cal1.R")
aic_models$delta_AIC <- aic_models$AIC - min(aic_models$AIC)
aic_models <- aic_models[order(aic_models$delta_AIC),]

# writing results of model AIC comparison
write.csv(aic_models, "all_comb_glms_aic1.csv", row.names = TRUE)

# searching for best models
best_aic <- aic_models[aic_models$delta_AIC <= 2, ]
best_aic <- best_aic[order(best_aic$delta_AIC), ]
best_aic$formula <- as.character(c(L2V4i$formula))

coefs <- list(L2V4i$coefficients)

coef_names <- rep(NA, length(unique(names(do.call(c, coefs)))))
names(coef_names) <- sort(unique(names(do.call(c, coefs))))

coef_lis <- list()
for (i in 1:length(coefs)) {
  for (j in 1:length(coef_names)) {
    coef_names[j] <- ifelse(names(coef_names)[j] %in% names(coefs[[i]]), 
                             coefs[[i]][names(coefs[[i]]) == names(coef_names)[j]], NA)
  }
  
  coef_lis[[i]] <- na.omit(coef_names)
}

## completing the table with coeficients 
best_aicc <- data.frame(best_aic, do.call(rbind, coef_lis))

## writing best models info
write.csv(best_aicc, "best_glms_aic1.csv", row.names = TRUE)

## exploring best models
L2V4i

x11()
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 1.5, 0.5), cex = 0.7)
plot(L2V4i)

## significance of effects
fmodel <- anova(L2V4i, test = "F")
sig_effecs <- as.data.frame(fmodel)

## writing significnce of effects
sink("significnace_best_glms_aic1.txt")

cat(row.names(best_aic), "\t", best_aic[, 4], "\n")
print(sig_effecs)

sink()


######################################
#####
# Exploring outlier exclussion

# which one is that outlier
prec <- raster("D:/Claudia_Marlon/Climate_change/Processed_rasters/prec_change.tif")
tempe <- raster("D:/Claudia_Marlon/Climate_change/Processed_rasters/tmean_change.tif") 
buff <- readOGR(dsn = "D:/Claudia_Marlon/WRS2_descending", layer = "hostpot_buffer") # hotspot buffers

x11()
par(mar = c(1.7, 1.7, 0.5, 3.5), cex = 0.9)
plot(prec)
#plot(tempe)
#map("world", xlim = c(-130, -63), ylim = c(15, 53))
plot(buff, add = TRUE, lwd = 3)
plot(buff[101, ], add = TRUE, border = "red" , lwd = 3)
legend("bottomright", c("All hotspots", "Outlier"), pch = 1, bty = "n", 
       col = c("black", "red"), inset = 0.09, lty = 0, lwd = 3, cex = 1.5)

# all analyses without outlier 
# testing for correlations
# only for losses
loss_for_corrw <- data.frame(all_data[-101, c("loss_", "Forest_gain_ratio_.G.GpPNF.",
                                             "Forest_loss_ratio_.L.LpPF.", "Mean_change.x", "Mean_change.y",
                                             "Prop_pix_increase.x", "Prop_pix_decrease.x", "Prop_pix_increase.y",
                                             "Prop_pix_decrease.y")])

# writing data with no outlier
write.csv(loss_for_corrw, "loss_other_vars_no_outlier.csv", row.names = TRUE)
 
correlations_lossw <- cor(loss_for_corrw)
x11()
pairs(correlations_lossw)

# correlation for loss and other variables and significance
cor_losscw <- c(loss_meanP = cor.test(loss_for_corrw$Mean_change.x, loss_for_corrw$loss_)$estimate, 
               loss_meanT = cor.test(loss_for_corrw$Mean_change.y, loss_for_corrw$loss_)$estimate, 
               loss_IP = cor.test(loss_for_corrw$Prop_pix_increase.x, loss_for_corrw$loss_)$estimate, 
               loss_DP = cor.test(loss_for_corrw$Prop_pix_decrease.x, loss_for_corrw$loss_)$estimate,
               loss_IT = cor.test(loss_for_corrw$Prop_pix_increase.y, loss_for_corrw$loss_)$estimate, 
               loss_DT = cor.test(loss_for_corrw$Prop_pix_decrease.y, loss_for_corrw$loss_)$estimate, 
               loss_FLR = cor.test(loss_for_corrw$Forest_loss_ratio_.L.LpPF., loss_for_corrw$loss_)$estimate, 
               loss_FGR = cor.test(loss_for_corrw$Forest_gain_ratio_.G.GpPNF., loss_for_corrw$loss_)$estimate)

cor_losspw <- c(loss_meanP = cor.test(loss_for_corrw$Mean_change.x, loss_for_corrw$loss_)$p.value, 
               loss_meanT = cor.test(loss_for_corrw$Mean_change.y, loss_for_corrw$loss_)$p.value, 
               loss_IP = cor.test(loss_for_corrw$Prop_pix_increase.x, loss_for_corrw$loss_)$p.value, 
               loss_DP = cor.test(loss_for_corrw$Prop_pix_decrease.x, loss_for_corrw$loss_)$p.value,
               loss_IT = cor.test(loss_for_corrw$Prop_pix_increase.y, loss_for_corrw$loss_)$p.value, 
               loss_DT = cor.test(loss_for_corrw$Prop_pix_decrease.y, loss_for_corrw$loss_)$p.value, 
               loss_FLR = cor.test(loss_for_corrw$Forest_loss_ratio_.L.LpPF., loss_for_corrw$loss_)$p.value, 
               loss_FGR = cor.test(loss_for_corrw$Forest_gain_ratio_.G.GpPNF., loss_for_corrw$loss_)$p.value)

correlaion_ressw <- rbind(cor_losscw, cor_losspw)
rownames(correlaion_ressw) <- c("Correlation", "p_value")

# write correlation results
write.csv(correlaion_ressw, "correlation_loss_no_outlier.csv", row.names = TRUE)

# regresions
## one predictor
Incidence_loss <- loss_for_corrw$loss_
Precip_increase <- loss_for_corrw$Prop_pix_increase.x
Temp_increase <- loss_for_corrw$Prop_pix_increase.y
Precip_change <- loss_for_corrw$Mean_change.x
Temp_change <- loss_for_corrw$Mean_change.y
Forest_gain <- loss_for_corrw$Forest_gain_ratio_.G.GpPNF.
Forest_loss <- loss_for_corrw$Forest_loss_ratio_.L.LpPF.
E1 <- eigens[-101, 1]
E2 <- eigens[-101, 2]
E3 <- eigens[-101, 3]

# *************

# all combinations
variables <- c("Precip_increase", "Temp_increase", "Precip_change",
               "Temp_change", "Forest_gain", "Forest_loss")
var_comba <- list()
var_combi <- list()
objectsa <- list()
objectsi <- list()


for (i in 1:length(variables)) {
  comb <- combn(x = variables, m = i)
  comb_va <- vector()
  comb_vi <- vector()
  objecta <- vector()
  objecti <- vector()
  
  for (j in 1:dim(comb)[2]) {
    comb_va[j] <- paste("L", i, "V", j, "aw", paste(paste(" <- glm(Incidence_loss", paste(comb[, j], collapse = " + "), sep = " ~ "),
                                                   " + E1 + E2 + E3, family = gaussian)\n", sep = ""), sep = "")
    comb_vi[j] <- paste("L", i, "V", j, "iw", paste(paste(" <- glm(Incidence_loss", paste(comb[, j], collapse = " * "), sep = " ~ "), 
                                                   " + E1 + E2 + E3, family = gaussian)\n", sep = ""), sep = "")
    objecta[j] <- paste("L", i, "V", j, "aw", sep = "")
    objecti[j] <- paste("L", i, "V", j, "iw", sep = "")
    
  }
  
  var_comba[[i]] <- comb_va
  var_combi[[i]] <- comb_vi
  objectsa[[i]] <- paste(paste(objecta, collapse = ", "), "\n", sep = "")
  objectsi[[i]] <- paste(paste(objecti, collapse = ", "), "\n", sep = "")
}

var_combar <- do.call(c, var_comba)
var_combir <- do.call(c, var_combi)

objectsr <- paste("aic_modelsw <- AIC(", paste(do.call(c, objectsa), collapse = ", "), ", \n",
                  paste(do.call(c, objectsi), collapse = ", "), ")", sep = "")

sink("model_cal_no_outlier1.R")
cat(var_combar, "\n")
cat(var_combir, "\n")
cat(objectsr, "\n")
sink()

source("model_cal_no_outlier1.R")
aic_modelsw$delta_AIC <- aic_modelsw$AIC - min(aic_modelsw$AIC)

# writing results of model AIC comparison
write.csv(aic_modelsw, "all_comb_glms_aic_no_outlier1.csv", row.names = TRUE)

# searching for best models
best_aicw <- aic_modelsw[aic_modelsw$delta_AIC <= 2, ]
best_aicw <- best_aicw[order(best_aicw$delta_AIC), ]
best_aicw$formula <- as.character(c(L1V1aw$formula, L1V1iw$formula, L2V5aw$formula, L2V5iw$formula, 
                                    L3V4aw$formula, L2V3aw$formula,  L2V2aw$formula, L3V7aw$formula,
                                    L3V10aw$formula, L2V4aw$formula, L2V1aw$formula))

coefsw <- list(L1V1aw$coefficients, L1V1iw$coefficients, L2V5aw$coefficients, L2V5iw$coefficients, 
               L3V4aw$coefficients, L2V3aw$coefficients,  L2V2aw$coefficients, L3V7aw$coefficients,
               L3V10aw$coefficients, L2V4aw$coefficients, L2V1aw$coefficients)


coef_wnames <- rep(NA, length(unique(names(do.call(c, coefsw)))))
names(coef_wnames) <- sort(unique(names(do.call(c, coefsw))))

coef_list <- list()
for (i in 1:length(coefsw)) {
  for (j in 1:length(coef_wnames)) {
    coef_wnames[j] <- ifelse(names(coef_wnames)[j] %in% names(coefsw[[i]]), 
                             coefsw[[i]][names(coefsw[[i]]) == names(coef_wnames)[j]], NA)
  }
  
  coef_list[[i]] <- coef_wnames
}


# completing the table with coeficients 
best_aicwc <- data.frame(best_aicw, do.call(rbind, coef_list))[-2, ] # exclude model 2, it is the same than 1

# writing best models info
write.csv(best_aicwc, "best_glms_aic_no_outlier1.csv", row.names = TRUE)

# exploring best models
L3V10aw

x11()
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 1.5, 0.5), cex = 0.7)
plot(L3V10aw)

## significance of effects
f_models <- list(f_m1 <- as.data.frame(anova(L1V1aw, test = "F")),
                 f_m2 <- as.data.frame(anova(L2V5aw, test = "F")),
                 f_m3 <- as.data.frame(anova(L2V5iw, test = "F")),
                 f_m4 <- as.data.frame(anova(L3V4aw, test = "F")),
                 f_m5 <- as.data.frame(anova(L2V3aw, test = "F")),
                 f_m6 <- as.data.frame(anova(L2V2aw, test = "F")),
                 f_m7 <- as.data.frame(anova(L3V7aw, test = "F")),
                 f_m8 <- as.data.frame(anova(L3V10aw, test = "F")),
                 f_m9 <- as.data.frame(anova(L2V4aw, test = "F")),
                 f_m10 <- as.data.frame(anova(L2V1aw, test = "F")))

## writing f tables for all
sink("significnace_best_glms_aic_no_outlier1.txt")

for (i in 1:length(f_models)) {
  cat(row.names(best_aicw)[-2][i], "\t", best_aicw[-2, 4][i], "\n")
  print(f_models[[i]])
  cat("\n\n")
}

sink()
################################################################################
# 
# Purpose: produce Supplementary Table D.1
#        perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the bioassay data  
# Input: load 'result/bioassay_sep.RData' and 'result/bioassay_nov.RData' file
# Output: print Supplementary Table D.1 (September part)
#         print Supplementary Table D.1 (November part)
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/IMC.r', sep=''))

# September ############################################
load(paste(MS.PATH, 'result/bioassay_sep.RData',sep=''))
fitN <- IM.MVN(estN$para.est, estN$EST, Yc)
fitT <- IM.MVT(estT$para.est, estT$EST, Yc)
fitS <- IM.MSL(estS$para.est, estS$EST, Yc)

fitNC <- IM.MVNC(estNC$para.est, estNC$EST, Yc, cen)
fitTC <- IM.MVTC(estTC$para.est, estTC$EST, Yc, cen)
fitSC <- IM.MSLC(estSC$para.est, estSC$EST, Yc, cen, ITER=20, per.iter=2)

# EST
sep_EST <- cbind(c(fitN$EST,0), fitT$EST, fitS$EST,
                 c(fitNC$EST,0), fitTC$EST, fitSC$EST)
sep_EST_data <- as.data.frame(sep_EST)
colnames(sep_EST_data) <- c("MVN.EST", "MVT.EST", "MSL.EST", "MVNC.EST", "MVTC.EST", "MSLC.EST")
row.names(sep_EST_data) <- c("mu1", "mu2", "mu3", 
                             "sigma11", "sigma21", "sigma22",
                             "sigma31", "sigma32", "sigma33", "nu")
# SD
sep_SD <- cbind(c(fitN$SD,0), fitT$SD, fitS$SD,
                c(fitNC$SD,0), fitTC$SD, fitSC$SD)
sep_SD_data <- as.data.frame(sep_SD)
colnames(sep_SD_data) <- c("MVN.SE", "MVT.SE", "MSL.SE", "MVNC.SE", "MVTC.SE", "MSLC.SE")
row.names(sep_SD_data) <- c("sd.mu1", "sd.mu2", "sd.mu3", 
                            "sd.sigma11", "sd.sigma21", "sd.sigma22",
                            "sd.sigma31", "sd.sigma32", "sd.sigma33", "sd.nu")
# combined and sort
sep_IM <- cbind(sep_EST_data, sep_SD_data)
n_vars <- ncol(sep_EST_data)
idx <- c(rbind(1:n_vars, (n_vars + 1):(2 * n_vars)))
re_sep <- sep_IM[, idx]
re_sep[] <- lapply(re_sep, function(x) if(is.numeric(x)) round(x, 3) else x)

print(re_sep)

# November ############################################
rm(list = ls())
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/IMC.r', sep=''))

load(paste(MS.PATH, 'result/bioassay_nov.RData',sep=''))
fitN <- IM.MVN(estN$para.est, estN$EST, Yc)
fitT <- IM.MVT(estT$para.est, estT$EST, Yc)
fitS <- IM.MSL(estS$para.est, estS$EST, Yc)

fitNC <- IM.MVNC(estNC$para.est, estNC$EST, Yc, cen)
fitTC <- IM.MVTC(estTC$para.est, estTC$EST, Yc, cen)
fitSC <- IM.MSLC(estSC$para.est, estSC$EST, Yc, cen, ITER=20, per.iter=2)

# EST
nov_EST <- cbind(c(fitN$EST,0), fitT$EST, fitS$EST,
                 c(fitNC$EST,0), fitTC$EST, fitSC$EST)
nov_EST_data <- as.data.frame(nov_EST)
colnames(nov_EST_data) <- c("MVN.EST", "MVT.EST", "MSL.EST", "MVNC.EST", "MVTC.EST", "MSLC.EST")
row.names(nov_EST_data) <- c("mu1", "mu2", "mu3", 
                             "sigma11", "sigma21", "sigma22",
                             "sigma31", "sigma32", "sigma33", "nu")
# SD
nov_SD <- cbind(c(fitN$SD,0), fitT$SD, fitS$SD,
                c(fitNC$SD,0), fitTC$SD, fitSC$SD)
nov_SD_data <- as.data.frame(nov_SD)
colnames(nov_SD_data) <- c("MVN.SE", "MVT.SE", "MSL.SE", "MVNC.SE", "MVTC.SE", "MSLC.SE")
row.names(nov_SD_data) <- c("sd.mu1", "sd.mu2", "sd.mu3", 
                            "sd.sigma11", "sd.sigma21", "sd.sigma22",
                            "sd.sigma31", "sd.sigma32", "sd.sigma33", "sd.nu")
# combined and sort
nov_IM <- cbind(nov_EST_data, nov_SD_data)
n_vars <- ncol(nov_EST_data)
idx <- c(rbind(1:n_vars, (n_vars + 1):(2 * n_vars)))
re_nov <- nov_IM[, idx]
re_nov[] <- lapply(re_nov, function(x) if(is.numeric(x)) round(x, 3) else x)

print(re_nov)

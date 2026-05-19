################################################################################
# 
# Purpose: produce Supplementary Table D.2
#      perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the Lake data  
# Input: load 'result/lpdata.RData' file
# Output: print Supplementary Table D.2
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/IMC.r', sep=''))

load(paste(MS.PATH, 'result/lpdata.RData',sep=''))
fitN <- IM.MVN(estN$para.est, estN$EST, Yc)
fitT <- IM.MVT(estT$para.est, estT$EST, Yc)
fitS <- IM.MSL(estS$para.est, estS$EST, Yc)

fitNC <- IM.MVNC(estNC$para.est, estNC$EST, Yc, cen)
fitTC <- IM.MVTC(estTC$para.est, estTC$EST, Yc, cen)
fitSC <- IM.MSLC(estSC$para.est, estSC$EST, Yc, cen, ITER=20, per.iter=2)

# EST
lp_EST <- cbind(c(fitN$EST,0), fitT$EST, fitS$EST,
                c(fitNC$EST,0), fitTC$EST, fitSC$EST)
lp_EST_data <- as.data.frame(lp_EST)
colnames(lp_EST_data) <- c("MVN.EST", "MVT.EST", "MSL.EST", "MVNC.EST", "MVTC.EST", "MSLC.EST")
row.names(lp_EST_data) <- c("mu1", "mu2", "mu3", "mu4", "mu5", 
                            "sigma11", "sigma21", "sigma22",
                            "sigma31", "sigma32", "sigma33", 
                            "sigma41", "sigma42", "sigma43", "sigma44",
                            "sigma51", "sigma52", "sigma53", "sigma54", "sigma55", "nu")
# SD
lp_SD <- cbind(c(fitN$SD,0), fitT$SD, fitS$SD,
               c(fitNC$SD,0), fitTC$SD, fitSC$SD)
lp_SD_data <- as.data.frame(lp_SD)
colnames(lp_SD_data) <- c("MVN.SE", "MVT.SE", "MSL.SE", "MVNC.SE", "MVTC.SE", "MSLC.SE")
row.names(lp_SD_data) <- c("sd.mu1", "sd.mu2", "sd.mu3", "sd.mu4", "sd.mu5", 
                           "sd.sigma11", "sd.sigma21", "sd.sigma22",
                           "sd.sigma31", "sd.sigma32", "sd.sigma33", 
                           "sd.sigma41", "sd.sigma42", "sd.sigma43", "sd.sigma44",
                           "sd.sigma51", "sd.sigma52", "sd.sigma53", "sd.sigma54", "sd.sigma55", "sd.nu")
# combined and sort
lp_IM <- cbind(lp_EST_data, lp_SD_data)
n_vars <- ncol(lp_EST_data)
idx <- c(rbind(1:n_vars, (n_vars + 1):(2 * n_vars)))
re_lp <- lp_IM[, idx]
re_lp[] <- lapply(re_lp, function(x) if(is.numeric(x)) round(x, 3) else x)

print(re_lp)
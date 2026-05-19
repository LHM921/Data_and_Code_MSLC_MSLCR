################################################################################
# 
# Purpose: produce Supplementary Table D.3
#      perform model fitting of the MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models to the Lake data  
# Input: load 'result/lake.RData' file
# Output: print Supplementary Table D.3
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/IMCR.r', sep=''))

load(paste(MS.PATH, 'result/lake.RData',sep=''))
fitNR <- IM.MVNR(estNR$para.est, estNR$EST, Yc, X)
fitTR <- IM.MVTR(estTR$para.est, estTR$EST, Yc, X)
fitSR <- IM.MSLR(estSR$para.est, estSR$EST, Yc, X)

fitNCR <- IM.MVNCR(estNCR$para.est, estNCR$EST, Yc, X, cen)
fitTCR <- IM.MVTCR(estTCR$para.est, estTCR$EST, Yc, X, cen)
fitSCR <- IM.MSLCR(estSCR$para.est, estSCR$EST, Yc, X, cen, ITER=20, per.iter=2)

# EST
lake_EST <- cbind(c(fitNR$EST,0), fitTR$EST, fitSR$EST, 
                  c(fitNCR$EST,0), fitTCR$EST, fitSCR$EST)
lake_EST_data <- as.data.frame(lake_EST)
colnames(lake_EST_data) <- c("MVNR.EST", "MVTR.EST", "MSLR.EST", "MVNCR.EST", "MVTCR.EST", "MSLCR.EST")
row.names(lake_EST_data) <- c("beta11", "beta12", "beta13", "beta14",
                              "beta21", "beta22", "beta23", "beta24",
                              "beta31", "beta32", "beta33", "beta34",
                              "sigma11", "sigma21", "sigma22",
                              "sigma31", "sigma32", "sigma33", "nu")
# SD
lake_SD <- cbind(c(fitNR$SD,0), fitTR$SD, fitSR$SD,
                 c(fitNCR$SD,0), fitTCR$SD, fitSCR$SD)
lake_SD_data <- as.data.frame(lake_SD)
colnames(lake_SD_data) <- c("MVNR.SE", "MVTR.SE", "MSLR.SE", "MVNCR.SE", "MVTCR.SE", "MSLCR.SE")
row.names(lake_SD_data) <- c("sd.beta11", "sd.beta12", "sd.beta13", "sd.beta14",
                             "sd.beta21", "sd.beta22", "sd.beta23", "sd.beta24",
                             "sd.beta31", "sd.beta32", "sd.beta33", "sd.beta34",
                             "sd.sigma11", "sd.sigma21", "sd.sigma22",
                             "sd.sigma31", "sd.sigma32", "sd.sigma33", "sd.nu")
# combined and sort
lake_IM <- cbind(lake_EST_data, lake_SD_data)
n_vars <- ncol(lake_EST_data)
idx <- c(rbind(1:n_vars, (n_vars + 1):(2 * n_vars)))
re_lake <- lake_IM[, idx]
re_lake[] <- lapply(re_lake, function(x) if(is.numeric(x)) round(x, 3) else x)


print(re_lake)
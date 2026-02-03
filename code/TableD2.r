################################################################################
# 
# Purpose: produce Supplementary Table D.2
#      perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the Lake data  
# Input: load 'result/lpdata.RData' file
# Output: print Supplementary Table D.2
# 
################################################################################

set.seed(1)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
load(paste(MS.PATH, 'result/lpdata.RData',sep=''))

Table1b = cbind(t(estN$IM$mu.hat), t(estT$IM$mu.hat), t(estS$IM$mu.hat),
                t(estNC$IM$mu.hat), t(estTC$IM$mu.hat), t(estSC$IM$mu.hat))
colnames(Table1b) <- c("MVN_mu", "MVN_sd.mu", "MVT_mu", "MVT_sd.mu", "MSL_mu", "MSL_sd.mu",
                       "MVNC_mu", "MVNC_sd.mu", "MVTC_mu", "MVTC_sd.mu", "MSLC_mu", "MSLC_sd.mu")
rownames(Table1b) <- c("mu1", "mu2", "mu3", "mu4", "mu5") 
print(round(Table1b, 3))

Table1c = cbind(c(estN$EST[-c(1:ncol(Yc))],0), c(estT$EST[-c(1:ncol(Yc))]), c(estS$EST[-c(1:ncol(Yc))]),
                c(estNC$EST[-c(1:ncol(Yc))],0), c(estTC$EST[-c(1:ncol(Yc))]), c(estSC$EST[-c(1:ncol(Yc))]))
colnames(Table1c) <- c("MVN", "MVT", "MSL", "MVNC", "MVTC", "MSLC")
rownames(Table1c) <- c("sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33",
                       "sigma41", "sigma42", "sigma43", "sigma44", 
                       "sigma51", "sigma52", "sigma53", "sigma54","sigma55", "nu")
print(round(Table1c, 3))

################################################################################
# 
# Purpose: produce Supplementary Table D.3
#      perform model fitting of the MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models to the Lake data  
# Input: load 'result/lake.RData' file
# Output: print Supplementary Table D.3
# 
################################################################################

set.seed(18)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
load(paste(MS.PATH, 'result/lake.RData',sep=''))

Table1b = cbind(t(estNR$IM$beta.hat), t(estTR$IM$beta.hat), t(estSR$IM$beta.hat),
                t(estNCR$IM$beta.hat), t(estTCR$IM$beta.hat), t(estSCR$IM$beta.hat))
colnames(Table1b) <- c("MVNR_beta", "MVNR_sd.beta", "MVTR_beta", "MVTR_sd.beta","MSLR_beta", "MSLR_sd.beta",
                       "MVNCR_beta", "MVNCR_sd.beta", "MVTCR_beta", "MVTCR_sd.beta","MSLCR_beta", "MSLCR_sd.beta")
rownames(Table1b) <- c("beta11", "beta12", "beta13", "beta14",
                       "beta21", "beta22", "beta23", "beta24",
                       "beta31", "beta32", "beta33", "beta34")
print(round(Table1b, 3))

Table1c = cbind(c(estNR$EST[-c(1:ncol(X))],0), c(estTR$EST[-c(1:ncol(X))]), c(estSR$EST[-c(1:ncol(X))]),
                c(estNCR$EST[-c(1:ncol(X))],0), c(estTCR$EST[-c(1:ncol(X))]), c(estSCR$EST[-c(1:ncol(X))]))
colnames(Table1c) <- c("MVNR", "MVTR", "MSLR", "MVNCR", "MVTCR", "MSLCR")
rownames(Table1c) <- c("sigma11", "sigma21", "sigma22", 
                       "sigma31", "sigma32", "sigma33", "nu")
print(round(Table1c, 3))
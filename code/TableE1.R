################################################################################
# 
# Purpose: produce Supplementary Table E.1
#      perform simulation results for the MVN, MVT, MSL, MVNC, MVTC and MSLC models  
# Input: read the 'Sest_all.csv' and 'SSE_all.csv' for 100 Monte Carlo datasets
# Output: print Supplementary Table E.1: Averages of estimated values (EST), 
#                                        Monte Carlo standared deivations (MCSD), 
#                                        and information-based standard errors (IMSE)
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")

# Read All Simulation results -------------- 
# EST -------------- 
Sest <- read.csv(paste(MS.PATH,'result/SIM/Sest_all.csv', sep=''))
Sest$Model <- factor(Sest$Model, levels = c('MVN', 'MVT', 'MSL', 'MVNC', 'MVTC', 'MSLC'))
est_cols <- grep("^EST", names(Sest), value = TRUE)[1:9]
table_EST <- aggregate(Sest[, est_cols],
                       by = list(N = Sest$N, Crate = Sest$Crate, Model = Sest$Model, nu = Sest$nu),
                       FUN = function(x) round(mean(x), 3))

EST_MSL_0 <- table_EST[table_EST$Crate == 0 & table_EST$Model == "MSL", ]
EST_MSLC_5 <- table_EST[table_EST$Crate == 0.05 & table_EST$Model == "MSLC",]
EST_MSLC_15 <- table_EST[table_EST$Crate == 0.15 & table_EST$Model == "MSLC", ]

colnames(EST_MSL_0) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                         "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
colnames(EST_MSLC_5) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                          "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
colnames(EST_MSLC_15) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                           "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
rownames(EST_MSL_0) <- NULL ; rownames(EST_MSLC_5) <- NULL ; rownames(EST_MSLC_15) <- NULL

print(list(EST = EST_MSL_0, EST = EST_MSLC_5, EST = EST_MSLC_15))

# MCSD -------------- 
Sest <- read.csv(paste(MS.PATH,'result/SIM/Sest_all.csv', sep=''))
Sest$Model <- factor(Sest$Model, levels = c('MVN', 'MVT', 'MSL', 'MVNC', 'MVTC', 'MSLC'))
est_cols <- grep("^EST", names(Sest), value = TRUE)[1:9] 
table_mcsd <- aggregate(Sest[, est_cols],
                       by = list(N = Sest$N, Crate = Sest$Crate, Model = Sest$Model, nu = Sest$nu),
                       FUN = function(x) round(sd(x), 3))

mcsd_MSL_0 <- table_mcsd[table_mcsd$Crate == 0 & table_mcsd$Model == "MSL", ]
mcsd_MSLC_5 <- table_mcsd[table_mcsd$Crate == 0.05 & table_mcsd$Model == "MSLC", ]
mcsd_MSLC_15 <- table_mcsd[table_mcsd$Crate == 0.15 & table_mcsd$Model == "MSLC", ]

colnames(mcsd_MSL_0) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                         "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
colnames(mcsd_MSLC_5) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                          "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
colnames(mcsd_MSLC_15) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3", 
                           "sigma11", "sigma21", "sigma22", "sigma31", "sigma32", "sigma33")
rownames(mcsd_MSL_0) <- NULL ; rownames(mcsd_MSLC_5) <- NULL ; rownames(mcsd_MSLC_15) <- NULL

print(list(MCSD = mcsd_MSL_0, MCSD = mcsd_MSLC_5, MCSD = mcsd_MSLC_15))

# IMSE -------------- 
SSE <- read.csv(paste(MS.PATH,'result/SIM/SSE_all.csv', sep=''))
SSE$Model <- factor(SSE$Model, levels = c('MVN', 'MVT', 'MSL', 'MVNC', 'MVTC', 'MSLC'))
sse_cols <- grep("^SE", names(SSE), value = TRUE)
table_sse <- aggregate(SSE[, sse_cols],
                        by = list(N = SSE$N, Crate = SSE$Crate, Model = SSE$Model, nu = SSE$nu),
                        FUN = function(x) round(mean(x), 3))

sse_MSL_0 <- table_sse[table_sse$Crate == 0 & table_sse$Model == "MSL", ]
sse_MSLC_5 <- table_sse[table_sse$Crate == 0.05 & table_sse$Model == "MSLC", ]
sse_MSLC_15 <- table_sse[table_sse$Crate == 0.15 & table_sse$Model == "MSLC", ]

colnames(sse_MSL_0) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3")
colnames(sse_MSLC_5) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3")
colnames(sse_MSLC_15) <- c("N", "Crate", "Model", "nu", "mu1", "mu2", "mu3")
rownames(sse_MSL_0) <- NULL ; rownames(sse_MSLC_5) <- NULL ; rownames(sse_MSLC_15) <- NULL

print(list(IMSE = sse_MSL_0, IMSE = sse_MSLC_5, IMSE = sse_MSLC_15))
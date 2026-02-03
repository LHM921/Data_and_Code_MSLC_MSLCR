################################################################################
# 
# Purpose: produce Supplementary Table E.2
#      perform simulation results for the MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models  
# Input: read the 'Sest_all.csv' and 'SSE_all.csv' for 100 Monte Carlo datasets
# Output: print Supplementary Table E.2: Averages of estimated values (EST), 
#                                      Monte Carlo standared deivations (MCSD), 
#                                      Information-based standard errors (IMSE),
#                                      and Mean squared errors (MSE)
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")

# Read All Simulation results -------------- 
# EST -------------- 
Sest <- read.csv(paste(MS.PATH,'result/RSIM/Sest_all.csv', sep=''))
Sest$Model <- factor(Sest$Model, levels = c('MVNR', 'MVTR', 'MSLR', 'MVNCR', 'MVTCR', 'MSLCR'))
est_cols <- grep("^EST", names(Sest), value = TRUE)[1:3]
table_EST <- aggregate(Sest[, est_cols],
                       by = list(N = Sest$N, Crate = Sest$Crate, Model = Sest$Model, nu = Sest$nu),
                       FUN = function(x) round(mean(x), 3))
table_EST <- table_EST[order(table_EST$N, table_EST$Crate, table_EST$nu, table_EST$Model), ]

EST_0 <- table_EST[table_EST$Crate == 0 & table_EST$Model == c('MVNR', 'MVTR', 'MSLR'), ]
EST_5 <- table_EST[table_EST$Crate == 0.05 & table_EST$Model == c('MVNCR', 'MVTCR', 'MSLCR'), ]

colnames(EST_0) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
colnames(EST_5) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
rownames(EST_0) <- NULL ; rownames(EST_5) <- NULL

print(list(EST = EST_0, EST = EST_5))


# MCSD -------------- 
Sest <- read.csv(paste(MS.PATH,'result/RSIM/Sest_all.csv', sep=''))
Sest$Model <- factor(Sest$Model, levels = c('MVNR', 'MVTR', 'MSLR', 'MVNCR', 'MVTCR', 'MSLCR'))
est_cols <- grep("^EST", names(Sest), value = TRUE)[1:3] 
table_mcsd <- aggregate(Sest[, est_cols],
                        by = list(N = Sest$N, Crate = Sest$Crate, Model = Sest$Model, nu = Sest$nu),
                        FUN = function(x) round(sd(x), 3))
table_mcsd <- table_mcsd[order(table_mcsd$N, table_mcsd$Crate, table_mcsd$nu, table_mcsd$Model), ]

mcsd_0 <- table_mcsd[table_mcsd$Crate == 0 & table_mcsd$Model == c('MVNR', 'MVTR', 'MSLR'), ]
mcsd_5 <- table_mcsd[table_mcsd$Crate == 0.05 & table_mcsd$Model == c('MVNCR', 'MVTCR', 'MSLCR'), ]

colnames(mcsd_0) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
colnames(mcsd_5) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
rownames(mcsd_0) <- NULL ; rownames(mcsd_5) <- NULL

print(list(MCSD = mcsd_0, MCSD = mcsd_5))

# IMSE -------------- 
SSE <- read.csv(paste(MS.PATH,'result/RSIM/SSE_all.csv', sep=''))
SSE$Model <- factor(SSE$Model, levels = c('MVNR', 'MVTR', 'MSLR', 'MVNCR', 'MVTCR', 'MSLCR'))
sse_cols <- grep("^SE", names(SSE), value = TRUE)[1:3]
table_sse <- aggregate(SSE[, sse_cols],
                       by = list(N = SSE$N, Crate = SSE$Crate, Model = SSE$Model, nu = SSE$nu),
                       FUN = function(x) round(mean(x), 3))
table_sse <- table_sse[order(table_sse$N, table_sse$Crate, table_sse$nu, table_sse$Model), ]

imse_0 <- table_sse[table_sse$Crate == 0 & table_sse$Model == c('MVNR', 'MVTR', 'MSLR'), ]
imse_5 <- table_sse[table_sse$Crate == 0.05 & table_sse$Model == c('MVNCR', 'MVTCR', 'MSLCR'), ]

colnames(imse_0) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
colnames(imse_5) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
rownames(imse_0) <- NULL ; rownames(imse_5) <- NULL

print(list(IMSE = imse_0, IMSE = imse_5))

# MSE -------------- 
Smse <- read.csv(paste(MS.PATH,'result/RSIM/Smse_all.csv', sep=''))
Smse <- subset(Smse, select = -c(Squared_Error19, MSE))
Smse$Model <- factor(Smse$Model, levels = c('MVNR', 'MVTR', 'MSLR', 'MVNCR', 'MVTCR', 'MSLCR'))
smse_cols <- grep("^Squared_Error", names(Smse), value = TRUE)[1:3] 
table_Smse <- aggregate(Smse[, smse_cols],
                       by = list(N = Smse$N, Crate = Smse$Crate, Model = Smse$Model, nu = Smse$nu),
                       FUN = function(x) round(mean(x), 3))

table_Smse <- table_Smse[order(table_Smse$N, table_Smse$Crate, table_Smse$nu, table_Smse$Model), ]

mse_0 <- table_Smse[table_Smse$Crate == 0 & table_Smse$Model == c('MVNR', 'MVTR', 'MSLR'), ]
mse_5 <- table_Smse[table_Smse$Crate == 0.05 & table_Smse$Model == c('MVNCR', 'MVTCR', 'MSLCR'), ]

colnames(mse_0) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
colnames(mse_5) <- c("N", "Crate", "Model", "nu", "beta1", "beat2", "beta3")
rownames(mse_0) <- NULL ; rownames(mse_5) <- NULL

print(list(MSE = mse_0, MSE = mse_5))

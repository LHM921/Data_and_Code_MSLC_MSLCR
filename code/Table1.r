################################################################################
# 
# Purpose: produce Table 1
#        perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the bioassay data  
# Input: read the 'bioassay_dl.csv' and 'bioassay_cen.csv' datasets
# Source function: 'TMSLmoment.r', 'MSL.EM.r', 'MSLC.EM.r', 'MVT.EM.r' and 'MVTC.EM.r' R functions
# Output: produce 'result/bioassay_sep.RData' and print Table 1 (September part)
#         produce 'result/bioassay_nov.RData' and print Table 1 (November part)
# 
################################################################################
set.seed(6)

# September ############################################
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/MSL.EM.r', sep=''))
source(paste(MS.PATH, 'function/MSLC.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVT.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVTC.EM.r', sep=''))

# Read data
data <- read.csv(paste(MS.PATH,'data/bioassay_dl.csv', sep=''))
cen_data <- read.csv(paste(MS.PATH,'data/bioassay_cen.csv', sep=''))
data <- data[data$Month == "September", ]
cen_data <- cen_data[cen_data$Month == "September", ]
set <-  c("Chl_b","Allox","Diad")

Yc <- as.matrix((data[, set]))
cen <- as.matrix(cen_data[, set])
n = nrow(data)

# initial values
mu <- as.vector(colMeans(Yc)); p <- length(mu)
Sigma <- matrix(c(5191.758, 4447.051, 7214.674,
                  4447.051, 4011.285, 4848.212,
                  7214.674, 4848.212, 21622.464), 3, 3, byrow = T)
nu = 4; m = p + p*(p+1)/2

par.true = c(mu, Sigma[vech.posi(p)], nu)
init.par = list(mu=mu, Sigma=Sigma, nu=nu) 

# Model fitting
estN =  MSL.EM(Y=Yc, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estT =  MVT.EM(Y=Yc, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estS =  MSL.EM(Y=Yc, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1)

estNC = MSLC.EM(Yc=Yc, cen=cen, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estTC = MVTC.EM(Yc=Yc, cen=cen, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estSC = MSLC.EM(Yc=Yc, cen=cen, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1, ITER=10, per.iter=2) 

save.image(paste(MS.PATH, 'result/bioassay_sep.RData', sep=''))

# load(paste(MS.PATH, 'result/bioassay_sep.RData',sep=''))
Table1a_sep = cbind(c(estN$model.inf$m, estN$model.inf$loglik, estN$model.inf$aic, estN$model.inf$bic), 
                c(estT$model.inf$m, estT$model.inf$loglik, estT$model.inf$aic, estT$model.inf$bic),
                c(estS$model.inf$m, estS$model.inf$loglik, estS$model.inf$aic, estS$model.inf$bic), 
                c(estNC$model.inf$m, estNC$model.inf$loglik, estNC$model.inf$aic, estNC$model.inf$bic), 
                c(estTC$model.inf$m, estTC$model.inf$loglik, estTC$model.inf$aic, estTC$model.inf$bic),
                c(estSC$model.inf$m, estSC$model.inf$loglik, estSC$model.inf$aic, estSC$model.inf$bic)) 
colnames(Table1a_sep) <- c("MVN", "MVT", "MSL", "MVNC", "MVTC", "MSLC")
rownames(Table1a_sep) <- c("m", "lmax", "AIC", "BIC")
# print(round(Table1a_sep, 3))

# November ############################################
set.seed(9)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/MSL.EM.r', sep=''))
source(paste(MS.PATH, 'function/MSLC.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVT.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVTC.EM.r', sep=''))

# Read data
data <- read.csv(paste(MS.PATH,'data/bioassay_dl.csv', sep=''))
cen_data <- read.csv(paste(MS.PATH,'data/bioassay_cen.csv', sep=''))
data <- data[data$Month == "November", ]
cen_data <- cen_data[cen_data$Month == "November", ]
set <-  c("Chl_b","Allox","Diad")

Yc <- as.matrix((data[, set]))
cen <- as.matrix(cen_data[, set])
n = nrow(data)

# initial values
mu <- as.vector(colMeans(Yc)); p <- length(mu)
Sigma <- cov(Yc)
nu = 1.1; m = p + p*(p+1)/2

par.true = c(mu, Sigma[vech.posi(p)], nu)
init.par = list(mu=mu, Sigma=Sigma, nu=nu) 

# Model fitting
estN =  MSL.EM(Y=Yc, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estT =  MVT.EM(Y=Yc, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estS =  MSL.EM(Y=Yc, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1)

estNC = MSLC.EM(Yc=Yc, cen=cen, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estTC = MVTC.EM(Yc=Yc, cen=cen, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estSC = MSLC.EM(Yc=Yc, cen=cen, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1, ITER=10, per.iter=2) 

save.image(paste(MS.PATH, 'result/bioassay_nov.RData', sep=''))

# load(paste(MS.PATH, 'result/bioassay_nov.RData',sep=''))
Table1a_nov = cbind(c(estN$model.inf$m, estN$model.inf$loglik, estN$model.inf$aic, estN$model.inf$bic), 
                c(estT$model.inf$m, estT$model.inf$loglik, estT$model.inf$aic, estT$model.inf$bic),
                c(estS$model.inf$m, estS$model.inf$loglik, estS$model.inf$aic, estS$model.inf$bic), 
                c(estNC$model.inf$m, estNC$model.inf$loglik, estNC$model.inf$aic, estNC$model.inf$bic), 
                c(estTC$model.inf$m, estTC$model.inf$loglik, estTC$model.inf$aic, estTC$model.inf$bic),
                c(estSC$model.inf$m, estSC$model.inf$loglik, estSC$model.inf$aic, estSC$model.inf$bic)) 
colnames(Table1a_nov) <- c("MVN", "MVT", "MSL", "MVNC", "MVTC", "MSLC")
rownames(Table1a_nov) <- c("m", "lmax", "AIC", "BIC")


################################
# print Table 1 (September part) 
print(round(Table1a_sep, 3))
# print Table 1 (November part)
print(round(Table1a_nov, 3))

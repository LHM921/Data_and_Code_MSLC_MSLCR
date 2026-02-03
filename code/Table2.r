################################################################################
# 
# Purpose: produce Table 2
#      perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the Lake data  
# Input: read the 'data/LPdata.csv' datasets
# Source function: 'TMSLmoment.r', 'MSL.EM.r', 'MSLC.EM.r', 'MVT.EM.r' and 'MVTC.EM.r' R functions
# Output: produce 'result/lpdata.RData' and print Table 2
# 
################################################################################

set.seed(1)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/MSL.EM.r', sep=''))
source(paste(MS.PATH, 'function/MSLC.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVT.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVTC.EM.r', sep=''))

# Read data
LPdata <- read.csv(paste(MS.PATH,'data/LPdata.csv', sep=''))
dataf <- function(df) {
  updated <- lapply(df, function(col) {
    if (is.numeric(col)) {
      min_val <- min(col[col != 0], na.rm = TRUE)
      col[col == 0] <- min_val
    }
    col
  })
  cen_ind <- lapply(df, function(col) {
    if (is.numeric(col)) ifelse(col == 0, 1, 0) else rep(0, length(col))
  })
  list(
    y = as.data.frame(updated),       
    CEN = as.data.frame(cen_ind) 
  )
}
result <- dataf(LPdata)
set <- c("Cr", "Y", "La", "Ce", "Ba")
Yc <- as.matrix(result$y[, set])
cen <- as.matrix(result$CEN[, set])
n = nrow(Yc)

# initial values
mu <- as.vector(colMeans(Yc)); p <- length(mu)
generate_sigma <- function(p) {
  A <- matrix(runif(p^2, min=11, max=30), nrow=p)  
  sigma <- A %*% t(A)                            
  diag(sigma) <- diag(sigma) + 1                
  return(sigma)
}
Sigma <- generate_sigma(p)
nu = 2; m = p + p*(p+1)/2

par.true = c(mu, Sigma[vech.posi(p)], nu)
init.par = list(mu=mu, Sigma=Sigma, nu=nu) 

# Model fitting
estN =  MSL.EM(Y=Yc, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estT =  MVT.EM(Y=Yc, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estS =  MSL.EM(Y=Yc, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1)

estNC = MSLC.EM(Yc=Yc, cen=cen, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estTC = MVTC.EM(Yc=Yc, cen=cen, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estSC = MSLC.EM(Yc=Yc, cen=cen, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1, ITER=10, per.iter=2) 

save.image(paste(MS.PATH, 'result/lpdata.RData', sep=''))

# load(paste(MS.PATH, 'result/lpdata.RData',sep=''))
Table1a = cbind(c(estN$model.inf$m, estN$model.inf$loglik, estN$model.inf$aic, estN$model.inf$bic), 
                c(estT$model.inf$m, estT$model.inf$loglik, estT$model.inf$aic, estT$model.inf$bic),
                c(estS$model.inf$m, estS$model.inf$loglik, estS$model.inf$aic, estS$model.inf$bic), 
                c(estNC$model.inf$m, estNC$model.inf$loglik, estNC$model.inf$aic, estNC$model.inf$bic), 
                c(estTC$model.inf$m, estTC$model.inf$loglik, estTC$model.inf$aic, estTC$model.inf$bic),
                c(estSC$model.inf$m, estSC$model.inf$loglik, estSC$model.inf$aic, estSC$model.inf$bic)) 
colnames(Table1a) <- c("MVN", "MVT", "MSL", "MVNC", "MVTC", "MSLC")
rownames(Table1a) <- c("m", "lmax", "AIC", "BIC")
print(round(Table1a, 3))
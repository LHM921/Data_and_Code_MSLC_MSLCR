################################################################################
# 
# Purpose: produce Table 3
#      perform model fitting of the MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models to the Lake data  
# Input: read the 'data/lake_ori.csv', 'data/lake_dl.csv' and 'data/lake_cen.csv' datasets
# Source function: 'TMSLmoment.r', 'MSLR.EM.r', 'MSLC.EM.r', 'MSLCR.EM.r', 'MVTR.EM.r' and 'MVTCR.EM.r' R functions
# Output: produce 'result/lake.RData' and print Table 3
# 
################################################################################

set.seed(18)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))
source(paste(MS.PATH, 'function/MSLR.EM.r', sep=''))    
source(paste(MS.PATH, 'function/MSLC.EM.r', sep=''))
source(paste(MS.PATH, 'function/MSLCR.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVTR.EM.r', sep=''))
source(paste(MS.PATH, 'function/MVTCR.EM.r', sep=''))

# Read data
df <- read.csv(paste(MS.PATH,'data/lake_ori.csv', sep=''))
dc <- read.csv(paste(MS.PATH,'data/lake_dl.csv', sep=''))
cen <- read.csv(paste(MS.PATH,'data/lake_cen.csv', sep=''))

setx <- c("Lat","Long","Depth_Site", "Depth_Smp")
sety <- c("TDP","PC","PN")
yori <- as.matrix(df[,sety])
y <- scale(yori)
Yc <- as.matrix(dc[, sety])
Yc <- scale(Yc)

p = ncol(Yc) ; n = nrow(Yc)
cen = as.matrix(cen[,sety]) 
X <-  as.matrix(scale(df[, setx]))
X = kronecker(diag(p),X)

# initial value
beta = c(solve(t(X)%*%X) %*% t(X) %*% as.vector(t(Yc)))
generate_sigma <- function(p) {
  A <- matrix(runif(p^2, min=0, max=1), nrow=p)  
  sigma <- A %*% t(A)                            
  diag(sigma) <- diag(sigma) +0.5
  return(sigma)
}
Sigma <- generate_sigma(p)
m = length(beta) + p*(p+1)/2
nu = 5
para = list(beta=beta, Sigma=Sigma, nu=nu)
par.true = c(beta, Sigma[vech.posi(p)], nu)
init.par = list(beta=beta, Sigma=Sigma, nu=nu) 

# Model fitting
estNR = MVTR.EM(Y=Yc, X, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)   
estTR = MVTR.EM(Y=Yc, X, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estSR = MSLR.EM(Y=Yc, X, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1) 

estNCR = MSLCR.EM(Yc=Yc, X, cen=cen, distr='MVN', init.par=init.par, tol=1e-3, max.iter=100, per=1)  
estTCR = MVTCR.EM(Yc=Yc, X, cen=cen, distr='MVT', init.par=init.par, tol=1e-3, max.iter=100, per=1)
estSCR = MSLCR.EM(Yc=Yc, X, cen=cen, distr='MSL', init.par=init.par, tol=1e-3, max.iter=100, per=1, ITER=10, per.iter=2)  

save.image(paste(MS.PATH, 'result/lake.RData', sep=''))

# load(paste(MS.PATH, 'result/lake.RData',sep=''))
Table1a = cbind(c(estNR$model.inf$m, estNR$model.inf$loglik, estNR$model.inf$aic, estNR$model.inf$bic), 
                c(estTR$model.inf$m, estTR$model.inf$loglik, estTR$model.inf$aic, estTR$model.inf$bic),
                c(estSR$model.inf$m, estSR$model.inf$loglik, estSR$model.inf$aic, estSR$model.inf$bic), 
                c(estNCR$model.inf$m, estNCR$model.inf$loglik, estNCR$model.inf$aic, estNCR$model.inf$bic), 
                c(estTCR$model.inf$m, estTCR$model.inf$loglik, estTCR$model.inf$aic, estTCR$model.inf$bic),
                c(estSCR$model.inf$m, estSCR$model.inf$loglik, estSCR$model.inf$aic, estSCR$model.inf$bic)) 
colnames(Table1a) <- c("MVNR", "MVTR", "MSLR", "MVNCR", "MVTCR", "MSLCR")
rownames(Table1a) <- c("m", "lmax", "AIC", "BIC")
print(round(Table1a, 3))
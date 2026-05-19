################################################################################ 
# Purpose: produce Figure 2
#   draw various scatter diagrams and marginal histograms of bivariate slash random samples
# Source function: 'TMSLmoment.r' R functions
# Output: produce 'result/fig2.pdf'
################################################################################

set.seed(112)
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))

##### 【fig2-1】under no truncation ##### 
################################################################ nu=4
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_2 = 4
ubd_2 = 2
a.low_2=rep(-Inf, 2); a.upp_2=rep(Inf, 2)
n=500
x1 = seq(-15, 15, 0.3)
len1 = length(x1)
d2y1 = d2y2 = d2y3 = matrix(NA, len1, len1)
d2y21 = d2y22 = d2y23 = d2y24 =d2y25 = d2y26 = numeric(len1)
# TMVN
for(i in 1: len1)for(j in 1: len1) d2y1[i,j] = dmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN') 
for(j in 1: len1){
  d2y21[j] = dmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVN')    
  d2y22[j] = dmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVN')    
}
# TMSL
y_2 = rtmsl(n=n, mu_2, Sigma_2, nu_2, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len1)for(j in 1: len1) d2y2[i,j] = dmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MSL') 
for(j in 1: len1){
  d2y23[j] = dmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MSL')    
  d2y24[j] = dmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MSL')   
}
# TMVT
for(i in 1: len1)for(j in 1: len1) d2y3[i,j] = dmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MVT') 
for(j in 1: len1){
  d2y25[j] = dmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVT')    
  d2y26[j] = dmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVT')    
}

y2=hist(y_2[,2], nclass=35, plot=FALSE)

################################################################ nu=10
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_4 = 10
ubd_2 = 2 
a.low_2=rep(-Inf, 2); a.upp_2=rep(Inf, 2)
n=500
x2 = seq(-15, 15, 0.3)
len2 = length(x2)
d4y1 = d4y2 = d4y3 = matrix(NA, len2, len2)
d4y21 = d4y22 = d4y23 = d4y24 =d4y25 = d4y26 = numeric(len2)
# TMVN
for(i in 1: len2)for(j in 1: len2) d4y1[i,j] = dmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN') 
for(j in 1: len2){
  d4y21[j] = dmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVN')    
  d4y22[j] = dmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVN')    
}
# TMSL
y_4 = rtmsl(n=n, mu_2, Sigma_2, nu_4, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len2)for(j in 1: len2) d4y2[i,j] = dmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MSL') 
for(j in 1: len2){
  d4y23[j] = dmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MSL')    
  d4y24[j] = dmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MSL')    
}
# TMVT
for(i in 1: len2)for(j in 1: len2) d4y3[i,j] = dmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MVT') 
for(j in 1: len2){
  d4y25[j] = dmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVT')    
  d4y26[j] = dmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVT')    
}

y4=hist(y_4[,2], nclass=22, plot=FALSE)

################################################################ Plot
cairo_pdf(file = file.path(MS.PATH, "result", "fig2.pdf"), width=10.5, height=7)

layout(matrix(c(0,0,0,0,0,0,0,0,1,0,7,0,13,0,19,0,3,4,9,10,15,16,21,22,2,0,8,0,14,0,20,0,5,6,11,12,17,18,23,24,0,0,0,0,0,0,0,0),6, 8, byrow=TRUE),heights=c(0.5,1,2,1,2,0.5),widths = c(2.7,1,2.7,1,2.7,1,2.7,1,2.7,1))

par(mar = c(0, 4.1, 1, 0))
hist(y_2[,1], nclass=35, prob=T, ylim=c(0, 0.4), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=4")),las="1",xlim=c(-8,10))
lines(x1, d2y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x1, d2y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x1, d2y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd = 0.75)

par(mar = c(0, 4.1, 1, 0))
hist(y_4[,1], nclass=21, prob=T, ylim=c(0, 0.4), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=10")),las="1",xlim=c(-7,9))
lines(x2, d4y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x2, d4y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x2, d4y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd = 0.75)


par(mar = c(4.6, 4.1, 0, 0))
plot(y_2, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",ylim=c(-5,9),xlim=c(-8,10))
contour(x1,x2,d2y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd=0.5)

par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n" ,xlim = c(0,0.4),ylim=c(-5,9),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y2$breaks[1:(length(y2$breaks) - 1)], y2$density, y2$breaks[2:length(y2$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d2y22,x1, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d2y24,x1, lty=1, col="red",lwd = 0.75)#MSL
lines(d2y26,x1, lty=4, col="green4",lwd = 0.65)#MVT

par(mar = c(4.6, 4.1, 0, 0))
plot(y_4, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",ylim=c(-5,8),xlim=c(-7,9))
contour(x2,x2,d4y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)


par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0,0.4),ylim=c(-5,8),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y4$breaks[1:(length(y4$breaks) - 1)], y4$density, y4$breaks[2:length(y4$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d4y22,x2, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d4y24,x2, lty=1, col="red",lwd = 0.75)#MSL
lines(d4y26,x2, lty=4, col="green4",lwd = 0.65)#MVT


##### 【fig2-2】under right truncation ##### 
################################################################ nu=4
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_2 = 4
ubd_2 = 2
a.low_2=c(-Inf,-Inf); a.upp_2=c(2,4)
n=500
x1 = seq(-15, 15, 0.3)
len1 = length(x1)
d2y1 = d2y2 = d2y3 = matrix(NA, len1, len1)
d2y21 = d2y22 = d2y23 = d2y24 =d2y25 = d2y26 = numeric(len1)
# TMVN
for(i in 1: len1)for(j in 1: len1) d2y1[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2)
for(j in 1: len1){
  d2y21[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])      
  d2y22[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])      
}
# TMSL
y_2 = rtmsl(n=n, mu_2, Sigma_2, nu_2, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len1)for(j in 1: len1) d2y2[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y23[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])  
  d2y24[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])   
}
# TMVT
for(i in 1: len1)for(j in 1: len1) d2y3[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MVT', a.low=a.low_2, a.upp=a.upp_2)
for(j in 1: len1){
  d2y25[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])     
  d2y26[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])      
}

y2=hist(y_2[,2], nclass=25, plot=FALSE)

################################################################ nu=10
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_4 = 10
ubd_2 = 2
a.low_2=c(-Inf,-Inf); a.upp_2=c(2,4)
n=500
x2 = seq(-15, 15, 0.3)
len2 = length(x2)
d4y1 = d4y2 = d4y3 = matrix(NA, len2, len2)
d4y21 = d4y22 = d4y23 = d4y24 =d4y25 = d4y26 = numeric(len2)
# TMVN
for(i in 1: len2)for(j in 1: len2) d4y1[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y21[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y22[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])   
}
# TMSL
y_4 = rtmsl(n=n, mu_2, Sigma_2, nu_4, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len2)for(j in 1: len2) d4y2[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y23[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y24[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMVT
for(i in 1: len2)for(j in 1: len2) d4y3[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MVT', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y25[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y26[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}

y4=hist(y_4[,2], nclass=20, plot=FALSE)


################################################################ Plot

par(mar = c(0, 4.1, 1, 0))
hist(y_2[,1], nclass=30, prob=T, ylim=c(0, 0.5), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=4")),las="1",xlim=c(-9,2.2))
lines(x1, d2y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x1, d2y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x1, d2y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)

par(mar = c(0, 4.1, 1, 0))
hist(y_4[,1], nclass=20, prob=T, ylim=c(0, 0.5), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=10")),las="1",xlim=c(-10,2.2))
lines(x2, d4y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x2, d4y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x2, d4y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)


par(mar = c(4.6, 4.1, 0, 0))
plot(y_2, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",xlim=c(-9,2.2),ylim=c(-5,4.5))
contour(x1,x2,d2y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)


par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0, 0.5), ylim=c(-5,4.5),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y2$breaks[1:(length(y2$breaks) - 1)], y2$density, y2$breaks[2:length(y2$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d2y22,x1, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d2y24,x1, lty=1, col="red",lwd = 0.75)#MSL
lines(d2y26,x1, lty=4, col="green4",lwd = 0.65)#MVT

par(mar = c(4.6, 4.1, 0, 0))
plot(y_4, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",xlim=c(-10,2.2),ylim=c(-4,4.6))
contour(x2,x2,d4y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)

par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0, 0.5), ylim=c(-4,4.6),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y4$breaks[1:(length(y4$breaks) - 1)], y4$density, y4$breaks[2:length(y4$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d4y22,x2, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d4y24,x2, lty=1, col="red",lwd = 0.75)#MSL
lines(d4y26,x2, lty=4, col="green4",lwd = 0.65)#MVT


##### 【fig2-3】 under left truncation ##### 
################################################################ nu=4
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_2 = 4
ubd_2 = 2
a.low_2=c(-ubd_2,1); a.upp_2=c(Inf,Inf)
n=500
x1 = seq(-15, 15, 0.3)
len1 = length(x1)
d2y1 = d2y2 = d2y3 = matrix(NA, len1, len1)
d2y21 = d2y22 = d2y23 = d2y24 =d2y25 = d2y26 = numeric(len1)
# TMVN
for(i in 1: len1)for(j in 1: len1) d2y1[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y21[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d2y22[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMSL
y_2 = rtmsl(n=n, mu_2, Sigma_2, nu_2, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len1)for(j in 1: len1) d2y2[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y23[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d2y24[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMVT
for(i in 1: len1)for(j in 1: len1) d2y3[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MVT', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y25[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d2y26[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}

y2=hist(y_2[,2], nclass=25, plot=FALSE)

################################################################ nu=10
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_4 = 10
ubd_2 = 2 
a.low_2=c(-ubd_2,1); a.upp_2=c(Inf,Inf)
n=500
x2 = seq(-15, 15, 0.3)
len2 = length(x2)
d4y1 = d4y2 = d4y3 = matrix(NA, len2, len2)
d4y21 = d4y22 = d4y23 = d4y24 =d4y25 = d4y26 = numeric(len2)
# TMVN
for(i in 1: len2)for(j in 1: len2) d4y1[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y21[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y22[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])   
}
# TMSL
y_4 = rtmsl(n=n, mu_2, Sigma_2, nu_4, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len2)for(j in 1: len2) d4y2[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y23[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y24[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMVT
for(i in 1: len2)for(j in 1: len2) d4y3[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MVT', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y25[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y26[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}

y4=hist(y_4[,2], nclass=20, plot=FALSE)


################################################################ Plot
par(mar = c(0, 4.1, 1, 0))
hist(y_2[,1], nclass=30, prob=T, ylim=c(0, 0.4), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=4")),las="1",xlim=c(-2.3,10))
lines(x1, d2y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x1, d2y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x1, d2y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)

par(mar = c(0, 4.1, 1, 0))
hist(y_4[,1], nclass=25, prob=T, ylim=c(0, 0.4), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=10")),las="1",xlim=c(-2.3,10))
lines(x2, d4y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x2, d4y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x2, d4y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA, c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)


par(mar = c(4.6, 4.1, 0, 0))
plot(y_2, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",xlim=c(-2.3,10),ylim=c(0.3,10.6))
contour(x1,x2,d2y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)


par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0, 0.4), ylim=c(0.3,10.6),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y2$breaks[1:(length(y2$breaks) - 1)], y2$density, y2$breaks[2:length(y2$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d2y22,x1, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d2y24,x1, lty=1, col="red",lwd = 0.75)#MSL
lines(d2y26,x1, lty=4, col="green4",lwd = 0.65)#MVT

par(mar = c(4.6, 4.1, 0, 0))
plot(y_4, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46",xlim=c(-2.3,10),ylim=c(0.3,10.6))
contour(x2,x2,d4y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)

par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0,0.4), ylim=c(0.3,10.6),bty="n", yaxt='n',xlab="",las="1") 
rect(0, y4$breaks[1:(length(y4$breaks) - 1)], y4$density, y4$breaks[2:length(y4$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d4y22,x2, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d4y24,x2, lty=1, col="red",lwd = 0.75)#MSL
lines(d4y26,x2, lty=4, col="green4",lwd = 0.65)#MVT


##### 【fig2-4】under double truncation ##### 
################################################################ nu=4
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_2 = 4
ubd_2 = 2
a.low_2=c(-2,1); a.upp_2=c(2,4)
n=500
x1 = seq(-15, 15, 0.3)
len1 = length(x1)
d2y1 = d2y2 = d2y3 = matrix(NA, len1, len1)
d2y21 = d2y22 = d2y23 = d2y24 =d2y25 = d2y26 = numeric(len1)
# TMVN
for(i in 1: len1)for(j in 1: len1) d2y1[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y21[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])   
  d2y22[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMSL
y_2 = rtmsl(n=n, mu_2, Sigma_2, nu_2, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len1)for(j in 1: len1) d2y2[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y23[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d2y24[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])   
}
# TMVT
for(i in 1: len1)for(j in 1: len1) d2y3[i,j] = dtmsl(x=c(x1[i], x1[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_2, distr='MVT', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len1){
  d2y25[j] = dtmsl(x=x1[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_2, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d2y26[j] = dtmsl(x=x1[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_2, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}

y2=hist(y_2[,2], nclass=20, plot=FALSE)

################################################################ nu=10
mu=c(1,2,3)
Sigma=1.5*matrix(c(3,2,1,2,2,1,1,1,1), 3,3)
mu_2 = mu[c(1,3)]
Sigma_2=Sigma[c(1,3),c(1,3)]

nu_4 = 10
ubd_2 = 2
a.low_2=c(-2,1); a.upp_2=c(2,4)
n=500
x2 = seq(-15, 15, 0.3)
len2 = length(x2)
d4y1 = d4y2 = d4y3 = matrix(NA, len2, len2)
d4y21 = d4y22 = d4y23 = d4y24 =d4y25 = d4y26 = numeric(len2)
# TMVN
for(i in 1: len2)for(j in 1: len2) d4y1[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, distr='MVN', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y21[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVN', a.low=a.low_2[1], a.upp=a.upp_2[1])   
  d4y22[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVN', a.low=a.low_2[2], a.upp=a.upp_2[2])   
}
# TMSL
y_4 = rtmsl(n=n, mu_2, Sigma_2, nu_4, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
for(i in 1: len2)for(j in 1: len2) d4y2[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MSL', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y23[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MSL', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y24[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MSL', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}
# TMVT
for(i in 1: len2)for(j in 1: len2) d4y3[i,j] = dtmsl(x=c(x2[i], x2[j]), mu=mu_2, Sigma=Sigma_2, nu=nu_4, distr='MVT', a.low=a.low_2, a.upp=a.upp_2) 
for(j in 1: len2){
  d4y25[j] = dtmsl(x=x2[j], mu=mu_2[1], Sigma=Sigma_2[1,1], nu=nu_4, distr='MVT', a.low=a.low_2[1], a.upp=a.upp_2[1])    
  d4y26[j] = dtmsl(x=x2[j], mu=mu_2[2], Sigma=Sigma_2[2,2], nu=nu_4, distr='MVT', a.low=a.low_2[2], a.upp=a.upp_2[2])    
}

y4=hist(y_4[,2], nclass=20, plot=FALSE)


################################################################ Plot
par(mar = c(0, 4.1, 1, 0))
hist(y_2[,1], nclass=20, prob=T, ylim=c(0,0.5), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=4")),las="1", xlim=c(-2.5, 2.5))
lines(x1, d2y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x1, d2y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x1, d2y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA,c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)

par(mar = c(0, 4.1, 1, 0))
hist(y_4[,1], nclass=20, prob=T, ylim=c(0, 0.5), xaxt='n', col="#96FED1", border = "#D0D0D0",main=expression(paste(nu,"=10")),las="1", xlim=c(-2.5, 2.5))
lines(x2, d4y21, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(x2, d4y23, lty=1, col="red",lwd = 0.75)#MSL
lines(x2, d4y25, lty=4, col="green4",lwd = 0.65)#MVT
legend(par('usr')[2], par('usr')[4], xpd=NA,c('MSL','MVN','MVT'), lty=c(1,2,4), col=c("red","blue","green4"), bty='n',lwd=0.75)


par(mar = c(4.6, 4.1, 0, 0))
plot(y_2, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46", xlim=c(-2.5, 2.5), ylim=c(0.5, 4.5))
contour(x1,x2,d2y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x1,x2,d2y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)

par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0, 0.5), ylim=c(0.5, 4.5),bty="n",xlab="",las="1", yaxt='n') 
rect(0, y2$breaks[1:(length(y2$breaks) - 1)], y2$density, y2$breaks[2:length(y2$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d2y22,x1, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d2y24,x1, lty=1, col="red",lwd = 0.75)#MSL
lines(d2y26,x1, lty=4, col="green4",lwd = 0.65)#MVT


par(mar = c(4.6, 4.1, 0, 0))
plot(y_4, pch=18, cex=0.3, xlab='Variable 1', ylab='Variable 3',las="1",col="gray46", xlim=c(-2.5, 2.5), ylim=c(0.5, 4.5))
contour(x2,x2,d4y1, lty=2, col="blue", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y2, lty=1, col="red", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)
contour(x2,x2,d4y3, lty=4, col="green4", nlevels=7, add=T,levels=c(0.02,0.008,0.0008),lwd = 0.5)

par(mar = c(4.6, 0, 0,0.1))
plot(NULL, type = "n", xlim = c(0, 0.5), ylim=c(0.5, 4.5),bty="n",xlab="",las="1", yaxt='n') 
rect(0, y4$breaks[1:(length(y4$breaks) - 1)], y4$density, y4$breaks[2:length(y4$breaks)],col="#96FED1", border = "#D0D0D0") 
lines(d4y22,x2, lty=2, col="blue" ,lwd = 0.5)#MVN
lines(d4y24,x2, lty=1, col="red",lwd = 0.75)#MSL
lines(d4y26,x2, lty=4, col="green4",lwd = 0.65)#MVT



dev.off()

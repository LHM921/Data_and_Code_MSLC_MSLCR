################################################################################ 
# Purpose: produce Figure 1
#   draw 3D surfaces and contour plots of bivariate slash distributions with various truncations 
# Source function: 'TMSLmoment.r' R functions
# Output: produce 'result/fig1.eps'
################################################################################

set.seed(24)
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
source(paste(MS.PATH, 'function/TMSLmoment.r', sep=''))

##########################################
# bivariate SL density
bivSL = function(x1, x2, mu=rep(0,2), Sigma=diag(2), nu=NULL)
{
  p = 2
  sig21 = Sigma[2,1]; sig11 = Sigma[1,1]; sig22 = Sigma[2,2]
  rho = sig21 / sqrt(sig11*sig22)
  del = ( (x1-mu[1])^2/sig11 + (x2-mu[2])^2/sig22 - 2*rho*(x1-mu[1])*(x2-mu[2])/sqrt(sig11*sig22) ) / (1-rho^2)
  if(del[1] == 0){
    den = nu / ((2*pi)^(p/2)*(p+nu)*sqrt(det(Sigma)))
  } else den = nu / (2*(2*pi)^(p/2)*sqrt(det(Sigma))) * IGfn(a=p/2+nu/2, b=del/2)
  return(den)
}

# truncated bivariate SL density
bivTSL = function(x1, x2, mu=rep(0,2), Sigma=diag(2), nu=NULL, lower=rep(-Inf, 2), upper=rep(Inf, 2))
{
  cdf = cubintegrate(dmsl, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
  Tden = bivSL(x1, x2, mu=mu, Sigma=Sigma, nu=nu) / cdf
  return(Tden)
}
mu=c(1,3)
Sigma=matrix(c(4.5,1.5,1.5,1.5),2,2)
nu = 4
n=500

##############################  No Truncation
a.low_1=rep(-Inf, 2); a.upp_1=rep(Inf, 2)
Y = rtmsl(n=n, mu, Sigma, nu, distr='MSL', lower=a.low_1, upper=a.upp_1)$Y
x = seq(min(Y[,1]), max(Y[,1]), length=50)
y = seq(min(Y[,2]), max(Y[,2]), length=50)

############################## Right Truncation
ubd = 2
lbd = 4
a.low_2=c(-Inf,-Inf); a.upp_2=c(ubd,lbd)
Yt1 = rtmsl(n=n, mu, Sigma, nu, distr='MSL', lower=a.low_2, upper=a.upp_2)$Y
x1 = seq(min(Yt1[,1]), max(Yt1[,1]), length=50)
y1 = seq(min(Yt1[,2]), max(Yt1[,2]), length=50)

############################## Left Truncation
ubd = 2
lbd = -1
a.low_3=c(-ubd,-lbd); a.upp_3=c(Inf,Inf)
Yt2 = rtmsl(n=n, mu, Sigma, nu, distr='MSL', lower=a.low_3, upper=a.upp_3)$Y
x2 = seq(min(Yt2[,1]), max(Yt2[,1]), length=50)
y2 = seq(min(Yt2[,2]), max(Yt2[,2]), length=50)

############################## Double Truncation
ubd = 2
lbd = 4
a.low_4=c(-2,1); a.upp_4=c(ubd,lbd)
Yt3 = rtmsl(n=n, mu, Sigma, nu, distr='MSL', lower=a.low_4, upper=a.upp_4)$Y
x3 = seq(min(Yt3[,1]), max(Yt3[,1]), length=50)
y3 = seq(min(Yt3[,2]), max(Yt3[,2]), length=50)

############################## Plot
cairo_pdf(file = file.path(MS.PATH, "result", "fig1.pdf"), width=7, height=7)

layout(matrix(c(1,2,3,4),2,2,byrow = TRUE))

par(mar = c(1, 0, 1, 0))
den = outer(x, y, FUN = bivSL, mu=mu, Sigma=Sigma, nu=nu)
tran = persp(x,y,den, xlim = range(x), ylim=range(y), zlim=c(-0.01, max(den, na.rm=T)), box=T, theta=45, lwd=.5, axes=T, expand=1.2, col="#C1FFE4", phi=20, border="#7B7B7B", xlab='y1', ylab='y2', zlab='Joint density',main="No Truncation")
points(trans3d(Y[,1], Y[,2], -0.01, tran), col = '#FF0000', pch=19, cex=0.02)
clines = contourLines(x, y, den, nlevels = 10)
contourline = lapply(clines, function(contour){ lines(trans3d(contour$x,contour$y, -0.01, tran),lwd=0.45)})

par(mar = c(1, 0, 1, 0))
Tden1 = outer(x1, y1, FUN = bivTSL, mu=mu, Sigma=Sigma, nu=nu, lower=a.low_2, upper=a.upp_2)
tran1 = persp(x1,y1,Tden1, xlim = range(x1)+c(-0.5,0.5), ylim=range(y1)+c(-0.5,0.5), zlim=c(-0.01, max(Tden1, na.rm=T)), box=T, theta=80, lwd=.5, axes=T, expand=1.2, col="#C1FFE4", phi=20, border="#7B7B7B", xlab='y1', ylab='y2', zlab='Joint density',main="Right Truncation")
points(trans3d(Yt1[,1], Yt1[,2], -0.01, tran1), col = '#FF0000', pch=19, cex=0.02)
clinesT1 = contourLines(x1, y1, Tden1, nlevels = 10)
contourlineT1 = lapply(clinesT1, function(contour){ lines(trans3d(contour$x,contour$y, -0.01, tran1),lwd=0.45)})

par(mar = c(1, 0, 1, 0))
Tden2 = outer(x2, y2, FUN = bivTSL, mu=mu, Sigma=Sigma, nu=nu, lower=a.low_3, upper=a.upp_3)
tran2 = persp(x2,y2,Tden2, xlim = range(x2)+c(-0.7,3), ylim=range(y2)+c(-0.7,3), zlim=c(-0.01, max(Tden2, na.rm=T)), box=T, theta=-25, lwd=.5, axes=T, expand=1.2, col="#C1FFE4", phi=20, border="#7B7B7B", xlab='y1', ylab='y2', zlab='Joint density',main="Left Truncation")
points(trans3d(Yt2[,1], Yt2[,2], -0.01, tran2), col = '#FF0000', pch=19, cex=0.02)
clinesT2 = contourLines(x2, y2, Tden2, nlevels = 10)
contourlineT2 = lapply(clinesT2, function(contour){ lines(trans3d(contour$x,contour$y, -0.01, tran2),lwd=0.45)})

par(mar = c(1, 0, 1, 0))
Tden3 = outer(x3, y3, FUN = bivTSL, mu=mu, Sigma=Sigma, nu=nu, lower=a.low_4, upper=a.upp_4)
tran3 = persp(x3,y3,Tden3, xlim = range(x3)+c(-0.5,0.5), ylim=range(y3)+c(-0.5,0.5), zlim=c(-0.01, max(Tden3, na.rm=T)), box=T, theta=-225, lwd=.5, axes=T, expand=1.2, col="#C1FFE4", phi=20, border="#7B7B7B", xlab='y1', ylab='y2', zlab='Joint density',main="Double Truncation")
points(trans3d(Yt3[,1], Yt3[,2], -0.01, tran3), col = '#FF0000', pch=19, cex=0.02)
clinesT3 = contourLines(x3, y3, Tden3, nlevels = 10)
contourlineT3 = lapply(clinesT3, function(contour){ lines(trans3d(contour$x,contour$y, -0.01, tran3),lwd=0.45)})

dev.off()
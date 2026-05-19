################################################################################ 
# Purpose: produce Figure 3
  #     Visualization of the compositional data and corresponding analysis results
# Input: load 'result/bioassay_sep.RData' and 'result/bioassay_nov.RData' files
# Output: produce 'result/fig3.pdf'
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig3.pdf"), width=20, height=12)

library(mvtnorm)

# September -------------- 
load(paste(MS.PATH, 'result/bioassay_sep.RData',sep=''))

nf = layout(matrix(c(0,1,0,0,2,0,
                     3,6,7,15,18,19,
                     9,4,8,21,16,20,
                     10,11,5,22,23,17,
                     12,13,14,24,25,26), 5, 6, byrow = T), c(5,4,4.5,5,4,4.5), c(1,5,4,4,2.9))
# Setting -------------- 
brk = c(13,8,10)
cenrate = round(colMeans(cen)*100, 2)
skew = c(1.08,0.18,0.67)
kurt = c(3.12, 1.56, 2.23)
var.name = set
Xlab = as.list(p)
Xlab[[1]] = seq(-78,100.5,39)
Xlab[[2]] = seq(-86,70.5,36)
Xlab[[3]] = seq(9.5,520.5,101.5)
Ychat = estSC$yhat
YNchat = estNC$yhat

xli = c(-38,40,300)
alim = c(0.025,0.03,0.005)
lim = c(0.04,0.05,0.01)
at = c(-0.006,-0.0075,-0.0013)
at1 = c(-0.005,-0.0065,-0.0014)
bx = c(0.011,0.014, 0.0025)

# Month --------------
par(mar = c(0, 0, 1, 0))  
plot.new()                
text(x = 0.5, y = 0.5, labels = "September Data",
     cex = 2.5, xpd = TRUE, font=2)
plot.new()                
text(x = 0.5, y = 0.5, labels = "November Data",
     cex = 2.5, xpd = TRUE, font=2)

# Diagonal Entry -------------- 
i = 1
 par(mar=c(0,2.5,2.5,0.25))
 plot(0:1,0:1, ylim=c(at[i]-0.002, lim[i]), xlim=c(min(Ychat[,i]), max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
 h1 = hist(Yc[,i], breaks=brk[i], plot=F)
 a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
 a2 = length(h1$counts) + 1
 hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)), col=c(rep(gray(.6),a2-a1), rep(0,a1)), xlab='', ylab='Density', add=T)
 boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
 text(xli[i], alim[i], var.name[i], cex=2.9)
 legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                     paste('skew = ', skew[i], sep=''),
                     paste('kurt = ', kurt[i], sep='')), col=1, bty='n', cex=1.5)
i = 2 
 par(mar=c(0,0.25,0.25,0.25))
 plot(0:1,0:1, ylim=c(at[i]-0.004, lim[i]), xlim=c(min(Ychat[,i]), max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
 h1 = hist(Yc[,i], breaks=brk[i], plot=F)
 a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
 a2 = length(h1$counts) + 1
 hist(Yc[,i], breaks=brk[i], prob=T, xaxt='n', ylim=c(0, max(h1$density)), xlim=c(min(Ychat[,i]), max(Ychat[,i])), col=c(rep(gray(.6),a2-a1), rep(0,a1)), xlab='', ylab='Density', add=T) 
 boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
 text(xli[i], alim[i], var.name[i], cex=2.9)
 legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                    paste('skew = ', skew[i], sep=''),
                    paste('kurt = ', kurt[i], sep='')), col=1, bty='n', cex=1.5)
i = 3
 par(mar=c(0,0.25,0.25,2.5))
 plot(0:1,0:1, ylim=c(at[i]-0.001, lim[i]), xlim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))+20), type='n', xlab='', ylab='', axes=T, main='', yaxt='n', xaxt='n')
 h1 = hist(Yc[,i], breaks=brk[i], plot=F)
 a1 = sum(h1$breaks>Yc[cen[,i]==1,i][1])
 a2 = length(h1$counts) + 1
 hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)),col=c(rep(0,a2)), xlab='', ylab='Density', add=T)
 boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
 text(xli[i], alim[i], var.name[i], cex=2.9)
 legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                     paste('skew = ', skew[i], sep=''),
                     paste('kurt = ', kurt[i], sep='')), col=1, bty='n', cex=1.5)

# Upper trianglar -------------- 
 for(i in 1)for(j in 2){
   par(mar=c(0,0.25,2.5,0.25))
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
   legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
   axis(3, Xlab[[j]], cex.lab=0.9, lwd=0.45)
 }
 for(i in 1)for(j in 3){
   par(mar=c(0,0.25,2.5,2.5))
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
   legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
   axis(3, Xlab[[j]], cex.lab=0.9, lwd=0.45)
   axis(4, Xlab[[i]], cex.lab=0.9, lwd=0.45)
 }
 for(i in 2)for(j in 3){
   par(mar=c(0,0.25,0.25,2.5))
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
   legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
   axis(4, Xlab[[i]], cex.lab=0.9, lwd=0.45)
 }

# Lower trianglar --------------
m=100
mu = estSC$para$mu
Sig = estSC$para$Sigma
nu = estSC$para$nu
muN = estNC$para$mu
SigN = estNC$para$Sigma

 par(mar=c(0.25,2.5,0.25,0.25))
 for(i in 2){
  oID = which(rowSums(cen[,c(1,i)]) == 0)
  plot(Ychat[oID,1], Ychat[oID,i], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.7)
  axis(2, Xlab[[i]], cex.lab=0.9, lwd=0.45)
  cID = which(rowSums(cen[,c(1,i)]) != 0)
  points(Ychat[cID, 1], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)
  points(YNchat[cID, 1], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
  abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
  legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,i]/(sqrt(Sig[1,1]*Sig[i,i])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[1,i]/(sqrt(SigN[1,1]*SigN[i,i])),3)))), col=1, bty='n', cex=1.4)
  tmp1=range(pretty(Ychat[,1]))
  tmp2=range(pretty(Ychat[,i]))
  xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
  den=numeric(m^2)
  for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,i)], Sigma=Sig[c(1,i), c(1,i)], nu=nu, distr='MSL')
  den = matrix(den, m, m)
  denN = dmsl(xx, mu=muN[c(1,i)], Sigma=SigN[c(1,i), c(1,i)], distr='MVN')
  denN = matrix(denN, m, m)
  contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00035,0.00009,0.00001,0.0000016))
  contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00035,0.00009,0.00001,0.0000016))
 } 
 par(mar=c(0,2.5,0.25,0.25))
 oID = which(rowSums(cen[,c(1,3)]) == 0)
 plot(Ychat[oID,1], Ychat[oID,3], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,3],YNchat[,3]), max(Ychat[,3],YNchat[,3])), xaxt='n', yaxt='n', pch=16, cex=0.7)
 axis(1, Xlab[[1]], cex.lab=0.9, lwd=0.45)
 axis(2, Xlab[[3]], cex.lab=0.9, lwd=0.45)
 cID = which(rowSums(cen[,c(1,3)]) != 0)
 points(Ychat[cID, 1], Ychat[cID, 3], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, 1], YNchat[cID, 3], pch=3, col="orange", cex=0.95, lwd=0.6)
 abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,3]==1,3][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,3]/(sqrt(Sig[1,1]*Sig[3,3])),3))),
                              bquote(hat(rho)[MVNC] == .(round(SigN[1,3]/(sqrt(SigN[1,1]*SigN[3,3])),3)))), col=1, bty='n', cex=1.4)
 tmp1=range(pretty(Ychat[,1]))
 tmp2=range(pretty(Ychat[,3]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,3)], Sigma=Sig[c(1,3), c(1,3)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(1,3)], Sigma=SigN[c(1,3), c(1,3)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00005,0.00003,0.000009,0.000001,0.00000015))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00005,0.00003,0.000009,0.000001,0.00000015))

# ----------
 par(mar=c(0,0.25,0.25,0.25))
 oID = which(rowSums(cen[,c(2,3)]) == 0)
 plot(Ychat[oID,2], Ychat[oID,3], xlim=c(min(Ychat[,2],YNchat[,2]), max(Ychat[,2],YNchat[,2])), ylim=c(min(Ychat[,3],YNchat[,3]), max(Ychat[,3],YNchat[,3])), xaxt='n', yaxt='n', pch=16, cex=0.7)
 cID = which(rowSums(cen[,c(2,3)]) != 0)
 points(Ychat[cID, 2], Ychat[cID, 3], pch=15, col=3, cex=0.9, lwd=0.6) 
 points(YNchat[cID, 2], YNchat[cID, 3], pch=3, col="orange", cex=0.95, lwd=0.6)                                                                       
 abline(v=Yc[cen[,2]==1,2][1], h=Yc[cen[,3]==1,3][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[2,3]/(sqrt(Sig[2,2]*Sig[3,3])),3))),
                              bquote(hat(rho)[MVNC] == .(round(SigN[2,3]/(sqrt(SigN[2,2]*SigN[3,3])),3)))), col=1, bty='n', cex=1.4)
 axis(1, Xlab[[2]], cex.lab=0.9, lwd=0.45)
 tmp1=range(pretty(Ychat[,2]))
 tmp2=range(pretty(Ychat[,3]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(2,3)], Sigma=Sig[c(2,3), c(2,3)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(2,3)], Sigma=SigN[c(2,3), c(2,3)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00008,0.00003,0.00001,0.000003,0.0000004,0.000000015))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00008,0.00003,0.00001,0.000003,0.0000004,0.000000015))
  
# ---------- 
 par(mar=c(0.1,0,0.5,0))
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('Observations','Censored values'), fill=c(0,gray(.6)), bty='n', cex=1.3)
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('Observations','Censored values','MSLC recovered censored values','MVNC recovered censored values'), pch=c(16,17,15,3), col=c(1,"red",3,"orange"), bty='n', cex=1.2)
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('MSLC contours','MVNC contours'), lty=c(1,2), col=c("deeppink2", "blue"), lwd=c(0.9, 0.75), bty='n', cex=1.3)

 
 
 
 
 
# November -------------- 
rm(list = ls())
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
load(paste(MS.PATH, 'result/bioassay_nov.RData',sep=''))

# Setting -------------- 
brk = c(10, 15, 15)
cenrate = round(colMeans(cen)*100, 2)
skew = c(2.16, 2.06, 2.94)
kurt = c(7.58, 7.06, 11.04)
var.name = set
Xlab = as.list(p)
Xlab[[1]] = seq(-290,360.5,150)
Xlab[[2]] = seq(19.5,240.5,56.5)
Xlab[[3]] = seq(49.5,1800.5,550.5)
Ychat = estSC$yhat
YNchat = estNC$yhat

xli = c(-140,150,1000)
alim = c(0.0065,0.015,0.0025)
lim = c(0.011,0.03,0.005)
at = c(-0.001,-0.0069,-0.00073)
at1 = c(-0.0012,-0.0057,-0.00055)
bx = c(0.003, 0.008,0.001)
lex = c(0,0,0)
ley = c(0,0,0)

# Diagonal Entry -------------- 
i = 1 
par(mar=c(0,2.5,2.5,0.25))
plot(0:1,0:1, ylim=c(at[i]-0.001, lim[i]), xlim=c(min(Ychat[,i]), max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
h1 = hist(Yc[,i], breaks=brk[i], plot=F)
a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
a2 = length(h1$counts) + 1
hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)), col=c(rep(gray(.6),a2-a1), rep(0,a1)), xlab='', ylab='Density', add=T)
boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
text(xli[i], alim[i], var.name[i], cex=2.9)
legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                    paste('skew = ', skew[i], sep=''),
                    paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1.5)
i = 2
par(mar=c(0,0.25,0.25,0.25))
plot(0:1,0:1, ylim=c(at[i]-0.001, lim[i]), xlim=c(min(Ychat[,i]), max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
h1 = hist(Yc[,i], breaks=brk[i], plot=F)
a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
a2 = length(h1$counts) + 1
hist(Yc[,i], breaks=brk[i], prob=T, xaxt='n', ylim=c(0, max(h1$density)), xlim=c(min(Ychat[,i]), max(Ychat[,i])),col=c(rep(0,a2)), xlab='', ylab='Density', add=T) 
boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
text(xli[i], alim[i], var.name[i], cex=2.9)
legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                    paste('skew = ', skew[i], sep=''),
                    paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1.5)
i = 3
par(mar=c(0,0.25,0.25,2.5))
plot(0:1,0:1, ylim=c(at[i], lim[i]), xlim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))), type='n', xlab='', ylab='', axes=T, main='', yaxt='n', xaxt='n')
h1 = hist(Yc[,i], breaks=brk[i], plot=F)
a1 = sum(h1$breaks>Yc[cen[,i]==1,i][1])
a2 = length(h1$counts) + 1
hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)),col=c(rep(0,a2)), xlab='', ylab='Density', add=T)
boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
text(xli[i], alim[i], var.name[i], cex=2.9)
legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                    paste('skew = ', skew[i], sep=''),
                    paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1.5)

# Upper trianglar -------------- 
for(i in 1)for(j in 2){
 par(mar=c(0,0.25,2.5,0.25))
 oID = which(rowSums(cen[,c(j,i)]) == 0)
 plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
 cID = which(rowSums(cen[,c(j,i)]) != 0)
 points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
 abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
 text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
 text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
 legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
 axis(3, Xlab[[j]], cex.lab=0.9, lwd=0.45)
}
for(i in 1)for(j in 3){
 par(mar=c(0,0.25,2.5,2.5))
 oID = which(rowSums(cen[,c(j,i)]) == 0)
 plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
 cID = which(rowSums(cen[,c(j,i)]) != 0)
 points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
 abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
 text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
 text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
 legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
 axis(3, Xlab[[j]], cex.lab=0.9, lwd=0.45)
 axis(4, Xlab[[i]], cex.lab=0.9, lwd=0.45)
}
for(i in 2)for(j in 3){
 par(mar=c(0,0.25,0.25,2.5))
 oID = which(rowSums(cen[,c(j,i)]) == 0)
 plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.7, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
 cID = which(rowSums(cen[,c(j,i)]) != 0)
 points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
 abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
 text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1),cex=1.4)
 text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5),cex=1.4)
 legend('topleft', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1.4)
 axis(4, Xlab[[i]], cex.lab=0.9, lwd=0.45)
}
# Lower trianglar --------------
m=100
mu = estSC$para$mu
Sig = estSC$para$Sigma
nu = estSC$para$nu
muN = estNC$para$mu
SigN = estNC$para$Sigma

par(mar=c(0.25,2.5,0.25,0.25))
for(i in 2){
 oID = which(rowSums(cen[,c(1,i)]) == 0)
 plot(Ychat[oID,1], Ychat[oID,i], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.7)
 axis(2, Xlab[[i]], cex.lab=0.9, lwd=0.45)
 cID = which(rowSums(cen[,c(1,i)]) != 0)
 points(Ychat[cID, 1], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, 1], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
 abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,i]/(sqrt(Sig[1,1]*Sig[i,i])),3))),
                              bquote(hat(rho)[MVNC] == .(round(SigN[1,i]/(sqrt(SigN[1,1]*SigN[i,i])),3)))), col=1, bty='n', cex=1.4)
 tmp1=range(pretty(Ychat[,1]))
 tmp2=range(pretty(Ychat[,i]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,i)], Sigma=Sig[c(1,i), c(1,i)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(1,i)], Sigma=SigN[c(1,i), c(1,i)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00002,0.000005,0.00000055,0.00000008,0.00000002,0.000000008))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00002,0.000005,0.00000055,0.00000008,0.00000002,0.000000008))
} 
par(mar=c(0,2.5,0.25,0.25))
oID = which(rowSums(cen[,c(1,3)]) == 0)
plot(Ychat[oID,1], Ychat[oID,3], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,3],YNchat[,3]), max(Ychat[,3],YNchat[,3])), xaxt='n', yaxt='n', pch=16, cex=0.7)
axis(1, Xlab[[1]], cex.lab=0.9, lwd=0.45)
axis(2, Xlab[[3]], cex.lab=0.9, lwd=0.45)
cID = which(rowSums(cen[,c(1,3)]) != 0)
points(Ychat[cID, 1], Ychat[cID, 3], pch=15, col=3, cex=0.9, lwd=0.6)
points(YNchat[cID, 1], YNchat[cID, 3], pch=3, col="orange", cex=0.95, lwd=0.6)
abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,3]==1,3][1], lty=4, col="lavenderblush4", lwd=0.8) 
legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,3]/(sqrt(Sig[1,1]*Sig[3,3])),3))),
                            bquote(hat(rho)[MVNC] == .(round(SigN[1,3]/(sqrt(SigN[1,1]*SigN[3,3])),3)))), col=1, bty='n', cex=1.4)
tmp1=range(pretty(Ychat[,1]))
tmp2=range(pretty(Ychat[,3]))
xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
den=numeric(m^2)
for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,3)], Sigma=Sig[c(1,3), c(1,3)], nu=nu, distr='MSL')
den = matrix(den, m, m)
denN = dmsl(xx, mu=muN[c(1,3)], Sigma=SigN[c(1,3), c(1,3)], distr='MVN')
denN = matrix(denN, m, m)
contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.000002,0.0000004,0.00000004,0.00000001,0.000000003,0.0000000007,0.00000000009))
contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.000002,0.0000004,0.00000004,0.00000001,0.000000003,0.0000000007,0.00000000009))

# ----------
par(mar=c(0,0.25,0.25,0.25))
oID = which(rowSums(cen[,c(2,3)]) == 0)
plot(Ychat[oID,2], Ychat[oID,3], xlim=c(min(Ychat[,2],YNchat[,2]), max(Ychat[,2],YNchat[,2])), ylim=c(min(Ychat[,3],YNchat[,3]), max(Ychat[,3],YNchat[,3])), xaxt='n', yaxt='n', pch=16, cex=0.7)
cID = which(rowSums(cen[,c(2,3)]) != 0)
points(Ychat[cID, 2], Ychat[cID, 3], pch=15, col=3, cex=0.9, lwd=0.6)                                                                        
points(YNchat[cID, 2], YNchat[cID, 3], pch=3, col="orange", cex=0.95, lwd=0.6)                                                                       
abline(v=Yc[cen[,2]==1,2][1], h=Yc[cen[,3]==1,3][1], lty=4, col="lavenderblush4", lwd=0.8) 
legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[2,3]/(sqrt(Sig[2,2]*Sig[3,3])),3))),
                            bquote(hat(rho)[MVNC] == .(round(SigN[2,3]/(sqrt(SigN[2,2]*SigN[3,3])),3)))), col=1, bty='n', cex=1.4)
axis(1, Xlab[[2]], cex.lab=0.9, lwd=0.45)
tmp1=range(pretty(Ychat[,2]))
tmp2=range(pretty(Ychat[,3]))
xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
den=numeric(m^2)
for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(2,3)], Sigma=Sig[c(2,3), c(2,3)], nu=nu, distr='MSL')
den = matrix(den, m, m)
denN = dmsl(xx, mu=muN[c(2,3)], Sigma=SigN[c(2,3), c(2,3)], distr='MVN')
denN = matrix(denN, m, m)
contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.000027,0.000015,0.000003,0.0000005,0.0000001,0.00000003,0.00000001,0.000000001,0.0000000003))
contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.000027,0.000015,0.000003,0.0000005,0.0000001,0.00000003,0.00000001,0.000000001,0.0000000003))

# ---------- 
par(mar=c(0.1,0,0.5,0))
plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
legend("center", c('Observations','Censored values'), fill=c(0,gray(.6)), bty='n', cex=1.3)
plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
legend("center", c('Observations','Censored values','MSLC recovered censored values','MVNC recovered censored values'), pch=c(16,17,15,3), col=c(1,"red",3,"orange"), bty='n', cex=1.2)
plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
legend("center", c('MSLC contours','MVNC contours'), lty=c(1,2), col=c("deeppink2", "blue"), lwd=c(0.9, 0.75), bty='n', cex=1.3)

dev.off()
 
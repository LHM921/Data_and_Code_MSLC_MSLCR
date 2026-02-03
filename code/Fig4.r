################################################################################ 
# Purpose: produce Figure 4
#         Visualization of composition data along with the analysis results
# Input: load 'result/lpdata.RData' file
# Output: produce 'result/fig4.pdf'
################################################################################

library(mvtnorm)
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig4.pdf"), width=8.5, height=8.5)

load(paste(MS.PATH, 'result/lpdata.RData',sep=''))

nf = layout(matrix(c(1,16,17,18,19,0,
                     6,2,20,21,22,26,
                     7,10,3,23,24,27,
                     8,11,13,4,25,28,
                     9,12,14,15,5,0), 6, 5), c(5, rep(4,3),4.5), c(5,rep(4,4),2.9))

# Setting -------------- 
brk = c(15, 15, 20, 15, 18)
cenrate = round(colMeans(cen)*100, 2)
skew = c(1.37, 3.81, 1.44, 1.86, 5.39)
kurt = c(4.53, 26.19, 4.95, 7.49, 34.53)
var.name = set
Xlab = as.list(p)
Xlab[[1]] = seq(-55,90,35)
Xlab[[2]] = seq(-11,105,41)
Xlab[[3]] = seq(-14,60,14)
Xlab[[4]] = seq(-31,920,375)
Xlab[[5]] = seq(10,1120,500)
Ychat = estSC$yhat
YNchat = estNC$yhat

xli = c(-22,70,46,555,550)
alim = c(0.035,0.048,0.03,0.004,0.0073)
lim = c(0.05,0.06,0.07,0.005,0.01)
at = c(-0.011,-0.008,-0.01,0.003,0.001)
at1 = c(-0.01,-0.008,-0.01,-0.0007,-0.0016)
bx = c(0.013, 0.014, 0.016, 0.0013, 0.0035)
lex = c(0, 0.02,0,0.015,0)
ley = c(0, 0.28,0,0.3,0)

# Diagonal Entry -------------- 
i = 1
 par(mar=c(0,2.5,2.5,0.25))
 plot(0:1,0:1, ylim=c(at[i]-0.004, lim[i]), xlim=c(min(Ychat[,i]), max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
 h1 = hist(Yc[,i], breaks=brk[i], plot=F)
 a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
 a2 = length(h1$counts) + 1
 hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)), col=c(rep(gray(.6),a2-a1), rep(0,a1)), xlab='', ylab='Density', add=T)
 boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
 text(xli[i], alim[i], var.name[i], cex=2.5)
 legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                     paste('skew = ', skew[i], sep=''),
                     paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1)
 for(i in c(2,3,4)){
  par(mar=c(0,0.25,0.25,0.25))
  plot(0:1,0:1, ylim=c(at[i]-0.004, lim[i]), xlim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))), type='n', xlab='', ylab='', axes=T, main='', yaxt='n', xaxt='n')
  h1 = hist(Yc[,i], breaks=brk[i], plot=F)
  a1 = sum(h1$breaks>Yc[cen[,i]==1,i][1])
  a2 = length(h1$counts) + 1
  hist(Yc[,i], breaks=brk[i], prob=T, ylim=c(0, max(h1$density)), col=c(rep(gray(.6),a2-a1), rep(0,a1)), xlab='', ylab='Density', add=T)
  boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
  text(xli[i], alim[i], var.name[i], cex=2.5)
  legend('topright', c(paste('r = ', cenrate[i], '%', sep=''),
                       paste('skew = ', skew[i], sep=''),
                       paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1)
}  
i = 5 
 par(mar=c(0,0.25,0.25,2.5))
 plot(0:1,0:1, ylim=c(at[i]-0.004, lim[i]), xlim=c(min(Ychat[,i])*1.05, max(Ychat[,i])), type='n', xlab='', ylab='', xaxt='n', yaxt='n',  main='', las=1)
 h1 = hist(Yc[,i], breaks=brk[i], plot=F)
 a1 = sum(h1$breaks > Yc[cen[,i]==1,i][1])
 a2 = length(h1$counts) + 1
 hist(Yc[,i], breaks=brk[i], prob=T, xaxt='n', ylim=c(0, max(h1$density)), xlim=c(min(Ychat[,i]), max(Ychat[,i])), col=c(rep(0,a2)), xlab='', ylab='Density', add=T) 
 boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at1[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
 text(xli[i], alim[i], var.name[i], cex=2.5)
 legend('center', c(paste('r = ', cenrate[i], '%', sep=''),
                     paste('skew = ', skew[i], sep=''),
                     paste('kurt = ', kurt[i], sep='')), inset = c(lex[i], ley[i]), col=1, bty='n', cex=1)
 
# Upper trianglar -------------- 
 for(i in 1)for(j in (i+1):(p-1)){
   par(mar=c(0,0.25,2.5,0.25))
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.5, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1))
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5))
   legend('topright', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1)
   if (i==1 & j %in% c(2,3,4)){axis(3, Xlab[[j]], cex.lab=0.5, lwd=0.45)}
 }
 for(i in 1)for(j in 5){
   par(mar=c(0,0.25,2.5,2.5))
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.5, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1))
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5))
   legend('topright', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1)
   axis(3, Xlab[[j]], cex.lab=0.5, lwd=0.45)
   axis(4, Xlab[[i]], cex.lab=0.5, lwd=0.45)
 }
 for(i in (1+1): (p-1))for(j in (i+1):p){
   if(j %in% c(5)){par(mar=c(0,0.25,0.25,2.5))}else{par(mar=c(0,0.25,0.25,0.25))}
   oID = which(rowSums(cen[,c(j,i)]) == 0)
   plot(Yc[oID,j], Yc[oID,i], col=1, pch=16, cex=0.5, xlim=c(min(c(Yc[,j],Ychat[,j])), max(c(Yc[,j],Ychat[,j]))), ylim=c(min(c(Yc[,i],Ychat[,i])), max(c(Yc[,i],Ychat[,i]))),xaxt="n",yaxt="n")
   cID = which(rowSums(cen[,c(j,i)]) != 0)
   points(Yc[cID, j], Yc[cID, i], pch=17, col="red", cex=1, lwd=0.7)
   abline(h=Yc[cen[,i]==1,i][1], v=Yc[cen[,j]==1,j][1],lty=4, col="lavenderblush4", lwd=0.8)
   text(x=par("usr")[2], y=Yc[cen[,i]==1,i][1], labels=paste0("DL=", Yc[cen[,i]==1,i][1]),adj=c(1,1))
   text(x=Yc[cen[,j]==1,j][1], y=par("usr")[3], labels=paste0("DL=", Yc[cen[,j]==1,j][1]),adj=c(-0.08,-0.5))
   legend('topright', paste('Corr = ', round(cor(Yc[,j],Yc[,i]),3), sep=''), col=4, bty='n', cex=1)
   if (i %in% c(1,2,3,4) & j==5){axis(4, Xlab[[i]], cex.lab=0.5, lwd=0.45)}
 }
 
# Lower trianglar --------------
m=96 
mu = estSC$para$mu
Sig = estSC$para$Sigma
nu = estSC$para$nu
muN = estNC$para$mu
SigN = estNC$para$Sigma

 par(mar=c(0.25,2.5,0.25,0.25))
 for(i in 2:3){
  oID = which(rowSums(cen[,c(1,i)]) == 0)
  plot(Ychat[oID,1], Ychat[oID,i], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.5)
  axis(2, Xlab[[i]], cex.lab=0.5, lwd=0.45)
  cID = which(rowSums(cen[,c(1,i)]) != 0)
  points(Ychat[cID, 1], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)
  points(YNchat[cID, 1], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
  abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
  legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,i]/(sqrt(Sig[1,1]*Sig[i,i])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[1,i]/(sqrt(SigN[1,1]*SigN[i,i])),3)))), col=1, bty='n', cex=0.8)
  tmp1=range(pretty(Ychat[,1]))
  tmp2=range(pretty(Ychat[,i]))
  xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
  den=numeric(m^2)
  for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,i)], Sigma=Sig[c(1,i), c(1,i)], nu=nu, distr='MSL')
  den = matrix(den, m, m)
  denN = dmsl(xx, mu=muN[c(1,i)], Sigma=SigN[c(1,i), c(1,i)], distr='MVN')
  denN = matrix(denN, m, m)
  contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.001,0.0005,0.00005,0.00001,0.000001))
  contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.001,0.0005,0.00005,0.00001,0.000001))
 } 
 for(i in 4){
   oID = which(rowSums(cen[,c(1,i)]) == 0)
   plot(Ychat[oID,1], Ychat[oID,i], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.5)
   axis(2, Xlab[[i]], cex.lab=0.5, lwd=0.45)
   cID = which(rowSums(cen[,c(1,i)]) != 0)
   points(Ychat[cID, 1], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)
   points(YNchat[cID, 1], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
   abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
   legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,i]/(sqrt(Sig[1,1]*Sig[i,i])),3))),
                                 bquote(hat(rho)[MVNC] == .(round(SigN[1,i]/(sqrt(SigN[1,1]*SigN[i,i])),3)))), col=1, bty='n', cex=0.8)
   tmp1=range(pretty(Ychat[,1]))
   tmp2=range(pretty(Ychat[,i]))
   xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
   den=numeric(m^2)
   for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,i)], Sigma=Sig[c(1,i), c(1,i)], nu=nu, distr='MSL')
   den = matrix(den, m, m)
   denN = dmsl(xx, mu=muN[c(1,i)], Sigma=SigN[c(1,i), c(1,i)], distr='MVN')
   denN = matrix(denN, m, m)
   contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00005,0.00001,0.000001,0.00000007,0.000000006))
   contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00005,0.00001,0.000001,0.00000007,0.000000006))
 } 
 par(mar=c(0,2.5,0.25,0.25))
 oID = which(rowSums(cen[,c(1,5)]) == 0)
 plot(Ychat[oID,1], Ychat[oID,5], xlim=c(min(Ychat[,1],YNchat[,1]), max(Ychat[,1],YNchat[,1])), ylim=c(min(Ychat[,5],YNchat[,5]), max(Ychat[,5],YNchat[,5])), xaxt='n', yaxt='n', pch=16, cex=0.5)
 axis(1, Xlab[[1]], cex.lab=0.5, lwd=0.45)
 axis(2, Xlab[[5]], cex.lab=0.5, lwd=0.45)
 cID = which(rowSums(cen[,c(1,5)]) != 0)
 points(Ychat[cID, 1], Ychat[cID, 5], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, 1], YNchat[cID, 5], pch=3, col="orange", cex=0.95, lwd=0.6)
 abline(v=Yc[cen[,1]==1,1][1], h=Yc[cen[,5]==1,5][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[1,5]/(sqrt(Sig[1,1]*Sig[5,5])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[1,5]/(sqrt(SigN[1,1]*SigN[5,5])),3)))), col=1, bty='n', cex=0.8)
 tmp1=range(pretty(Ychat[,1]))
 tmp2=range(pretty(Ychat[,5]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(1,5)], Sigma=Sig[c(1,5), c(1,5)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(1,5)], Sigma=SigN[c(1,5), c(1,5)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.0003,0.00005,0.00001,0.000001,0.0000001))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.0003,0.00005,0.00001,0.000001,0.0000001))

# ----------
 par(mar=c(0.25,0.25,0.25,0.25))
 for(i in 3){
  oID = which(rowSums(cen[,c(2,i)]) == 0)
  plot(Ychat[oID,2], Ychat[oID,i], xlim=c(min(Ychat[,2],YNchat[,2]), max(Ychat[,2],YNchat[,2])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.5)
  cID = which(rowSums(cen[,c(2,i)]) != 0)
  points(Ychat[cID, 2], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)                                                                        
  points(YNchat[cID, 2], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
  abline(v=Yc[cen[,2]==1,2][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
  legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[2,i]/(sqrt(Sig[2,2]*Sig[i,i])),3))),
                                bquote(hat(rho)[MVNC] == .(round(SigN[2,i]/(sqrt(SigN[2,2]*SigN[i,i])),3)))), col=1, bty='n', cex=0.8)
  tmp1=range(pretty(Ychat[,2]))
  tmp2=range(pretty(Ychat[,i]))
  xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
  den=numeric(m^2)
  for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(2,i)], Sigma=Sig[c(2,i), c(2,i)], nu=nu, distr='MSL')
  den = matrix(den, m, m)
  denN = dmsl(xx, mu=muN[c(2,i)], Sigma=SigN[c(2,i), c(2,i)], distr='MVN')
  denN = matrix(denN, m, m)
  contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.0007,0.0001,0.00001,0.000001))
  contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.0007,0.0001,0.00001,0.000001))
 } 
 par(mar=c(0.25,0.25,0.25,0.25))
 for(i in 4){
   oID = which(rowSums(cen[,c(2,i)]) == 0)
   plot(Ychat[oID,2], Ychat[oID,i], xlim=c(min(Ychat[,2],YNchat[,2]), max(Ychat[,2],YNchat[,2])), ylim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), xaxt='n', yaxt='n', pch=16, cex=0.5)
   cID = which(rowSums(cen[,c(2,i)]) != 0)
   points(Ychat[cID, 2], Ychat[cID, i], pch=15, col=3, cex=0.9, lwd=0.6)                                                                        
   points(YNchat[cID, 2], YNchat[cID, i], pch=3, col="orange", cex=0.95, lwd=0.6)
   abline(v=Yc[cen[,2]==1,2][1], h=Yc[cen[,i]==1,i][1], lty=4, col="lavenderblush4", lwd=0.8) 
   legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[2,i]/(sqrt(Sig[2,2]*Sig[i,i])),3))),
                                 bquote(hat(rho)[MVNC] == .(round(SigN[2,i]/(sqrt(SigN[2,2]*SigN[i,i])),3)))), col=1, bty='n', cex=0.8)
   tmp1=range(pretty(Ychat[,2]))
   tmp2=range(pretty(Ychat[,i]))
   xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
   den=numeric(m^2)
   for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(2,i)], Sigma=Sig[c(2,i), c(2,i)], nu=nu, distr='MSL')
   den = matrix(den, m, m)
   denN = dmsl(xx, mu=muN[c(2,i)], Sigma=SigN[c(2,i), c(2,i)], distr='MVN')
   denN = matrix(denN, m, m)
   contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.0005,0.0001,0.00001,0.000001,0.0000001,0.00000002))
   contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.0005,0.0001,0.00001,0.000001,0.0000001,0.00000002))
 } 
 par(mar=c(0,0.25,0.25,0.25))
 oID = which(rowSums(cen[,c(2,5)]) == 0)
 plot(Ychat[oID,2], Ychat[oID,5], xlim=c(min(Ychat[,2],YNchat[,2]), max(Ychat[,2],YNchat[,2])), ylim=c(min(Ychat[,5],YNchat[,5]), max(Ychat[,5],YNchat[,5])), xaxt='n', yaxt='n', pch=16, cex=0.5)
 cID = which(rowSums(cen[,c(2,5)]) != 0)
 points(Ychat[cID, 2], Ychat[cID, 5], pch=15, col=3, cex=0.9, lwd=0.6)                                                                        
 points(YNchat[cID, 2], YNchat[cID, 5], pch=3, col="orange", cex=0.95, lwd=0.6)
 abline(v=Yc[cen[,2]==1,2][1], h=Yc[cen[,5]==1,5][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[2,5]/(sqrt(Sig[2,2]*Sig[5,5])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[2,5]/(sqrt(SigN[2,2]*SigN[5,5])),3)))), col=1, bty='n', cex=0.8)
 axis(1, Xlab[[2]], cex.lab=0.5, lwd=0.45)
 tmp1=range(pretty(Ychat[,2]))
 tmp2=range(pretty(Ychat[,5]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(2,5)], Sigma=Sig[c(2,5), c(2,5)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(2,5)], Sigma=SigN[c(2,5), c(2,5)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.0001,0.000007,0.000001,0.0000001))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.0001,0.000007,0.000001,0.0000001))
  
# ---------- 
 par(mar=c(0.25,0.25,0.25,0.25))
 oID = which(rowSums(cen[,c(3,4)]) == 0)
 plot(Ychat[oID,3], Ychat[oID,4], xlim=c(min(Ychat[,3],YNchat[,3]), max(Ychat[,3],YNchat[,3])), ylim=c(min(Ychat[,4],YNchat[,4]), max(Ychat[,4],YNchat[,4])), xaxt='n', yaxt='n', pch=16, cex=0.5)
 cID = which(rowSums(cen[,c(3,4)]) != 0)
 points(Ychat[cID, 3], Ychat[cID, 4], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, 3], YNchat[cID, 4], pch=3, col="orange", cex=0.95, lwd=0.6)
 abline(v=Yc[cen[,3]==1,3][1], h=Yc[cen[,4]==1,4][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[3,4]/(sqrt(Sig[3,3]*Sig[4,4])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[3,4]/(sqrt(SigN[3,3]*SigN[4,4])),3)))), col=1, bty='n', cex=0.8)
 tmp1=range(pretty(Ychat[,3]))
 tmp2=range(pretty(Ychat[,4]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(3,4)], Sigma=Sig[c(3,4), c(3,4)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(3,4)], Sigma=SigN[c(3,4), c(3,4)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.0001,0.00001,0.0000001,0.00000001))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.0001,0.00001,0.0000001,0.00000001))

 par(mar=c(0,0.25,0.25,0.25))
 for(i in 3){
 oID = which(rowSums(cen[,c(i,5)]) == 0)
 plot(Ychat[oID,i], Ychat[oID,5], xlim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), ylim=c(min(Ychat[,5],YNchat[,5]), max(Ychat[,5],YNchat[,5])), xaxt='n', yaxt='n', pch=16, cex=0.5)
 cID = which(rowSums(cen[,c(i,5)]) != 0)
 points(Ychat[cID, i], Ychat[cID, 5], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, i], YNchat[cID, 5], pch=3, col="orange", cex=0.95, lwd=0.6)
 axis(1, Xlab[[i]], cex.lab=0.5, lwd=0.45) 
 abline(v=Yc[cen[,i]==1,i][1], h=Yc[cen[,5]==1,5][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[i,5]/(sqrt(Sig[i,i]*Sig[5,5])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[i,5]/(sqrt(SigN[i,i]*SigN[5,5])),3)))), col=1, bty='n', cex=0.8)
 tmp1=range(pretty(Ychat[,i]))
 tmp2=range(pretty(Ychat[,5]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(i,5)], Sigma=Sig[c(i,5), c(i,5)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(i,5)], Sigma=SigN[c(i,5), c(i,5)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00007,0.00001,0.000001,0.0000001))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00007,0.00001,0.000001,0.0000001))
 }
 par(mar=c(0,0.25,0.25,0.25))
 for(i in 4){
 oID = which(rowSums(cen[,c(i,5)]) == 0)
 plot(Ychat[oID,i], Ychat[oID,5], xlim=c(min(Ychat[,i],YNchat[,i]), max(Ychat[,i],YNchat[,i])), ylim=c(min(Ychat[,5],YNchat[,5]), max(Ychat[,5],YNchat[,5])), xaxt='n', yaxt='n', pch=16, cex=0.5)
 cID = which(rowSums(cen[,c(i,5)]) != 0)
 points(Ychat[cID, i], Ychat[cID, 5], pch=15, col=3, cex=0.9, lwd=0.6)
 points(YNchat[cID, i], YNchat[cID, 5], pch=3, col="orange", cex=0.95, lwd=0.6)
 axis(1, Xlab[[i]], cex.lab=0.5, lwd=0.45) 
 abline(v=Yc[cen[,i]==1,i][1], h=Yc[cen[,5]==1,5][1], lty=4, col="lavenderblush4", lwd=0.8) 
 legend('topleft', legend = c(bquote(hat(rho)[MSLC] == .(round(Sig[i,5]/(sqrt(Sig[i,i]*Sig[5,5])),3))),
                               bquote(hat(rho)[MVNC] == .(round(SigN[i,5]/(sqrt(SigN[i,i]*SigN[5,5])),3)))), col=1, bty='n', cex=0.8)
 tmp1=range(pretty(Ychat[,i]))
 tmp2=range(pretty(Ychat[,5]))
 xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
 den=numeric(m^2)
 for(k in 1: m^2) den[k] = dmsl(as.numeric(c(xx[k, ])), mu=mu[c(i,5)], Sigma=Sig[c(i,5), c(i,5)], nu=nu, distr='MSL')
 den = matrix(den, m, m)
 denN = dmsl(xx, mu=muN[c(i,5)], Sigma=SigN[c(i,5), c(i,5)], distr='MVN')
 denN = matrix(denN, m, m)
 contour(x1, x2, den, labcex = 0.3, add=T, lty=1, col="deeppink2", lwd=0.4, drawlabels = T, levels=c(0.00001,0.000001,0.0000001,0.00000001,0.000000001))
 contour(x1, x2, denN, labcex = 0.3, add=T, lty=2, col="blue", lwd=0.5, drawlabels = T, levels=c(0.00001,0.000001,0.0000001,0.00000001,0.000000001))
}
 par(mar=c(0.1,0,2,0))
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('Observations','Censored values'), fill=c(0,gray(.6)), bty='n', cex=0.9)
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('Observations','Censored values','MSLC recovered censored values','MVNC recovered censored values'), pch=c(16,17,15,3), col=c(1,"red",3,"orange"), bty='n', cex=0.8)
 plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
 legend("center", c('MSLC contours','MVNC contours'), lty=c(1,2), col=c("deeppink2", "blue"), lwd=c(0.9, 0.75), bty='n', cex=0.9)
 
dev.off()

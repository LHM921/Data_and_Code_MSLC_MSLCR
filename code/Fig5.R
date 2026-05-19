################################################################################ 
# Purpose: produce Figure 4
#     3D scatter plots and fitted regression planes for lake data
# Package: scatterplot3d
# Input: load 'result/lake.RData' file
# Output: produce 'result/fig5.pdf'
################################################################################

#install.packages("scatterplot3d")
library(scatterplot3d)

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig5.pdf"), width=13.8, height=10)
load(paste(MS.PATH, 'result/lake.RData',sep=''))

nf = layout(matrix(c(0,0,0,
                     1,2,3,
                     4,5,6,
                     0,0,0), 4, 3, byrow = T), c(5,5,5),c(0.5,6,6,0.5))

# Setting -------------- 
cenrate = round(colMeans(cen)*100, 2)
skew = c(1.75, 1.27, 0.62)
kurt = c(9.08, 6.70, 2.87)

Ychat = estSCR$yoyc
YNchat = estNCR$yoyc
# YTchat = estTCR$yoyc
beta = estSCR$para$beta
Sig = estSCR$para$Sigma
nu = estSCR$para$nu
betaN = estNCR$para$beta
SigN = estNCR$para$Sigma
# betaT = estTCR$para$beta
# SigT = estTCR$para$Sigma

results_EST = as.data.frame(rbind(c('MSLR', estSR$EST), c('MVNR', estNR$EST, 0), c('MVTR', estTR$EST), 
                                  c('MSLCR', estSCR$EST), c('MVNCR', estNCR$EST, 0), c('MVTCR', estTCR$EST)))
colnames(results_EST) <- c("Model", paste0("EST", 1:length(estSR$EST)))
results_SE = as.data.frame(rbind(c('MSLR', estSR$IM$SD), c('MVNR', estNR$IM$SD),c('MVTR', estTR$IM$SD),
                                 c('MSLCR', estSCR$IM$SD), c('MVNCR', estNCR$IM$SD), c('MVTCR', estTCR$IM$SD)))
colnames(results_SE) <- c("Model", paste0("SE", 1:length(estSR$IM$SD)))

# Wald -------------- 
X <-  as.matrix(scale(df[, setx]))
colnames(X) <- c("Lat","Long","DepthSite","DepthSmp")
setx <- c("Lat","Long","DepthSite","DepthSmp")
wald_score <- as.data.frame(
  lapply(1:(length(setx)*length(sety))+1, function(i) abs(as.numeric(results_EST[[i]])) / as.numeric(results_SE[[i]]))
)
colnames(wald_score) <- 1:(length(setx)*length(sety))
wald_score_top2 <- function(idx) {
  cols <- colMeans(wald_score)[idx]
  names(sort(cols, decreasing = TRUE))[1:2]
}
y1x <- wald_score_top2(1:length(setx)); y2x <- wald_score_top2((length(setx)+1):(length(setx)*2)); y3x <- wald_score_top2((length(setx)*2+1):(length(setx)*3)) 
selectedx <- X[,rep(c(setx),3)[as.integer(c(y1x, y2x, y3x))]]


# Ori Yc --------------
i=1
oID <- which(cen[,i] == 0) ; cID <- which(cen[,i] != 0)        
col_vec <- rep("black", nrow(X)) ; col_vec[cID] <- "red"
pch_vec <- rep(16, nrow(X)) ; pch_vec[cID] <- 17
s3d <- scatterplot3d(Yc[,1], Yc[,2], Yc[,3], cex.lab=1.5, cex.axis = 1.3,
                     xlab=sety[1], ylab=sety[2], zlab=sety[3],
                     pch=pch_vec, color = col_vec, grid=F, cex=1.3)

# Text --------------
text_info <- list(
  TDP = c(paste(sety[1], ':', sep=''), paste('r = ', cenrate[1], '%', sep=''), paste('skew = ', skew[1], sep=''),paste('kurt = ', kurt[1], sep='')),
  PC = c(paste(sety[2], ':', sep=''), paste('r = ', cenrate[2], '%', sep=''), paste('skew = ', skew[2], sep=''),paste('kurt = ', kurt[2], sep='')),
  PN = c(paste(sety[3], ':', sep=''), paste('r = ', cenrate[3], '%', sep=''), paste('skew = ', skew[3], sep=''),paste('kurt = ', kurt[3], sep=''))
)
plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n',xlim=c(0,3), ylim=c(0,0.1))
x_pos <- c(0.3, 1.3, 2.3)  
y_pos <- c(0.095,0.085,0.075,0.065)  
for(i in 1:length(text_info)){
  text(x=x_pos[i], y=y_pos, labels=text_info[[i]], adj=c(0,0.4), cex=1.5)
}

# MSLCR sigma_hat
text(0.86, 0.042, bquote(hat(Sigma)[MSLCR] == "["), cex=3)
for(i in 1:nrow(Sig)){
  for(j in i:ncol(Sig)){ 
    text(1.27+0.46*j, 0.05-0.0053*i, sprintf("%.3f", Sig[i, j]), cex=1.5)
  }
}
text(2.89, 0.04, expression("]"), cex=3)

# MVNCR sigma_hat
text(0.86, 0.022, bquote(hat(Sigma)[MVNCR] == "["), cex=3)
for(i in 1:nrow(SigN)){
  for(j in i:ncol(SigN)){ 
    text(1.27+0.46*j, 0.03-0.0053*i, sprintf("%.3f", SigN[i, j]), cex=1.5)
  }
}
text(2.89, 0.02, expression("]"), cex=3)

# MVTCR sigma_hat
# text(0.9, 0.02, bquote(hat(Sigma)[MVTCR] == "["), cex=3)
# for(i in 1:nrow(SigT)){
#   for(j in i:ncol(SigT)){ 
#     text(1.14+0.3*j, 0.03-0.0053*i, sprintf("%.3f", SigT[i, j]), cex=1.5)
#   }
# }
# text(2.22, 0.02, expression("]"), cex=3)


plot(0:1,0:1, type='n', xlab='', ylab='', axes=F,  main='', bty='n')
legend("center", c('Observations','Censored values','MSLCR recovered censored values','MVNCR recovered censored values'),
       pch=c(16,17,15,3), col=c(1,"red",3,"orange"), bty='n', cex=1.6)
legend(x=0.025,y=0.4, c('MSLCR fitted plane','MVNCR fitted plane'),
       fill=c(rgb(1, 0.7, 0.8, 0.58), rgb(0.2, 0.6, 0.8, 0.4)), border = "white",bty='n', cex=1.6, x.intersp=0.67)

# Model fitting --------------
i=1
oID <- which(cen[,i] == 0) ; cID <- which(cen[,i] != 0)        
col_vec <- rep("black", nrow(X)) ; col_vec[cID] <- 3
pch_vec <- rep(16, nrow(X)) ; pch_vec[cID] <- 15
s3d <- scatterplot3d(selectedx[,(2*i-1)], selectedx[,(2*i)], Ychat[,i], cex.lab=1.5, cex.axis = 1.3,
              xlab=colnames(selectedx)[2*i-1], ylab=colnames(selectedx)[2*i], zlab=sety[i],
              pch=pch_vec, color = col_vec, grid=F, cex=1.3, angle = 185)
s3d$points3d(selectedx[cID,(2*i-1)], selectedx[cID,(2*i)], YNchat[cID,i],
             pch=3, col = "orange", cex=1.3)
s3d$plane3d(Intercept = 0, betaN[(2*i-1),], betaN[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(0.2, 0.6, 0.8, 0.2), border = NA))
# s3d$plane3d(Intercept = 0, betaT[(2*i-1),], betaT[(2*i),], draw_polygon = T, draw_lines = F,
#             polygon_args = list(col = rgb(1, 1, 0.6, 0.4), border = NA))
s3d$plane3d(Intercept = 0, beta[(2*i-1),], beta[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(1, 0.7, 0.8, 0.4), border = NA))

i=2
s3d <- scatterplot3d(selectedx[,(2*i-1)], selectedx[,(2*i)], Ychat[,i], cex.lab=1.5, cex.axis = 1.3,
                     xlab=colnames(selectedx)[2*i-1], ylab=colnames(selectedx)[2*i], zlab=sety[i],
                     pch=16, color = "black", grid=F, cex=1.3, angle = 225)
s3d$plane3d(Intercept = 0, betaN[(2*i-1),], betaN[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(0.2, 0.6, 0.8, 0.2), border = NA))
# s3d$plane3d(Intercept = 0, betaT[(2*i-1),], betaT[(2*i),], draw_polygon = T, draw_lines = F,
#             polygon_args = list(col = rgb(1, 1, 0.6, 0.4), border = NA))
s3d$plane3d(Intercept = 0, beta[(2*i-1),], beta[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(1, 0.7, 0.8, 0.4), border = NA))

i=3
s3d <- scatterplot3d(selectedx[,(2*i-1)], selectedx[,(2*i)], Ychat[,i], cex.lab=1.5, cex.axis = 1.3,
                     xlab=colnames(selectedx)[2*i-1], ylab=colnames(selectedx)[2*i], zlab=sety[i],
                     pch=16, color = "black", grid=F, cex=1.3, angle = 5)
s3d$plane3d(Intercept = 0, betaN[(2*i-1),], betaN[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(0.2, 0.6, 0.8, 0.2), border = NA))
# s3d$plane3d(Intercept = 0, betaT[(2*i-1),], betaT[(2*i),], draw_polygon = T, draw_lines = F,
#             polygon_args = list(col = rgb(1, 1, 0.6, 0.4), border = NA))
s3d$plane3d(Intercept = 0, beta[(2*i-1),], beta[(2*i),], draw_polygon = T, draw_lines = F,
            polygon_args = list(col = rgb(1, 0.7, 0.8, 0.4), border = NA))

dev.off()

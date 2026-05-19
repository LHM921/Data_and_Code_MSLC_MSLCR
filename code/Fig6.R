################################################################################
# 
# Purpose: produce boxplots of MSE for MSL and MSLC models under nu=4 and 10, as the sample size increases
# Input: read the 'result/SIM/Smse_all.csv' for 100 Monte Carlo datasets
# Package: dplyr
# Output: produce 'result/fig6_nu4.pdf' and 'result/fig6_nu10.pdf'
# 
################################################################################

#install.packages("dplyr")
library(dplyr)

##### 【Plot nu=4】##### 
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig6_nu4.pdf"),width=13.5, height=8)

layout(
  matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
           0,1:10,0,
           0,11:20,0,
           0,21:30,0,
           0,31:40,0,
           0,0,0,0,0,0,0,0,0,0,0,0), 6, 12, byrow = TRUE),
  heights=c(0.3,1,4,4,4,0.8),widths = c(0.04,0.33,rep(1.01,9),0.05)
)
# layout.show(40)

#####【nu=4】 #####
par(mar = c(0, 0, 1, 0))  
plot.new()  
text(x = 1, y =1.16, labels = expression(paste("MSE(", nu, "=4)")),cex =1.6, xpd = NA)
legend(0.35, 1.16, 
       xpd=NA, c('MSL','MSLC'),
       fill = c("pink", "skyblue3"), 
       bty='n', cex=1.2, x.intersp = 0.4)

par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[1]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[2]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[3]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[11]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[21]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[22]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[31]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[32]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[33]), cex.main = 3) 



##### 【Crate=0】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "0%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0
df_0 <- df_mse %>%
  filter(Crate==0)
# nu=4
df4 <- df_0 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df4_MSL <- df4 %>%
  filter(Model %in% c('MSL')) 
df4_MSLC <- df4 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error1~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error1~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu2 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error2~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error2~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu3 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error3~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error3~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma11 #####
par(mar = c(0.1, 0.1, 0.1, 0.1))
bp <- boxplot(Squared_Error4~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma21 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error5~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error5~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma22 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error6~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error6~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma31 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error7~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error7~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.01))

##### sigma32 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error8~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error8~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.023))

##### sigma33 #####
par(mar = c(0.1, 0, 0.1, 0.4))
bp <- boxplot(Squared_Error9~N, data=df4_MSL, plot = FALSE)
boxplot(Squared_Error9~N, data=df4_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))




##### 【Crate=0.05】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "5%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.05
df_5 <- df_mse %>%
  filter(Crate==0.05)
# nu=4
df4 <- df_5 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df4_MSL <- df4 %>%
  filter(Model %in% c('MSL')) 
df4_MSLC <- df4 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error1~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error1~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu2 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error2~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error2~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.02))

##### mu3 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error3~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error3~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma11 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error4~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma21 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error5~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error5~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma22 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error6~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error6~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))


##### sigma31 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error7~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error7~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma32 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error8~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error8~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma33 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error9~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error9~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))





##### 【Crate=0.15】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "15%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.15
df_15 <- df_mse %>%
  filter(Crate==0.15)
# nu=4
df4 <- df_15 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df4_MSL <- df4 %>%
  filter(Model %in% c('MSL')) 
df4_MSLC <- df4 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error1~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error1~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### mu2 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error2~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error2~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.2))
mtext("Sample Size", side = 1, line = 2.2)

##### mu3 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error3~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error3~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma11 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error4~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma21 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error5~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error5~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma22 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error6~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error6~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma31 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error7~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error7~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma32 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error8~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error8~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma33 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error9~N, data=df4_MSLC, plot = FALSE)
boxplot(Squared_Error9~N, data=df4_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

dev.off()





##### 【Plot nu=10】##### 
rm(list = ls())
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig6_nu10.pdf"),width=13.5, height=8)

layout(
  matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
           0,1:10,0,
           0,11:20,0,
           0,21:30,0,
           0,31:40,0,
           0,0,0,0,0,0,0,0,0,0,0,0), 6, 12, byrow = TRUE),
  heights=c(0.3,1,4,4,4,0.8),widths = c(0.04,0.33,rep(1.01,9),0.05)
)
# layout.show(40)

##### 【nu=10】 #####
par(mar = c(0, 0, 1, 0))  
plot.new()  
text(x = 1, y =1.16, labels = expression(paste("MSE(", nu, "=10)")),cex =1.6, xpd = NA)
legend(0.35, 1.16, 
       xpd=NA, c('MSL','MSLC'),
       fill = c("pink", "skyblue3"), 
       bty='n', cex=1.2, x.intersp = 0.4)

par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[1]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[2]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(mu[3]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[11]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[21]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[22]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[31]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[32]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[33]), cex.main = 3) 



##### 【Crate=0】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "0%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0
df_0 <- df_mse %>%
  filter(Crate==0)
# nu=10
df10 <- df_0 %>% 
  filter(nu ==10)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df10_MSL <- df10 %>%
  filter(Model %in% c('MSL')) 
df10_MSLC <- df10 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error1~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error1~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu2 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error2~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error2~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.03))

##### mu3 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error3~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error3~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.03))

##### sigma11 #####
par(mar = c(0.1, 0.1, 0.1, 0.1))
bp <- boxplot(Squared_Error4~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error4~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma21 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error5~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error5~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma22 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error6~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error6~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma31 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error7~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error7~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma32 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error8~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error8~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma33 #####
par(mar = c(0.1, 0, 0.1, 0.4))
bp <- boxplot(Squared_Error9~N, data=df10_MSL, plot = FALSE)
boxplot(Squared_Error9~N, data=df10_MSL, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))




##### 【Crate=0.05】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "5%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.05
df_5 <- df_mse %>%
  filter(Crate==0.05)
# nu=10
df10 <- df_5 %>% 
  filter(nu ==10)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df10_MSL <- df10 %>%
  filter(Model %in% c('MSL')) 
df10_MSLC <- df10 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error1~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error1~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu2 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error2~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error2~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### mu3 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error3~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error3~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma11 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error4~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error4~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma21 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error5~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error5~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma22 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error6~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error6~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))


##### sigma31 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error7~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error7~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma32 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error8~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error8~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma33 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error9~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error9~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))





##### 【Crate=0.15】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "15%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/SIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error10, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.15
df_15 <- df_mse %>%
  filter(Crate==0.15)
# nu=10
df10 <- df_15 %>% 
  filter(nu ==10)%>% 
  filter(Model %in% c('MSL', 'MSLC'))
df10_MSL <- df10 %>%
  filter(Model %in% c('MSL')) 
df10_MSLC <- df10 %>%
  filter(Model %in% c('MSLC'))


##### mu1 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error1~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error1~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### mu2 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error2~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error2~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### mu3 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error3~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error3~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma11 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error4~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error4~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma21 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error5~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error5~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma22 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error6~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error6~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma31 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error7~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error7~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,2]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma32 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error8~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error8~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma33 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error9~N, data=df10_MSLC, plot = FALSE)
boxplot(Squared_Error9~N, data=df10_MSLC, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

dev.off()


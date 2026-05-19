################################################################################
# 
# Purpose: produce boxplots of MSE for MSLR and MSLCR models under nu=4, as the sample size increases
# Input: read the 'result/RSIM/Smse_all.csv' for 100 Monte Carlo datasets
# Package: dplyr
# Output: produce 'result/fig7_1.pdf' and 'result/fig7_2.pdf'
# 
################################################################################

# install.packages("dplyr")
library(dplyr)

##### 【Plot beta1~sigma33】##### 
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig7_1.pdf"), width=13.5, height=8)

layout(
  matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
           0,1:10,0,
           0,11:20,0,
           0,21:30,0,
           0,0,0,0,0,0,0,0,0,0,0,0), 5, 12, byrow = TRUE),
  heights=c(1.8,1,4,4,1.8),widths = c(0.05,0.34,rep(1,9),0.05)
)
# layout.show(30)

#####【nu=4】 #####
par(mar = c(0, 0, 1, 0))  
plot.new()  
text(x = 1.15, y =1.05, labels = expression(paste("MSE(", nu, "=4)")),cex =1.6, xpd = NA)
legend(0.4, 1.05, 
       xpd=NA, c('MSLR','MSLCR'),
       fill = c("pink", "skyblue3"), 
       bty='n', cex=1.2, x.intersp = 0.5)


par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(beta[1]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(beta[2]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(beta[3]), cex.main = 3) 
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

data <- read.csv(paste(MS.PATH,'result/RSIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error19, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0
df_0 <- df_mse %>%
  filter(Crate==0)
# nu=4
df4 <- df_0 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSLR', 'MSLCR'))
df4_MSLR <- df4 %>%
  filter(Model %in% c('MSLR')) 
df4_MSLCR <- df4 %>%
  filter(Model %in% c('MSLCR'))


##### beta1 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error1~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error1~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### beta2 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error2~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error2~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### beta3 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error3~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error3~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma11 #####
par(mar = c(0.1, 0.1, 0.1, 0.1))
bp <- boxplot(Squared_Error4~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma21 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error5~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error5~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma22 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error6~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error6~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma31 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error7~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error7~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma32 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error8~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error8~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma33 #####
par(mar = c(0.1, 0, 0.1, 0.4))
bp <- boxplot(Squared_Error9~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error9~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))




##### 【Crate=0.05】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "5%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/RSIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error19, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.05
df_5 <- df_mse %>%
  filter(Crate==0.05)
# nu=4
df4 <- df_5 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSLR', 'MSLCR'))
df4_MSLR <- df4 %>%
  filter(Model %in% c('MSLR')) 
df4_MSLCR <- df4 %>%
  filter(Model %in% c('MSLCR'))


##### beta1 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error1~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error1~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### beta2 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error2~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error2~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### beta3 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error3~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error3~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma11 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error4~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma21 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error5~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error5~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma22 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error6~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error6~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma31 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error7~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error7~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma32 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error8~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error8~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma33 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error9~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error9~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

dev.off()



##### 【Plot sigma41~sigma55】##### 
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
cairo_pdf(file = file.path(MS.PATH, "result", "fig7_2.pdf"), width=13.5, height=8)

layout(
  matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
           0,1:10,0,
           0,11:20,0,
           0,21:30,0,
           0,0,0,0,0,0,0,0,0,0,0,0), 5, 12, byrow = TRUE),
  heights=c(1.8,1,4,4,1.8),widths = c(0.05,0.34,rep(1,9),0.05)
)
# layout.show(30)

#####【nu=4】 #####
par(mar = c(0, 0, 1, 0))  
plot.new()  
text(x = 1.15, y =1.05, labels = expression(paste("MSE(", nu, "=4)")),cex =1.6, xpd = NA)
legend(0.4, 1.05, 
       xpd=NA, c('MSLR','MSLCR'),
       fill = c("pink", "skyblue3"), 
       bty='n', cex=1.2, x.intersp = 0.5)


par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[41]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[42]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[43]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[44]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[51]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[52]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[53]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[54]), cex.main = 3) 
par(mar = c(0, 0, 4, 0))  
plot.new()                
title(expression(sigma[55]), cex.main = 3) 



##### 【Crate=0】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "0%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/RSIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error19, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0
df_0 <- df_mse %>%
  filter(Crate==0)
# nu=4
df4 <- df_0 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSLR', 'MSLCR'))
df4_MSLR <- df4 %>%
  filter(Model %in% c('MSLR')) 
df4_MSLCR <- df4 %>%
  filter(Model %in% c('MSLCR'))


##### sigma41 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error10~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error10~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma42 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error11~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error11~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma43 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error12~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error12~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma44 #####
par(mar = c(0.1, 0.1, 0.1, 0.1))
bp <- boxplot(Squared_Error13~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error4~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma51 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error14~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error14~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma52 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error15~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error15~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))


##### sigma53 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error16~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error16~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma54 #####
par(mar = c(0.1, 0, 0.1, 0.1))
bp <- boxplot(Squared_Error17~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error17~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))

##### sigma55 #####
par(mar = c(0.1, 0, 0.1, 0.4))
bp <- boxplot(Squared_Error18~N, data=df4_MSLR, plot = FALSE)
boxplot(Squared_Error18~N, data=df4_MSLR, col = "pink"
        , border = "red4",  xlab = "", ylab = "", xaxt = "n", boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))




##### 【Crate=0.05】 #####
par(mar = c(0, 0, 7, 0.1))  
plot.new()                
text(x = 0.2, y = 0.9, labels = "5%",
     cex = 2, srt = 90, xpd = TRUE)

data <- read.csv(paste(MS.PATH,'result/RSIM/Smse_all.csv', sep=''))
df_mse <- subset(data, select = -c(Squared_Error19, MSE))
df_mse$Crate <- as.numeric(df_mse$Crate)
# censoring =0.05
df_5 <- df_mse %>%
  filter(Crate==0.05)
# nu=4
df4 <- df_5 %>% 
  filter(nu ==4)%>% 
  filter(Model %in% c('MSLR', 'MSLCR'))
df4_MSLR <- df4 %>%
  filter(Model %in% c('MSLR')) 
df4_MSLCR <- df4 %>%
  filter(Model %in% c('MSLCR'))

##### sigma41 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error10~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error10~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma42 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error11~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error11~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma43 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error12~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error12~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma44 #####
par(mar = c(0, 0.1, 0.2, 0.1))
bp <- boxplot(Squared_Error13~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error13~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma51 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error14~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error14~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma52 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error15~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error15~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma53 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error16~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error16~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma54 #####
par(mar = c(0, 0, 0.2, 0.1))
bp <- boxplot(Squared_Error17~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error17~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

##### sigma55 #####
par(mar = c(0, 0, 0.2, 0.4))
bp <- boxplot(Squared_Error18~N, data=df4_MSLCR, plot = FALSE)
boxplot(Squared_Error18~N, data=df4_MSLCR, col = "skyblue3"
        , border = "blue4",  xlab = "", ylab = "",  boxwex = 0.66, xlim=c(0.5,4.5)
        , ylim = c(bp$stats[1,1]*0.995, bp$stats[5,1]*1.005))
mtext("Sample Size", side = 1, line = 2.2)

dev.off()
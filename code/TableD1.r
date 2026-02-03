################################################################################
# 
# Purpose: produce Supplementary Table D.1
#        perform model fitting of the MVN, MVT, MSL, MVNC, MVTC and MSLC models to the bioassay data  
# Input: load 'result/bioassay_sep.RData' and 'result/bioassay_nov.RData' file
# Output: print Supplementary Table D.1 (September part)
#         print Supplementary Table D.1 (November part)
# 
################################################################################

# September ############################################
MS.PATH = paste(getwd(),"/Data and Code/",sep="")
load(paste(MS.PATH, 'result/bioassay_sep.RData',sep=''))

Table1b = cbind(t(estN$IM$mu.hat), t(estT$IM$mu.hat), t(estS$IM$mu.hat),
                t(estNC$IM$mu.hat), t(estTC$IM$mu.hat), t(estSC$IM$mu.hat))
print(round(Table1b, 3))

Table1c = cbind(c(estN$EST[-c(1:ncol(Yc))],0), c(estT$EST[-c(1:ncol(Yc))]), c(estS$EST[-c(1:ncol(Yc))]),
                c(estNC$EST[-c(1:ncol(Yc))],0), c(estTC$EST[-c(1:ncol(Yc))]), c(estSC$EST[-c(1:ncol(Yc))]))
print(round(Table1c, 3))


# November ############################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")
load(paste(MS.PATH, 'result/bioassay_nov.RData',sep=''))

Table1b = cbind(t(estN$IM$mu.hat), t(estT$IM$mu.hat), t(estS$IM$mu.hat),
                t(estNC$IM$mu.hat), t(estTC$IM$mu.hat), t(estSC$IM$mu.hat))
print(round(Table1b, 3))

Table1c = cbind(c(estN$EST[-c(1:ncol(Yc))],0), c(estT$EST[-c(1:ncol(Yc))]), c(estS$EST[-c(1:ncol(Yc))]),
                c(estNC$EST[-c(1:ncol(Yc))],0), c(estTC$EST[-c(1:ncol(Yc))]), c(estSC$EST[-c(1:ncol(Yc))]))
print(round(Table1c, 3))

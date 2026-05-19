################################################################################
# 
# Purpose: produce Table 5
#      perform simulation results for the MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models  
# Input: read the 'Sfit_all.csv' for 100 Monte Carlo datasets
# Output: print Table 5: Averages of AIC and BIC scores and frequencies
# 
################################################################################

MS.PATH = paste(getwd(),"/Data and Code/",sep="")

# Read All Simulation results -------------- 
Sfit <- read.csv(paste(MS.PATH,'result/RSIM/Sfit_all.csv', sep=''))
Sfit$Model <- factor(Sfit$Model, levels = c('MVNR', 'MVTR', 'MSLR', 'MVNCR', 'MVTCR', 'MSLCR'))

# Sfit -------------- 
table_AIC <- aggregate(AIC ~ N + Crate + Model + nu, data = Sfit, FUN = function(x) round(mean(x), 3))
table_AIC <- table_AIC[order(table_AIC$Crate, table_AIC$nu, table_AIC$N, table_AIC$Model), ]
table_BIC <- aggregate(BIC ~ N + Crate + Model + nu, data = Sfit, FUN = function(x) round(mean(x), 3))
table_BIC <- table_BIC[order(table_BIC$Crate, table_BIC$nu, table_BIC$N, table_BIC$Model), ]

print(cbind(table_AIC[table_AIC$Crate == 0 & table_AIC$Model == c('MVNR', 'MVTR', 'MSLR'), ],
            table_BIC[table_BIC$Crate == 0 & table_BIC$Model == c('MVNR', 'MVTR', 'MSLR'), "BIC", drop = FALSE],
            table_AIC[table_AIC$Crate == 0.05 & table_AIC$Model == c('MVNCR', 'MVTCR', 'MSLCR'), c("Model","Crate", "AIC")],
            table_BIC[table_BIC$Crate == 0.05 & table_BIC$Model == c('MVNCR', 'MVTCR', 'MSLCR'), "BIC", drop = FALSE]))

# Sfit freq--------------
freq_model <- function(Sfit, models, value, colname) {
  Sfit_sub <- Sfit[Sfit$Model %in% models, ]
  Sfit_sub$Model <- factor(Sfit_sub$Model, levels = models)
  comb <- expand.grid(
    N = unique(Sfit_sub$N),
    Crate = unique(Sfit_sub$Crate),
    nu = unique(Sfit_sub$nu),
    Model = models
  )
  comb[[colname]] <- 0
  
  for(r in unique(Sfit_sub$Rep)) {
    Sfit_r <- Sfit_sub[Sfit_sub$Rep == r, ]
    groups <- split(Sfit_r, list(Sfit_r$Crate, Sfit_r$N, Sfit_r$nu), drop = TRUE)
    
    for(g in groups) {
      if(nrow(g) == 0) next
      min_value <- min(g[[value]], na.rm = TRUE)
      models_min <- g$Model[g[[value]] == min_value]
      
      for(m in models_min) {
        idx <- which(comb$N == g$N[1] &
                       comb$Crate == g$Crate[1] &
                       comb$nu == g$nu[1] &
                       comb$Model == m)
        comb[[colname]][idx] <- comb[[colname]][idx] + 1
      }
    }
  }
  comb <- comb[order(comb$Crate, comb$nu, comb$N, comb$Model), ]
  
  return(comb)
}

freq_AIC <- freq_model(Sfit, c('MVNR', 'MVTR', 'MSLR'), "AIC", "AIC_freq")
freq_BIC <- freq_model(Sfit, c('MVNR', 'MVTR', 'MSLR'), "BIC", "BIC_freq")
freqc_AIC <- freq_model(Sfit, c('MVNCR', 'MVTCR', 'MSLCR'), "AIC", "AIC_freq")
freqc_BIC <- freq_model(Sfit, c('MVNCR', 'MVTCR', 'MSLCR'), "BIC", "BIC_freq")

print(cbind(freq_AIC[freq_AIC$Crate == 0 & freq_AIC$Model == c('MVNR', 'MVTR', 'MSLR'), ],
            freq_BIC[freq_BIC$Crate == 0 & freq_BIC$Model == c('MVNR', 'MVTR', 'MSLR'), "BIC_freq", drop = FALSE],
            freqc_AIC[freqc_AIC$Crate == 0.05 & freqc_AIC$Model == c('MVNCR', 'MVTCR', 'MSLCR'), c("Model","Crate", "AIC_freq")],
            freqc_BIC[freqc_BIC$Crate == 0.05 & freqc_BIC$Model == c('MVNCR', 'MVTCR', 'MSLCR'), "BIC_freq", drop = FALSE]))

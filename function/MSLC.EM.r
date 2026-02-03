# Basic functions:
tr = function(M)  sum(diag(M))
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
as.vech = function(M)
{
    Dim = dim(M)[1]
    posi = vech.posi(Dim)
    temp = paste('M', posi[,1], posi[,2], sep = '')
    vech.M = matrix(M[posi], nrow=1)
    colnames(vech.M) = temp
    return(vech.M)
}

# generate random sample from multivariate normal independent distribution
rmsl = function(n, mu, Sigma, nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN','MSL'))
{
  distr = distr[1]
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
  Lam = diag(sqrt(diag(Sigma)))
  Lam.inv = diag(1/sqrt(diag(Sigma)))
  R = Lam.inv %*%Sigma%*% Lam.inv 
  }
  if(distr=='MVN') tau = rep(1, n)  
  if(distr=='MSL') tau = rbeta(n, shape1=nu/2, shape2=1)
  tau.sq = sqrt(tau)
  if(p==1){
    Y = mu + 1/tau.sq * Sd * rnorm(n, mean=0, sd=1) 
  } else  Y = matrix(rep(mu, n), ncol=n) + matrix(rep(1/tau.sq, each=p), ncol=n) * (Lam %*% t(rmvnorm(n, mean=rep(0,p), sigma=R))) 
  return(list(Y=t(Y), tau=tau))
}

rcen = function(Y, cen.rate)
{
  n = nrow(Y); p = ncol(Y)
  ubd = apply(Y, 2, quantile, prob=cen.rate)
  CEN = matrix(0, nrow=n, ncol=p)
  for(i in 1:p){
    Y[,i] = ifelse(Y[,i]<=ubd[i], ubd[i], Y[,i])
    CEN[,i] = ifelse(Y[,i]==ubd[i], 1, 0)
  } 
  return(list(Y=Y, CEN=CEN))
}

rcen1 = function(Y, cen.rate)
{
  n = nrow(Y); p = ncol(Y)
  ubd = apply(Y, 2, quantile, prob=cen.rate)[1]
  CEN = matrix(0, nrow=n, ncol=p)
  Y[,1] = ifelse(Y[,1]<=ubd, ubd, Y[,1])
  CEN[,1] = ifelse(Y[,1]==ubd, 1, 0)
  return(list(Y=Y, CEN=CEN))
}

rcen2 = function(Y, cen.rate)
{
  n = nrow(Y); p = ncol(Y)
  ubd = apply(Y, 2, quantile, prob=cen.rate)
  CEN = matrix(0, nrow=n, ncol=p)
  for(i in 1:2){
    Y[,i] = ifelse(Y[,i]<=ubd[i], ubd[i], Y[,i])
    CEN[,i] = ifelse(Y[,i]==ubd[i], 1, 0)
  } 
  return(list(Y=Y, CEN=CEN))
}

gaqh.SL = function(x, dd, a.ast2, mu2.1, R22.1)
{
  di = length(mu2.1)
  if(di==1) del.i = (x-mu2.1)^2/c(R22.1)
  else del.i = t(x-mu2.1)%*%solve(R22.1)%*%(x-mu2.1)
  tag = IGfn(a=dd, b=c(a.ast2+del.i)/2)
  return(tag)
}

MSLC.EM = function(Yc, cen, distr=c('MVN','MSL'), init.par=NULL, tol=1e-5, max.iter=1000, per=10, ITER=20, per.iter=2)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  begin = proc.time()[1]
  n = nrow(Yc)
  p = ncol(Yc)
  cho = seq(1, ITER, per.iter)
  no.cho = length(cho)
  vechS = vech.posi(p)
# initial values
  if(length(init.par) != 0){
   mu = init.par$mu
   Sigma = init.par$Sigma
   old.par = c(mu, Sigma[vechS])
   if(distr=='MSL'){
    nu = init.par$nu
    old.par = c(mu, Sigma[vechS], nu)
  }} else{
   mu = colMeans(Yc)
   Sigma = cov(Yc)#*(n-1)/n
   old.par = c(mu, Sigma[vechS])
   if(distr=='MSL'){
    nu = runif(1,1,10)
    old.par = c(mu, Sigma[vechS], nu)
  }}
  
  Ip = diag(p)
  po = p - rowSums(cen)
  cen.subj = which(rowSums(cen) != 0)
  Nc = length(cen.subj)
  obs.subj = which(rowSums(cen) == 0)
  No = n - Nc
  ind.cen = colSums(t(cen) * 2 ^ (1:p - 1))
  num.class.cen = length(unique(ind.cen))
  row.posi = O.list = C.list = as.list(numeric(num.class.cen))
  uni.ind = unique(ind.cen)
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  
# old observed-data log-likelihood
  w.den = matrix(NA, n, 1)
  det.o = rep(NA, n)
  delta.o = rep(NA, n)
  pdf = cdf = w.cdf.pdf = rep(1, n)
  cent = t(Yc) - mu
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind  = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = matrix(rep(mu, no.ind), nrow=p)
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% (mu + Sigma %*% Soo %*% cent.ind)
      Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
    }
    if(distr=='MVN'){
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s]] = pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o)
    }}
    if(distr=='MSL'){
      if(pjc == 0){
        pdf[ind] = IGfn(a=(p+nu)/2, b=delta.o[ind]/2)
        idx0 = which(delta.o==0)
        pdf[ind[idx0]] = 2/(p+nu)
      } 
      if(pjc == 1){  
        for(s in 1: no.ind) cdf[ind[s]] = integrate(gaqh.SL, lower=-Inf, upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=mu.co[,s], R22.1=Scc.o)$value
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = cubintegrate(gaqh.SL, lower=rep(-Inf, pjc), upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=c(mu.co[,s]), R22.1=Scc.o)$integral 
    }}
    }
    if(distr=='MVN') w.den = -.5*po*log(2*pi) -.5*log(det.o) -.5*delta.o + log(cdf)
    if(distr=='MSL') w.den = log(cdf) + log(pdf) + log(nu) - log(2) - .5*p*log(2*pi) -.5*log(det.o)
    
    loglik.old = iter.lnL = sum(w.den)
    iter = 0
    if(distr == 'MVN') iter.EST = c(iter, mu, Sigma[vech.posi(p)])
    if(distr == 'MSL') iter.EST = c(iter, mu, Sigma[vech.posi(p)], nu)

    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the MNIC with ", distr, "distribution with censoring proportion ", mean(cen)*100, '%.', "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")
    tau = rep(1,n)
    yhat = matrix(NA, nrow=n, ncol=p)         # add yhat
    tauyhat = matrix(NA, n, p)
    tauy2hat = array(NA, dim=c(p, p, n))
    repeat
    {
      iter = iter + 1 
#E-step
      for(j in 1:num.class.cen){
        Cj = C.list[[j]]
        ind = row.posi[[j]]               
        no.ind = length(ind)
        pjc = nrow(Cj)     #p-pjo
        
        if(pjc == 0){ 
          yhat[ind, ] = Yc[ind, ]   # add yhat
          cent.ind = cent[,ind]
          d = diag(t(cent.ind)%*%solve(Sigma)%*%cent.ind)
          if(distr=='MSL'){
           tau[ind] = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
           idx0 = which(d==0)
           tau[ind[idx0]] = (p+nu)/(p+nu+2)
          }
          tauyhat[ind,] = tau[ind]*Yc[ind, ]
          for(s in 1: no.ind) tauy2hat[,,ind[s]] = tau[ind[s]]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
        } else{

          yc.hat = tau.yc.hat = matrix(NA, nrow=pjc, ncol=no.ind)
          yc2.hat = tau.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind))          
          if(pjc == p){
            mu.co = matrix(rep(mu, length(ind)), nrow=p)
            Scc.o = Sigma

            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
               yc.hat[,s] = EX$EY
               tauy2hat[,,ind[s]] = EX$EYY  
             }
             tauyhat[ind,] = t(t(Cj) %*% yc.hat)
            }
            if(distr=='MSL'){
              for(s in 1: no.ind){
               yjc = c(Cj%*%Yc[ind[s], ]) 
               tau[ind[s]] = nu/(nu+2) * pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu+2, distr='MSL')/pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')
               EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, nu=nu+2, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
               yc.hat[,s] = EX$EY
               tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY  
              }
              tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat)) 
            }  
            yhat[ind, ] = t(yc.hat)      # add yhat
            }else{
            O = O.list[[j]]
            OSO = O %*% Sigma %*% t(O)
            Soo = t(O) %*% solve(OSO) %*% O
            cent.ind = cent[, ind]
            mu.co = Cj %*% (mu + Sigma %*% Soo %*% cent.ind)
            Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
            d=diag(t(cent.ind)%*%Soo%*%cent.ind)
            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
               yc.hat[,s] = EX$EY
               yc2.hat[,,s] = EX$EYY
               yjo = O %*% Yc[ind[s],]
               tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s]) %*% Cj + t(Cj) %*% yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s] %*% Cj)
              }
            tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)) 
            yhat[ind, ] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)       # add yhat
            } 
            if(distr=='MSL'){
              nu.rate = nu/(nu+2)
              for(s in 1: no.ind){
                yjc = c(Cj%*%Yc[ind[s], ]) 
                if(pjc == 1){ 
                 int.sl = integrate(SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value       # Modify
                 int.fgsl = integrate(Fgr.SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value # Modify 
               } else{
                 int.sl = cubintegrate(SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral       # Modify
                 int.fgsl = cubintegrate(Fgr.SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral # Modify
               } 
               tau[ind[s]] = nu.rate * int.fgsl / int.sl
               yjo = O %*% Yc[ind[s],]
               M = 2*ITER
               yjc.gen = matrix(NA, nrow=pjc, ncol=M)
               yjc.gen[,1] = yjc
               # Modify
               for(m in 2: M) yjc.gen[,m] = MH.MSL.yjc(yjc.gen[,(m-1)], mu.co.j=mu.co[,s], Scc.o=Scc.o, nu=nu, yj=Yc[ind[s],], yjo=yjo, yjc=yjc, do.j=d[s])
               yjc.samp = matrix(yjc.gen[,-c(1:ITER)], nrow=pjc)[, cho] 
               cent.yc = matrix((yjc.samp - mu.co[,s]), nrow=pjc)
               d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
               bi = (d[s]+d.co)/2 
               IGrate = IGfn(a=((p+nu)/2)+1, b=bi)/IGfn(a=(p+nu)/2, b=bi)
               if(pjc == 1){ 
                yc.hat[,s] = mean(yjc.samp)            # add yhat
                tau.yc.hat[,s] = mean(IGrate * yjc.samp)
                tau.yc2.hat[,,s] = mean(IGrate * yjc.samp^2)
               } else{
                yc.hat[,s] = rowMeans(yjc.samp)        # add yhat
                tau.yc.hat[,s] = rowMeans(rep(IGrate, each=pjc) * yjc.samp)
                tauyc2 = array(NA, dim=c(pjc, pjc, no.cho))
                for(m in 1: no.cho) tauyc2[,,m] = IGrate[m] * yjc.samp[,m] %*% t(yjc.samp[,m])
                tau.yc2.hat[,,s] = apply(tauyc2, 1:2, mean)
               }
               tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(tau.yc.hat[,s]) %*% Cj) + (t(Cj) %*% tau.yc.hat[,s] %*% t(yjo) %*% O) + (t(Cj) %*% tau.yc2.hat[,,s] %*% Cj)
              }
            tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)))) + t(t(Cj) %*% tau.yc.hat)
            yhat[ind, ] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)         # add yhat
            } 
          }
      }}    
# M-step
      tauyhat = matrix(as.numeric(tauyhat),n,p)
      mu = colSums(tauyhat)/sum(tau)
      sum.Y2 = apply(tauy2hat, 1:2, sum) 
      Sigma = (sum.Y2 - colSums(tauyhat)%*%t(mu) - mu%*%t(colSums(tauyhat)) + sum(tau)*mu%*%t(mu)) / n
      new.par = c(mu, Sigma[vechS])
       # Modify
      if(distr=='MSL'){
        nu = nlminb(start = nu, objective = MSLC.nu.fn, lower = 2+1e-6, upper =Inf, mu=mu, Sig=Sigma, Yc=Yc, cen=cen)$par  
        new.par = c(mu, Sigma[vechS], nu)
      }

# new observed-data log-likelihood
  cent = t(Yc) - mu
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind  = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = matrix(rep(mu, no.ind), nrow=p)
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% (mu + Sigma %*% Soo %*% cent.ind)
      Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
    }
    if(distr=='MVN'){
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s]] = pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o)
    }}
    if(distr=='MSL'){
      if(pjc == 0){
        pdf[ind] = IGfn(a=((p+nu)/2), b=delta.o[ind]/2)
        idx0 = which(delta.o==0)
        pdf[ind[idx0]] = 2/(p+nu)
      } 
      if(pjc == 1){  
        for(s in 1: no.ind) cdf[ind[s]] = integrate(gaqh.SL, lower=-Inf, upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=mu.co[,s], R22.1=Scc.o)$value
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = cubintegrate(gaqh.SL, lower=rep(-Inf, pjc), upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=c(mu.co[,s]), R22.1=Scc.o)$integral 
      }}
    }
    if(distr=='MVN') w.den = -.5*po*log(2*pi) -.5*log(det.o) -.5*delta.o + log(cdf)
    if(distr=='MSL') w.den = log(cdf) + log(pdf) + log(nu) -log(2) - .5*p*log(2*pi) -.5*log(det.o)
    loglik.new = sum(w.den)
    diff = loglik.new - loglik.old
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)]))
    if(distr == 'MSL') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)], nu))

    # if (iter%%per == 0){
    #   if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t mu=', mu, '\t sigma=', Sigma[vech.posi(p)], "\n")
    #   if(distr=='MSL') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t nu=', nu, '\t mu=', mu, '\t sigma=', Sigma[vech.posi(p)], "\n")
    # }
    if(abs(diff) < tol | diff.par < tol | iter > max.iter) break
    loglik.old = loglik.new
    old.par = new.par
    }
    end = proc.time()[1]
    cat(paste(rep("-", 50), collapse = ""), "\n")
    run.sec = end - begin
    cat("The CPU time takes", run.sec, "seconds.\n")
    m = ncol(iter.EST)-1 
    m1 = p*(p+1)/2
    if(diff < 0){ 
      aic = 2 * m - 2 * loglik.old
      bic = m * log(n) - 2 * loglik.old
      model.inf = list(m = m, loglik = loglik.old, aic = aic, bic = bic, iter=iter)
      EST = iter.EST[iter, -1]
      mu.hat = iter.EST[iter, 2:(p+1)]
      signu.hat = iter.EST[iter, -c(1:(p+1))]
      vechS = vech.posi(p); vechS1 = vechS[,c(2,1)]
      Sigma[vechS] = Sigma[vechS1] = signu.hat[1:m1]
      if(distr == 'MVN') para.est = list(mu=mu.hat, Sigma=Sigma)
      if(distr == 'MSL') para.est = list(mu=mu.hat, Sigma=Sigma, nu=signu.hat[(m1+1)])
      cat("iter",iter, ":", "observed-data log-likelihood =", loglik.old, "\t diff =", diff, '\t mu=', mu.hat, '\t Sigma nu=', signu.hat, "\n") 
    } else{
     cat("iter",iter, ":", "observed-data log-likelihood =", loglik.new, "\t diff =", diff, '\t mu=', mu, '\t Sigma=', Sigma[vech.posi(p)], "\n")
     aic = 2 * m - 2 * loglik.new
     bic = m * log(n) - 2 * loglik.new
     model.inf = list(m = m, loglik = loglik.new, aic = aic, bic = bic, iter=iter)
     if(distr == 'MVN'){
       EST = c(mu, Sigma[vech.posi(p)])
       para.est = list(mu=mu, Sigma=Sigma) 
     } 
     if(distr == 'MSL'){
       EST = c(mu, Sigma[vech.posi(p)], nu)
       para.est = list(mu=mu, Sigma=Sigma, nu=nu)
     }
   }
   IM = I.MSLC(para.est, Yc, cen=cen, distr=distr)
   cat(paste(rep("-", 50), collapse = ""), "\n")
   return(list(run.sec = run.sec, model.inf = model.inf,tau=tau, para.est = para.est, EST = EST, yhat = yhat, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST))
}

# MSLC log-likelihood function for nu
MSLC.nu.fn = function(par, mu, Sig, Yc, cen)
{
  nu = par
  n = nrow(Yc)
  p = ncol(Yc)
  Ip = diag(p)
  po = p - rowSums(cen)
  ind.cen = colSums(t(cen) * 2 ^ (1:p - 1))
  num.class.cen = length(unique(ind.cen))
  row.posi = O.list = C.list = as.list(numeric(num.class.cen))
  uni.ind = unique(ind.cen)
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  
  w.den = matrix(NA, n, 1)
  det.o = rep(NA, n)
  delta.o = rep(NA, n)
  pdf = cdf = w.cdf.pdf = rep(1, n)
  cent = t(Yc) - mu
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind  = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = matrix(rep(mu, no.ind), nrow=p)
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% (mu + Sigma %*% Soo %*% cent.ind)
      Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
    }
    if(pjc == 0){
        pdf[ind] = IGfn(a=(p+nu)/2, b=delta.o[ind]/2)
        idx0 = which(delta.o==0)
        pdf[ind[idx0]] = 2/(p+nu)
    } 
    if(pjc == 1){  
        for(s in 1: no.ind) cdf[ind[s]] = integrate(gaqh.SL, lower=-Inf, upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=mu.co[,s], R22.1=Scc.o)$value
    }
    if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = cubintegrate(gaqh.SL, lower=rep(-Inf, pjc), upper=c(Cj%*%Yc[ind[s], ]), dd=((p+nu)/2), a.ast2=delta.o[ind[s]], mu2.1=c(mu.co[,s]), R22.1=Scc.o)$integral 
    }}
    w.den = log(cdf) + log(pdf) + log(nu) - log(2) - .5*p*log(2*pi) -.5*log(det.o)
    return(-sum(w.den))
}

# Louis' information matrix
I.MSLC = function(para.est, Yc, cen, distr=c('MVN','MSL'), ITER=20, per.iter=2)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  cho = seq(1, ITER, per.iter)
  no.cho = length(cho)
  n = nrow(Yc)
  p = ncol(Yc)
# parameter estimates
  mu = para.est$mu
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MSL') nu = para.est$nu 

  Ip = diag(p)
  po = p - rowSums(cen)
  cen.subj = which(rowSums(cen) != 0)
  Nc = length(cen.subj)
  obs.subj = which(rowSums(cen) == 0)
  No = n - Nc
  ind.cen = colSums(t(cen) * 2 ^ (1:p - 1))
  num.class.cen = length(unique(ind.cen))
  row.posi = O.list = C.list = as.list(numeric(num.class.cen))
  uni.ind = unique(ind.cen)
  cent = t(Yc) - mu

  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }

  tau = tau2 = rep(1,n)
  tauyhat = tau2yhat = matrix(NA, n, p)
  tauy2hat = tau2y2hat = array(NA, dim=c(p, p, n))
  J2 = array(NA, dim=c(p, p, n))
  J1J2 = array(NA, dim=c(p, p, n))
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]               
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
      
    if(pjc == 0){ 
      cent.ind = cent[,ind]
      d = diag(t(cent.ind)%*%Sig.inv%*%cent.ind)
      if(distr=='MSL'){
        tau[ind] = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
        idx0 = which(d==0)
        tau[ind[idx0]] = (p+nu)/(p+nu+2)
        tau2[ind] = tau[ind]^2
        G2G = IGfn(a=((p+nu)/2)+2, b=d/2) / IGfn(a=(p+nu)/2, b=d/2)
        for(s in 1: no.ind) J2[,,ind[s]] = (G2G[s] - tau[ind[s]]^2) * cent.ind[,s] %*% t(cent.ind[,s]) 
      }
      tauyhat[ind,] = tau[ind]*Yc[ind, ]
      tau2yhat[ind,] = tau2[ind]*Yc[ind, ] 
      for(s in 1: no.ind){
        tauy2hat[,,ind[s]] = tau[ind[s]]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
        tau2y2hat[,,ind[s]] = tau2[ind[s]]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
      }
    } else{
      yc.hat = tau.yc.hat = tau2.yc.hat = matrix(NA, nrow=pjc, ncol=no.ind)
      yc2.hat = tau.yc2.hat = tau2.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind))          
      if(pjc == p){
        mu.co = matrix(rep(mu, length(ind)), nrow=p)
        Scc.o = Sigma

        if(distr=='MVN'){
          for(s in 1: no.ind){
          EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
          yc.hat[,s] = EX$EY
          tauy2hat[,,ind[s]] = tau2y2hat[,,ind[s]] = EX$EYY  
          }
          tauyhat[ind,] = tau2yhat[ind,] = t(t(Cj) %*% yc.hat)
          }
        if(distr=='MSL'){
          for(s in 1: no.ind){
           yjc = c(Cj%*%Yc[ind[s], ]) 
           yjc.gen = t(rtmsl(ITER, mu, Sigma, nu=nu+2, distr='MSL', lower=rep(-Inf,pjc), upper=yjc)$Y)
           cent.yc = matrix((yjc.gen - mu), nrow=pjc)
           dc = diag(t(cent.yc)%*%Sig.inv%*%cent.yc)
           IGrate = IGfn(a=((p+nu)/2)+1, b=dc/2)/IGfn(a=(p+nu)/2, b=dc/2)
           IG2rate = IGfn(a=((p+nu)/2)+2, b=dc/2)/IGfn(a=(p+nu)/2, b=dc/2)             
           tau[ind[s]] = nu/(nu+2) * pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu+2, distr='MSL')/pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')
           tau2[ind[s]] = tau[ind[s]] * mean(IGrate) 
           EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, nu=nu+2, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
           yc.hat[,s] = EX$EY
           tauyhat[ind[s],] = tau[ind[s]]*(t(t(Cj) %*% yc.hat[,s])) 
           tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY  
           tau2yhat[ind[s],] = tau[ind[s]] * rowMeans(rep(IGrate, each=pjc)*yjc.gen)
           tau2yc2 = J2m = array(NA, dim=c(pjc, pjc, ITER))
           for(m in 1: ITER){
             tau2yc2[,,m] = IGrate[m] * yjc.gen[,m] %*% t(yjc.gen[,m])
             J2m[,,m] = (IG2rate[m]-IGrate[m]^2) * cent.yc[,m] %*% t(cent.yc[,m])
           } 
           tau2y2hat[,,ind[s]] = tau[ind[s]] * apply(tau2yc2, 1:2, mean)
           J2[,,ind[s]] = apply(J2m, 1:2, mean) 
        }}  
        }else{
          O = O.list[[j]]
          OSO = O %*% Sigma %*% t(O)
          Soo = t(O) %*% solve(OSO) %*% O
          cent.ind = cent[, ind]
          mu.co = Cj %*% (mu + Sigma %*% Soo %*% cent.ind)
          Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
          d = diag(t(cent.ind)%*%Soo%*%cent.ind)
          if(distr=='MVN'){
            for(s in 1: no.ind){
             EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
             yc.hat[,s] = EX$EY
             yc2.hat[,,s] = EX$EYY
             yjo = O %*% Yc[ind[s],]
             tauy2hat[,,ind[s]] = tau2y2hat[,,ind[s]] = (t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s]) %*% Cj + t(Cj) %*% yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s] %*% Cj)
            }
          tauyhat[ind,] = tau2yhat[ind,] = (t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)) 
          } 
          if(distr=='MSL'){
            nu.rate = nu/(nu+2)
            for(s in 1: no.ind){
              yjc = c(Cj%*%Yc[ind[s], ]) 
              if(pjc == 1){ 
               int.sl = integrate(SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value        # Modify
               int.fgsl = integrate(Fgr.SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value  # Modify
             } else{
               int.sl = cubintegrate(SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral       # Modify
               int.fgsl = cubintegrate(Fgr.SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral # Modify
             } 
             tau[ind[s]] = nu.rate * int.fgsl / int.sl
             yjo = O %*% Yc[ind[s],]
             M = 2*ITER
             yjc.gen = matrix(NA, nrow=pjc, ncol=M)
             yjc.gen[,1] = yjc
             # Modify
             for(m in 2: M) yjc.gen[,m] = MH.MSL.yjc(yjc.gen[,(m-1)], mu.co.j=mu.co[,s], Scc.o=Scc.o, nu=nu, yj=Yc[ind[s],], yjo=yjo, yjc=yjc, do.j=d[s])
             yjc.samp = matrix(yjc.gen[,-c(1:ITER)], nrow=pjc)[, cho] 
             cent.yc = matrix((yjc.samp - mu.co[,s]), nrow=pjc)
             d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
             bi = (d[s]+d.co)/2 
             IGrate = IGfn(a=((p+nu)/2)+1, b=bi)/IGfn(a=(p+nu)/2, b=bi)
             IG2rate = IGfn(a=((p+nu)/2)+2, b=bi)/IGfn(a=(p+nu)/2, b=bi)
             tau2[ind[s]] = mean(IGrate^2)
             if(pjc == 1){ 
              tau.yc.hat[,s] = mean(IGrate * yjc.samp)
              tau2.yc.hat[,s] = mean(IGrate^2 * yjc.samp)
              tau.yc2.hat[,,s] = mean(IGrate * yjc.samp^2)
              tau2.yc2.hat[,,s] = mean(IGrate^2 * yjc.samp^2) 
              J2m = array(NA, dim=c(p, p, no.cho))
              for(m in 1: no.cho) J2m[,,m] = (IG2rate[m]-IGrate[m]^2) * (t(O)%*%yjo+t(Cj)%*%yjc.samp[m]-mu) %*% t(t(O)%*%yjo+t(Cj)%*%yjc.samp[m]-mu) 
            } else{
              tau.yc.hat[,s] = rowMeans(rep(IGrate, each=pjc) * yjc.samp) 
              tau2.yc.hat[,s] = rowMeans(rep(IGrate^2, each=pjc) * yjc.samp)
              tauyc2 = tau2yc2 = array(NA, dim=c(pjc, pjc, no.cho))
              J2m = array(NA, dim=c(p, p, no.cho))
              for(m in 1: no.cho){
               tauyc2[,,m] = IGrate[m] * yjc.samp[,m] %*% t(yjc.samp[,m])
               tau2yc2[,,m] = IGrate[m]^2 * yjc.samp[,m] %*% t(yjc.samp[,m])
               J2m[,,m] = (IG2rate[m]-IGrate[m]^2) * (t(O)%*%yjo+t(Cj)%*%yjc.samp[,m]-mu) %*% t(t(O)%*%yjo+t(Cj)%*%yjc.samp[,m]-mu) 
              }
              tau.yc2.hat[,,s] = apply(tauyc2, 1:2, mean)
              tau2.yc2.hat[,,s] = apply(tau2yc2, 1:2, mean)
             }
             tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(tau.yc.hat[,s]) %*% Cj) + (t(Cj) %*% tau.yc.hat[,s] %*% t(yjo) %*% O) + (t(Cj) %*% tau.yc2.hat[,,s] %*% Cj)
             tau2y2hat[,,ind[s]] = tau2[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(tau2.yc.hat[,s]) %*% Cj) + (t(Cj) %*% tau2.yc.hat[,s] %*% t(yjo) %*% O) + (t(Cj) %*% tau2.yc2.hat[,,s] %*% Cj)
             J2[,,ind[s]] = apply(J2m, 1:2, mean)
            }
          tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)))) + t(t(Cj) %*% tau.yc.hat)
          tau2yhat[ind,] = tau2[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)))) + t(t(Cj) %*% tau2.yc.hat) 
          } 
  }}}    
  
  if(distr=='MVN'){
    for(j in 1: n) J1J2[,,j] = (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu) - mu%*%t(tau2yhat[j,]) + mu%*%t(mu)) - (tauyhat[j,]-mu) %*% t(tauyhat[j,]-mu) 
  }
  if(distr=='MSL'){
    for(j in 1: n) J1J2[,,j] = (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu) - mu%*%t(tau2yhat[j,]) + tau2[j]*mu%*%t(mu)) - (tauyhat[j,]-tau[j]*mu) %*% t(tauyhat[j,]-tau[j]*mu) + J2[,,j]  
  }
  
  I.mu = sum(tau)*Sig.inv - Sig.inv %*% apply(J1J2, 1:2, sum) %*% Sig.inv  
  sd.mu = sqrt(diag(solve(I.mu))) 
  return(list(I.mu = I.mu, SD=sd.mu, mu.hat=rbind(mu, sd.mu)))
}

# MSL integrate functions
SL = function(x, Sigma, nu, a.ast2, mu2.1, R22.1)
{
  p = nrow(Sigma)
  di = length(mu2.1)
  if(di==1) del.i = (x-mu2.1)^2/c(R22.1)
  else del.i = t(x-mu2.1)%*%solve(R22.1)%*%(x-mu2.1)
  tag = nu/(2*(2*pi)^(p/2)*sqrt(det(Sigma))) * IGfn(a=(p+nu)/2, b=c(a.ast2+del.i)/2)
  return(tag)
}

Fgr.SL = function(x, Sigma, nu, a.ast2, mu2.1, R22.1)
{
  p = nrow(Sigma)
  di = length(mu2.1)
  if(di==1) del.i = (x-mu2.1)^2/c(R22.1)
  else del.i = t(x-mu2.1)%*%solve(R22.1)%*%(x-mu2.1)
  delta = c(a.ast2+del.i)
  tag = nu/(2*(2*pi)^(p/2)*sqrt(det(Sigma))) * IGfn(a=(p+nu)/2, b=delta/2)
  Fg.tag = tag * pgamma(1, shape=(p+nu)/2, rate=delta)/pgamma(1, shape=((p+nu)/2)+1, rate=delta) 
  return(Fg.tag)
}

# Metropolis-Hastings
MH.MSL.yjc = function(yjc.old, mu.co.j, Scc.o, nu, yj, yjo, yjc, do.j)    # Modify 
{  
   p = length(yj)
   pjc = length(yjc)
   qyc.old = dtmsl(yjc.old, mu.co.j, Scc.o, nu, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
   pco.old = dyco.MSL(yjc.old, mu2.1=mu.co.j, R22.1=Scc.o, nu=nu, do.j=do.j, p=p)

   yjc.new = c(rtmsl(n=1, mu.co.j, Scc.o, nu, distr='MSL', lower=rep(-Inf, pjc), upper=yjc)$Y)
   qyc.new = dtmsl(yjc.new, mu.co.j, Scc.o, nu, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
   pco.new = dyco.MSL(yjc.new, mu2.1=mu.co.j, R22.1=Scc.o, nu=nu, do.j=do.j, p=p)
   
   log.accept = min(1, (qyc.old*pco.new) / (qyc.new*pco.old))
   a = runif(1)
   if(a < log.accept){ 
     yjc.new = yjc.new 
   } else yjc.new = yjc.old
   return(yjc.new)
}           

dyco.MSL = function(yjc, mu2.1, R22.1, nu, do.j, p)
{
  pjc = length(yjc)
  if(pjc==1){
   del.j = (yjc-mu2.1)^2/c(R22.1)
   intG = integrate(Gyc, lower=-Inf, upper=yjc, nu=nu, a.ast2=do.j, mu2.1=mu2.1, R22.1=R22.1, p=p)$value
  } else{
   del.j = t(yjc-mu2.1)%*%solve(R22.1)%*%(yjc-mu2.1)
   intG = cubintegrate(Gyc, lower=rep(-Inf, pjc), upper=yjc, nu=nu, a.ast2=do.j, mu2.1=mu2.1, R22.1=R22.1, p=p)$integral
  }
  tag = IGfn(a=(p+nu)/2, b=c(do.j+del.j)/2) / intG
  return(tag)
}

Gyc = function(x, nu, a.ast2, mu2.1, R22.1, p)
{
  di = length(mu2.1)
  if(di==1) del.i = (x-mu2.1)^2/c(R22.1)
  else del.i = t(x-mu2.1)%*%solve(R22.1)%*%(x-mu2.1)
  tag = IGfn(a=(p+nu)/2, b=c(a.ast2+del.i)/2)
  return(tag)
}


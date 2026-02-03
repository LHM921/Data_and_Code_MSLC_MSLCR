library(mvtnorm)
library(tmvtnorm)
library(mclust)

# generate MSLR data          
rmsl.reg = function(n, para, Xi)
{
  beta = para$beta
  Sigma = para$Sigma
  nu = para$nu
  p = nrow(Xi)
  tau = rbeta(n, shape1=nu/2, shape2=1)
  tau.sq = sqrt(tau)
  Ymat = matrix(rep(Xi %*% beta, n), nrow=n, byrow=T) + matrix(rep(1/tau.sq, each=p), nrow=n, byrow=T) * rmvnorm(n, mean=rep(0,p), sigma=Sigma)
  return(Ymat)
}

MSLCR.EM = function(Yc, X, cen, distr=c('MVN','MSL'), init.par=NULL, tol=1e-5, max.iter=1000, per=10, ITER=20, per.iter=2)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  begin = proc.time()[1]
  n = nrow(Yc)
  p = ncol(Yc)
  cho = seq(1, ITER, per.iter)
  no.cho = length(cho)
  vechS = vech.posi(p)
  yc = as.vector(t(Yc))

# initial values
  if(length(init.par) != 0){
   beta = init.par$beta
   Sigma = init.par$Sigma
   old.par = c(beta, Sigma[vechS])
   if(distr=='MSL'){
    nu = init.par$nu
    old.par = c(beta, Sigma[vechS], nu)
  }} else{
   beta = solve(t(X)%*%X) %*% t(X) %*% yc
   Sigma = cov(Yc)#*(n-1)/n
   old.par = c(beta, Sigma[vechS])
   if(distr=='MSL'){
    nu = runif(1,1,10)
    old.par = c(beta, Sigma[vechS], nu)
  }}
  
  Sig.inv = solve(Sigma)
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
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  
  cent = t(Yc) - mu
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind  = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = mu[, ind]
      if(no.ind==1) mu.co=as.matrix(mu.co) 
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
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

    if(distr == 'MVN') {
      iter.EST = c(iter, beta, Sigma[vech.posi(p)])
      old.par = c(beta, Sigma[vech.posi(p)])
    }
    if(distr == 'MSL') {
      iter.EST = c(iter, beta, Sigma[vech.posi(p)], nu)
      old.par = c(beta, Sigma[vech.posi(p)], nu)
    }

    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the MCR with ", distr, "distribution with censoring proportion ", mean(cen)*100, '%.', "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")
    tau = rep(1,n)
    yhat = matrix(NA, nrow=n, ncol=p)         
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
          yhat[ind, ] = Yc[ind, ]   
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
            mu.co = mu[, ind]
            if(no.ind==1) mu.co=as.matrix(mu.co) 
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
               
               tau[ind[s]] = nu/(nu+2) * pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu.co[,s], Sigma=Sigma, nu=nu+2, distr='MSL')/pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu.co[,s], Sigma=Sigma, nu=nu, distr='MSL')
               EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, nu=nu+2, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
               yc.hat[,s] = EX$EY
               tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY  
              }
              tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat)) 
            }  
            yhat[ind, ] = t(yc.hat)      
            }else{
            O = O.list[[j]]
            OSO = O %*% Sigma %*% t(O)
            Soo = t(O) %*% solve(OSO) %*% O
            cent.ind = cent[, ind]
            mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
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
            yhat[ind, ] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)       
            } 
            if(distr=='MSL'){
              nu.rate = nu/(nu+2)
              for(s in 1: no.ind){
                yjc = c(Cj%*%Yc[ind[s], ]) 
                if(pjc == 1){ 
                 int.sl = integrate(SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value        
                 int.fgsl = integrate(Fgr.SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value  
               } else{
                 int.sl = cubintegrate(SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral       
                 int.fgsl = cubintegrate(Fgr.SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral 
               } 
               tau[ind[s]] = nu.rate * int.fgsl / int.sl
               yjo = O %*% Yc[ind[s],]
               M = 2*ITER
               yjc.gen = matrix(NA, nrow=pjc, ncol=M)
               yjc.gen[,1] = yjc
               for(m in 2: M) yjc.gen[,m] = MH.MSL.yjc(yjc.gen[,(m-1)], mu.co.j=mu.co[,s], Scc.o=Scc.o, nu=nu, yj=Yc[ind[s],], yjo=yjo, yjc=yjc, do.j=d[s])
               yjc.samp = matrix(yjc.gen[,-c(1:ITER)], nrow=pjc)[, cho] 
               cent.yc = matrix((yjc.samp - mu.co[,s]), nrow=pjc)
               d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
               bi = (d[s]+d.co)/2 
               IGrate = IGfn(a=((p+nu)/2)+1, b=bi)/IGfn(a=(p+nu)/2, b=bi)
               if(pjc == 1){ 
                yc.hat[,s] = mean(yjc.samp)            
                tau.yc.hat[,s] = mean(IGrate * yjc.samp)
                tau.yc2.hat[,,s] = mean(IGrate * yjc.samp^2)
               } else{
                yc.hat[,s] = rowMeans(yjc.samp)        
                tau.yc.hat[,s] = rowMeans(rep(IGrate, each=pjc) * yjc.samp)
                tauyc2 = array(NA, dim=c(pjc, pjc, no.cho))
                for(m in 1: no.cho) tauyc2[,,m] = IGrate[m] * yjc.samp[,m] %*% t(yjc.samp[,m])
                tau.yc2.hat[,,s] = apply(tauyc2, 1:2, mean)
               }
               tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(tau.yc.hat[,s]) %*% Cj) + (t(Cj) %*% tau.yc.hat[,s] %*% t(yjo) %*% O) + (t(Cj) %*% tau.yc2.hat[,,s] %*% Cj)
              }
            tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)))) + t(t(Cj) %*% tau.yc.hat)
            yhat[ind, ] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)        
            } 
          }
      }}    
# M-step
      Ttauyhat = as.vector(t(tauyhat))
      TSig = kronecker(diag(n), Sigma)
      TSig.inv = kronecker(diag(n), Sig.inv)
      Tga = rep(tau, each=p)
      beta = solve(t(Tga * X) %*% TSig.inv %*% X) %*% (t(X) %*% TSig.inv %*% Ttauyhat)
      vec.mu = X %*% beta
      mu = matrix(vec.mu, nrow=p, ncol=n)
      
      E = array(NA, dim=c(p, p, n))
      for(i in 1: n) E[,,i] = tauy2hat[,,i] - tauyhat[i, ]%*%t(mu[,i]) - mu[,i]%*%t(tauyhat[i,]) + tau[i]*mu[,i]%*%t(mu[,i])
      Sigma = apply(E, 1:2, sum) / n
      Sig.inv = solve(Sigma)
      new.par = c(beta, Sigma[vechS])
      if(distr=='MSL'){
        nu =  nlminb(start = nu, objective = MSLCR.nu.fn, lower = 2+1e-6, upper =Inf, beta=beta, Sig=Sigma, Yc=Yc, X=X, cen=cen)$par     
        new.par = c(beta, Sigma[vechS], nu)
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
      mu.co = mu[, ind]
      if(no.ind==1) mu.co=as.matrix(mu.co) 
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
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
    # attention point
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vech.posi(p)]))
    if(distr == 'MSL') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vech.posi(p)], nu))
    
    # if (iter%%per == 0){
    #   if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t beta=', beta, "\n")
    #   if(distr=='MSL') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t nu=', nu, '\t beta=', beta, "\n")
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
    
    cat("iter",iter, ":", "observed-data log-likelihood =", loglik.new, "\t diff =", diff, '\t beta=', beta, '\t Sigma=', Sigma[vechS], "\n")
    aic = 2 * m - 2 * loglik.new
    bic = m * log(n) - 2 * loglik.new
    model.inf = list(m = m, loglik = loglik.new, aic = aic, bic = bic)
    if(distr == 'MVN'){
       EST = c(beta, Sigma[vechS])
       para.est = list(beta=beta, Sigma=Sigma) 
    } 
    if(distr == 'MSL'){
       EST = c(beta, Sigma[vechS], nu)
       para.est = list(beta=beta, Sigma=Sigma, nu=nu)
    }
    
    Yfit = X %*% beta
    e.hat = yc - Yfit
    Mtau = matrix(rep(1/tau, each=p), nrow=n, ncol=p, byrow=T)
    yhat = tauyhat * Mtau
    ychat = as.vector(t(yhat))[which(as.vector(t(cen))==1)]
    IM = I.MSLCR(para.est, Yc, X, cen=cen, distr=distr)
    cat(paste(rep("-", 50), collapse = ""), "\n")
    return(list(run.sec = run.sec, iter=iter, model.inf = model.inf, para.est = para.est, EST = EST, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST, yfit = Yfit, yoyc = yhat, ychat = ychat, error = e.hat))
}

# MSLCR log-likelihood function for nu
MSLCR.nu.fn = function(par, beta, Sig, Yc, X, cen)
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
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  
  cent = t(Yc) - mu
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind  = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = mu[, ind]
      if(no.ind==1) mu.co=as.matrix(mu.co) 
      Scc.o = Sigma
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
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
I.MSLCR = function(para.est, Yc, X, cen, distr=c('MVN','MSL'), ITER=20, per.iter=2)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  cho = seq(1, ITER, per.iter)
  no.cho = length(cho)
  n = nrow(Yc)
  p = ncol(Yc)
# parameter estimates
  
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MSL') nu = para.est$nu 
  
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)

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
  
  K12 = array(NA, dim=c(p, p, n))
  J2 = array(NA, dim=c(p, p, n))
  # J1J2 = array(NA, dim=c(p, p, n))
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
        
        mu.co = mu[, ind]
        if(no.ind==1) mu.co=as.matrix(mu.co) 
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
          
          mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
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
               int.sl = integrate(SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value                 
               int.fgsl = integrate(Fgr.SL, lower=-Inf, upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$value           
             } else{
               int.sl = cubintegrate(SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral        
               int.fgsl = cubintegrate(Fgr.SL, lower=rep(-Inf, pjc), upper=yjc, Sigma=Sigma, nu=nu, a.ast2=d[s], mu2.1=mu.co[,s], R22.1=Scc.o)$integral  
             } 
             tau[ind[s]] = nu.rate * int.fgsl / int.sl
             yjo = O %*% Yc[ind[s],]
             M = 2*ITER
             yjc.gen = matrix(NA, nrow=pjc, ncol=M)
             yjc.gen[,1] = yjc
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
              # for(m in 1: no.cho) J2m[,,m] = (IG2rate[m]-IGrate[m]^2) * (t(O)%*%yjo+t(Cj)%*%yjc.samp[m]-mu) %*% t(t(O)%*%yjo+t(Cj)%*%yjc.samp[m]-mu) 
            } else{
              tau.yc.hat[,s] = rowMeans(rep(IGrate, each=pjc) * yjc.samp) 
              tau2.yc.hat[,s] = rowMeans(rep(IGrate^2, each=pjc) * yjc.samp)
              tauyc2 = tau2yc2 = array(NA, dim=c(pjc, pjc, no.cho))
              J2m = array(NA, dim=c(p, p, no.cho))
              for(m in 1: no.cho){
               tauyc2[,,m] = IGrate[m] * yjc.samp[,m] %*% t(yjc.samp[,m])
               tau2yc2[,,m] = IGrate[m]^2 * yjc.samp[,m] %*% t(yjc.samp[,m])
               # J2m[,,m] = (IG2rate[m]-IGrate[m]^2) * (t(O)%*%yjo+t(Cj)%*%yjc.samp[,m]-mu) %*% t(t(O)%*%yjo+t(Cj)%*%yjc.samp[,m]-mu) 
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
    for(j in 1: n) K12[,,j] = (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu[,j]) - mu[,j]%*%t(tau2yhat[j,]) + mu[,j]%*%t(mu[,j])) - (tauyhat[j,]-mu[,j]) %*% t(tauyhat[j,]-mu[,j]) 
  }
  if(distr=='MSL'){
    for(j in 1: n) K12[,,j] = (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu[,j]) - mu[,j]%*%t(tau2yhat[j,]) + tau2[j]*mu[,j]%*%t(mu[,j])) - (tauyhat[j,]-tau[j]*mu[,j]) %*% t(tauyhat[j,]-tau[j]*mu[,j]) + J2[,,j]  
  }
  
  TSig = kronecker(diag(n), Sigma)
  TSig.inv = kronecker(diag(n), Sig.inv)
  Tga = rep(tau, each=p)
  np = n*p
  cumsum.np=seq(0, np, p)
  #  TK = matrix(0, np, np)
  #  for(i in 1: n){
  #   np.idx = (cumsum.np[i]+1): cumsum.np[(i+1)] 
  #   TK[np.idx, np.idx] = K12[,,i]
  #  } 
  #  I.beta = t(Tga * X) %*% TSig.inv %*% X - t(X) %*% TSig.inv %*% TK %*% TSig.inv %*% X  
  sumK2 = 0
  for(i in 1: n){
    np.idx = (cumsum.np[i]+1): cumsum.np[(i+1)]
    K12[,,i] <- replace(K12[,,i], is.na(K12[,,i]), 0) # add 
    sumK2 = sumK2 + t(X[np.idx, ]) %*% Sig.inv %*% K12[,,i] %*% Sig.inv %*% X[np.idx, ]
  }
  I.beta = t(Tga * X) %*% TSig.inv %*% X - sumK2 
  sd.beta = sqrt(diag(solve(I.beta))) 
  return(list(I.beta = I.beta, SD=sd.beta, beta.hat=rbind(c(beta), sd.beta)))
}

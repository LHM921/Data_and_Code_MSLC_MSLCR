library(mvtnorm)
library(tmvtnorm)
library(mclust)

# generate MVTR data
rmvt.reg = function(n, para, Xi)
{
  beta = para$beta 
  Sigma = para$Sigma 
  nu = para$nu
  p = nrow(Xi)
  tau = rgamma(n, shape=nu/2, rate=nu/2)
  tau.sq = sqrt(tau)
  Ymat = matrix(rep(Xi %*% beta, n), nrow=n, byrow=T) + matrix(rep(1/tau.sq, each=p), nrow=n, byrow=T) * rmvnorm(n, mean=rep(0,p), sigma=Sigma)
  return(Ymat)
}

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

MVTCR.EM = function(Yc, X, cen, distr=c('MVN','MVT'), init.par=NULL, tol=1e-5, max.iter=1000, per=10)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  begin = proc.time()[1]
  n = nrow(Yc)
  p = ncol(Yc)
  vechS = vech.posi(p)
  yc = as.vector(t(Yc))
# initial values
  if(length(init.par) != 0){
   beta = init.par$beta
   Sigma = init.par$Sigma
   old.par = c(beta, Sigma[vechS])
   if(distr=='MVT'){
    nu = init.par$nu
    old.par = c(beta, Sigma[vechS], nu)
  }} else{
   yc = as.vector(t(Yc))
   beta = solve(t(X)%*%X) %*% t(X) %*% yc
   Sigma = cov(Yc)#*(n-1)/n
   old.par = c(beta, Sigma[vechS])
   if(distr=='MVT'){
    nu = runif(1,2,50)
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
    if(distr=='MVT'){
      DF = round(nu+po[ind][1])
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s]] = pt(((Cj%*%Yc[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = ptmvt(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o,df=DF)
    }}
    }
    if(distr=='MVN') w.den = -.5*po*log(2*pi) -.5*log(det.o) -.5*delta.o + log(cdf)
    if(distr=='MVT') w.den = lgamma((nu+po)/2) - lgamma(nu/2) - .5*po*log(pi*nu)- .5*log(det.o)  - .5*(nu+po)*log(1+delta.o/nu) +log(cdf)
    
    loglik.old = iter.lnL = sum(w.den)
    iter = 0
    if(distr == 'MVN'){
      iter.EST = c(iter, beta, Sigma[vechS])
      old.par = c(beta, Sigma[vechS])
    }
    if(distr == 'MVT'){
     iter.EST = c(iter, beta, Sigma[vechS], nu)
     old.par = c(beta, Sigma[vechS], nu)
    }
    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the MCR with ", distr, "distribution under censoring proportion ", mean(cen)*100, '%.', "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")
    tau = rep(1,n)
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
          cent.ind = cent[,ind]
          d = diag(t(cent.ind) %*% Sig.inv %*% cent.ind)
          if(distr=='MVT') tau[ind] = (nu+p)/(nu+d)
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
            if(distr=='MVT'){
              for(s in 1: no.ind){
               yjc = c(Cj%*%Yc[ind[s], ])
               tau[ind[s]] = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=nu/(nu+2)*Sigma, df=round(nu+2))/ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=Sigma, df=round(nu))
               EX = TMSL.moment(mu=mu.co[,s], Sigma=nu/(nu+2)*Sigma, nu=nu+2, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
               yc.hat[,s] = EX$EY
               tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY  
              }
              tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat)) 
            }}else{
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
            } 
            if(distr=='MVT'){
              tmp = (nu + po[ind]) / (nu + d)
              DF = round(nu+po[ind][1])
              for(s in 1: no.ind){
               yjc = c(Cj%*%Yc[ind[s], ])
               newScc.o = rep(1/tmp[s], pjc) * Scc.o 
               if(pjc == 1){ cdf.rate = pt((yjc-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2) / pt((yjc-c(mu.co[,s]))/sqrt(newScc.o), df=DF)
               } else cdf.rate = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=rep((DF/(DF+2)),pjc)*newScc.o, df=DF+2) / ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=newScc.o, df=DF)
               tau[ind[s]] = tmp[s] * cdf.rate
               EX = TMSL.moment(mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
               yc.hat[,s] = EX$EY
               yc2.hat[,,s] = EX$EYY
               yjo = O %*% Yc[ind[s],]
               tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s]) %*% Cj + t(Cj) %*% yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s] %*% Cj)
              }
            tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)) 
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
      if(distr=='MVT'){
#        nu = optim(par = nu, fn = MVTCR.nu.fn, method = "L-BFGS-B", lower = 2+1e-6, upper =Inf, beta=beta, Sig=Sigma, Yc=Yc, X=X, cen=cen)$par
        nu = nlminb(start = nu, objective = MVTCR.nu.fn, lower = 2+1e-6, upper =Inf, beta=beta, Sig=Sigma, Yc=Yc, X=X, cen=cen)$par     
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
    if(distr=='MVT'){
      DF = round(nu+po[ind][1])
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s]] = pt(((Cj%*%Yc[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s]] = ptmvt(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o,df=DF)
      }}
    }
    if(distr=='MVN') w.den = -.5*po*log(2*pi) -.5*log(det.o) -.5*delta.o + log(cdf)
    if(distr=='MVT') w.den = lgamma((nu+po)/2) - lgamma(nu/2) - .5*po*log(pi*nu)- .5*log(det.o)  - .5*(nu+po)*log(1+delta.o/nu) +log(cdf)
    loglik.new = sum(w.den)
    diff = loglik.new - loglik.old
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vechS]))
    if(distr == 'MVT') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vechS], nu))

    # if (iter%%per == 0){
    # if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t beta=', beta, "\n")
    # if(distr=='MVT') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t nu=', nu, '\t beta=', beta, "\n")
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
    if(distr == 'MVT'){
      EST = c(beta, Sigma[vechS], nu)
      para.est = list(beta=beta, Sigma=Sigma, nu=nu)
    }
    Yfit = X %*% beta
    e.hat = yc - Yfit
    Mtau = matrix(rep(1/tau, each=p), nrow=n, ncol=p, byrow=T)
    yhat = tauyhat * Mtau
    ychat = as.vector(t(yhat))[which(as.vector(t(cen))==1)]
    IM = I.MVTCR(para.est, Yc, X, cen=cen, distr=distr)
    cat(paste(rep("-", 50), collapse = ""), "\n")
    return(list(run.sec = run.sec, iter=iter, model.inf = model.inf, para.est = para.est, EST = EST, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST, yfit = Yfit, yoyc = yhat, ychat = ychat, error = e.hat))
}

# MVT log-likelihood function for nu
MVTCR.nu.fn = function(par, beta, Sig, Yc, X, cen)
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
    no.ind = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){          
      det.o[ind] = 1
      delta.o[ind] = 0
      mu.co = mu[, ind]
      if(no.ind==1) mu.co=as.matrix(mu.co) 
      Scc.o = Sig
    } else{
      O = O.list[[j]]
      OSO = O %*% Sig %*% t(O)
      cent.ind = cent[, ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
      Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
    }
    DF = round(nu+po[ind][1])
    if(pjc == 1){
     for(s in 1: no.ind) cdf[ind[s]] = pt(((Cj%*%Yc[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
    }
    if(pjc >= 2){
      for(s in 1: no.ind) cdf[ind[s]] = ptmvt(lowerx=rep(-Inf, pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o, df=DF)
      # for(s in 1: length(ind)) cdf[ind[s]] = pmvt(lower=-Inf, upper=c(Cj%*%Yc[ind[s], ]), delta=c(mu.co[,s]), sigma=Scc.o, df=DF)
    }
    par.log.den = lgamma((nu+po)/2) - lgamma(nu/2) - .5*po*log(pi*nu) - .5*log(det.o) - .5*(nu+po)*log(1+delta.o/nu) +log(cdf) 
  }
  return(-sum(par.log.den))
}

# Louis' information matrix
I.MVTCR = function(para.est, Yc, X, cen, distr=c('MVN','MVT'))
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc)
  p = ncol(Yc)
# parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MVT') nu = para.est$nu 
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
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]               
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
      
    if(pjc == 0){ 
      cent.ind = cent[,ind]
      d = diag(t(cent.ind)%*%Sig.inv%*%cent.ind)
      if(distr=='MVT'){
        tau[ind] = (nu+p)/(nu+d)
        tau2[ind] = ((nu+p)/(nu+d))^2
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
        if(distr=='MVT'){
          for(s in 1: no.ind){
          yjc = c(Cj%*%Yc[ind[s], ])
          tau[ind[s]] = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=nu/(nu+2)*Sigma, df=round(nu+2))/ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=Sigma, df=round(nu))
          tau2[ind[s]] = ((nu+p)*(nu+2))/(nu*(p+nu+2)) * ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=nu/(nu+4)*Sigma, df=round(nu+4))/ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu.co[,s], sigma=Sigma, df=round(nu))
          EX = TMSL.moment(mu=mu.co[,s], Sigma=nu/(nu+2)*Sigma, nu=nu+2, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
          yc.hat[,s] = EX$EY
          EX2 = TMSL.moment(mu=mu.co[,s], Sigma=nu/(nu+4)*Sigma, nu=nu+4, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
          tau2.yc.hat[,s] = EX2$EY
          tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY 
          tau2y2hat[,,ind[s]] = tau2[ind[s]] * EX2$EYY 
         }
          tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat)) 
          tau2yhat[ind,] = tau2[ind]*(t(t(Cj) %*% tau2.yc.hat))
        }}else{
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
          if(distr=='MVT'){
            tmp = (nu + po[ind]) / (nu + d)
            tmp2 = (tmp * (p+nu) * (po[ind]+nu+2))/ ((p+nu+2)*(nu + d))
            DF = round(nu+po[ind][1])
            for(s in 1: no.ind){
             yjc = c(Cj%*%Yc[ind[s], ])
             newScc.o = rep(1/tmp[s], pjc) * Scc.o 
             if(pjc == 1){ 
               cdf.rate = pt((yjc-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2) / pt((yjc-c(mu.co[,s]))/sqrt(Scc.o), df=DF)
               cdf.rate2 = pt((yjc-c(mu.co[,s]))/sqrt((DF/(DF+4))*newScc.o), df=DF+4) / pt((yjc-c(mu.co[,s]))/sqrt(Scc.o), df=DF)
             } else{
               cdf.rate = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=rep((DF/(DF+2)),pjc)*newScc.o, df=DF+2) / ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=newScc.o, df=DF)
               cdf.rate2 = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=rep((DF/(DF+4)),pjc)*newScc.o, df=DF+4) / ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=c(mu.co[,s]), sigma=newScc.o, df=DF)
             }
             tau[ind[s]] = tmp[s] * cdf.rate
             tau2[ind[s]] = tmp2[s] * cdf.rate2 
             EX = TMSL.moment(mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
             yc.hat[,s] = EX$EY
             yc2.hat[,,s] = EX$EYY
             EX2 = TMSL.moment(mu=mu.co[,s], Sigma=rep(DF/(DF+4), pjc)*newScc.o, nu=DF+4, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
             tau2.yc.hat[,s] = EX2$EY
             tau2.yc2.hat[,,s] = EX2$EYY
             yjo = O %*% Yc[ind[s],]
             tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s]) %*% Cj + t(Cj) %*% yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s] %*% Cj)
             tau2y2hat[,,ind[s]] = tau2[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(tau2.yc.hat[,s]) %*% Cj + t(Cj) %*% tau2.yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% tau2.yc2.hat[,,s] %*% Cj)
          }
          tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat))
          tau2yhat[ind,] = tau2[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% tau2.yc.hat))
          }
  }}}    
  
  if(distr=='MVN'){
    for(j in 1: n) K12[,,j] = (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu[,j]) - mu[,j]%*%t(tau2yhat[j,]) + mu[,j]%*%t(mu[,j])) - (tauyhat[j,]-mu[,j]) %*% t(tauyhat[j,]-mu[,j]) 
  }
  if(distr=='MVT'){
    a1 = 1 + 2/(nu+p)
    for(j in 1: n) K12[,,j] = a1 * (tau2y2hat[,,j] - tau2yhat[j,]%*%t(mu[,j]) - mu[,j]%*%t(tau2yhat[j,]) + tau2[j]*mu[,j]%*%t(mu[,j])) - (tauyhat[j,]-tau[j]*mu[,j]) %*% t(tauyhat[j,]-tau[j]*mu[,j])  
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
   sumK2 = sumK2 + t(X[np.idx, ]) %*% Sig.inv %*% K12[,,i] %*% Sig.inv %*% X[np.idx, ]
  }
  I.beta = t(Tga * X) %*% TSig.inv %*% X - sumK2 
  sd.beta = sqrt(diag(solve(I.beta))) 
  return(list(I.beta = I.beta, SD=sd.beta, beta.hat=rbind(c(beta), sd.beta)))
}

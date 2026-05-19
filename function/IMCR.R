# ============================
# Louis' information matrix
# ============================
# MVNR
IM.MVNR = function(para.est, EST, Y, X)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Y)
  p = ncol(Y)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  
  # individual expected score vector
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = t(X[idx, ]) %*% Sig.inv %*% (Y[j, ] - X[idx, ] %*% beta)
    A2 = Y[j,] %*% t(Y[j,]) - Y[j,] %*% t(X[idx, ] %*% beta) - X[idx, ] %*% beta %*% t(Y[j,]) + X[idx, ] %*% beta %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}

# MVNCR
IM.MVNCR = function(para.est, EST, Yc, X, cen)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc)
  p = ncol(Yc)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p))
  Sig.inv = solve(Sigma)
  
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
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  yhat = array(NA, dim=c(n, p))
  y2hat = array(NA, dim=c(p, p, n))
  cent = array(NA, dim = c(p, n))
  cent = t(Yc) - mu
  
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]              
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
    
    if(pjc == 0){
      cent.ind = cent[, ind]
      d = diag(t(cent.ind)%*%solve(Sigma)%*%cent.ind)
      yhat[ind, ] = Yc[ind, ]
      for(s in 1: no.ind) y2hat[,,ind[s]] = Yc[ind[s], ] %*% t(Yc[ind[s], ])
    } else{
      yc.hat = array(NA, dim = c(pjc, no.ind))
      yc2.hat = array(NA, dim=c(pjc, pjc, no.ind))          
      if(pjc == p){
        mu.co = mu[, ind]
        Scc.o = Sigma
        for(s in 1: no.ind){
          EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
          yc.hat[,s] = EX$EY
          y2hat[,,ind[s]] = EX$EYY  
        }
        yhat[ind, ] = t(t(Cj) %*% yc.hat)
      }else{
        O = O.list[[j]]
        OSO = O %*% Sigma %*% t(O)
        Soo = t(O) %*% solve(OSO) %*% O
        cent.ind = cent[, ind]
        mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
        Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
        d=diag(t(cent.ind)%*%Soo%*%cent.ind)
        for(s in 1: no.ind){
          EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
          yc.hat[,s] = EX$EY
          yc2.hat[,,s] = EX$EYY
          yjo = O %*% Yc[ind[s],]
          y2hat[,,ind[s]] = t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s]) %*% Cj + t(Cj) %*% yc.hat[,s] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s] %*% Cj
        }
        yhat[ind,] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)
      }}} 
  
  # individual expected score vector
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = t(X[idx, ]) %*% Sig.inv %*% (yhat[j, ] - X[idx, ] %*% beta)
    A2 = y2hat[,,j] - yhat[j,] %*% t(X[idx, ] %*% beta) - X[idx, ] %*% beta %*% t(yhat[j,]) + X[idx, ] %*% beta %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}

# MVTR
IM.MVTR = function(para.est, EST, Y, X)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Y)
  p = ncol(Y)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p))
  Sig.inv = solve(Sigma)
  nu = para.est$nu
  
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  tau = ka = matrix(1, n)
  tauyhat = array(NA, dim=c(n, p))
  tauy2hat = array(NA, dim=c(p, p, n))
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  tau = (nu + p)/(nu + d)
  ka = digamma((nu + p)/2) - log((nu + d)/2)
  
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = 1 + xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = tau[j] * t(X[idx, ]) %*% Sig.inv %*% (Y[j, ]  - X[idx, ] %*% beta)
    A2 = (tau[j] * Y[j,]) %*% t(Y[j,]) - tau[j]*Y[j,] %*% t(X[idx, ] %*% beta) - (X[idx, ] %*% beta * tau[j]) %*% t(Y[j,]) + (tau[j] * X[idx, ] %*% beta) %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  EIs[(1+xp+q1), ] = .5 * (ka - tau - digamma(nu/2) - log(nu/2) -1)
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}

# MVTCR
IM.MVTCR = function(para.est, EST, Yc, X, cen)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc)
  p = ncol(Yc)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p))
  Sig.inv = solve(Sigma)
  nu = para.est$nu
  
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
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  tau = ka = matrix(1, n)
  tauyhat = array(NA, dim=c(n, p))
  tauy2hat = array(NA, dim=c(p, p, n))
  cent = array(NA, dim = c(p, n))
  cent = t(Yc) - mu
  
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]              
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
    
    if(pjc == 0){
      cent.ind = cent[, ind]
      d = diag(t(cent.ind)%*%solve(Sigma)%*%cent.ind)
      tau[ind] = (nu + p)/(nu + d)
      ka[ind] = digamma((nu + p)/2) - log((nu + d)/2)
      tauyhat[ind, ] = tau[ind] * Yc[ind, ]
      for(s in 1: no.ind) tauy2hat[,,ind[s]] = tau[ind[s]]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
    } else{
      yc.hat =  array(NA, dim = c(pjc, no.ind))
      yc2.hat = array(NA, dim=c(pjc, pjc, no.ind))          
      if(pjc == p){
        mu.co = mu[, ind]
        Scc.o = Sigma
        for(s in 1: no.ind){
          yjc = c(Cj%*%Yc[ind[s], ])
          tau[ind[s]] = ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu, sigma=nu/(nu+2)*Sigma, df=round(nu+2))/ptmvt(lowerx=rep(-Inf, pjc), upperx=yjc, mean=mu, sigma=Sigma, df=round(nu))
          EX = TMSL.moment(mu=mu.co[,s], Sigma=nu/(nu+2)*Sigma, nu=nu+2, distr='MVT', a.low=rep(-Inf, pjc), a.upp=yjc)
          yc.hat[,s] = EX$EY
          tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY 
          yjc.gen = t(rtmsl(10, mu=mu.co[,s], Sigma=nu/(nu+2)*Sigma, nu=round(nu+2), distr='MVT', lower=rep(-Inf,pjc), upper=yjc)$Y)
          cent.yc = matrix((yjc.gen - mu), nrow=pjc)
          dc = diag(t(cent.yc)%*%Sig.inv%*%cent.yc)
          ka[ind[s]] = digamma((nu+p)/2) - mean(log((nu+dc)/2))
        }
        tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat))
      }else{
        O = O.list[[j]]
        OSO = O %*% Sigma %*% t(O)
        Soo = t(O) %*% solve(OSO) %*% O
        cent.ind = cent[, ind]
        mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
        Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
        d=diag(t(cent.ind)%*%Soo%*%cent.ind)
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
          yjc.gen = t(rtmsl(10, mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', lower=rep(-Inf,pjc), upper=yjc)$Y)
          cent.yc = matrix((yjc.gen - mu.co[,s]), nrow=pjc)
          d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
          ka[ind[s]] = digamma((nu+p)/2) - mean(log((nu+d[s]+d.co)/2))  
        }
        tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat))
      }}}  
  # individual expected score vector
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = 1 + xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = t(X[idx, ]) %*% Sig.inv %*% (tauyhat[j, ] - tau[j]*X[idx, ] %*% beta)
    A2 = tauy2hat[,,j] - tauyhat[j,] %*% t(X[idx, ] %*% beta) - X[idx, ] %*% beta %*% t(tauyhat[j,]) + (tau[j] * X[idx, ] %*% beta) %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  EIs[(1+xp+q1), ] = .5 * (ka - tau - digamma(nu/2) - log(nu/2) -1)
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}

# MSLR
IM.MSLR = function(para.est, EST, Y, X)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Y)
  p = ncol(Y)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p))
  Sig.inv = solve(Sigma)
  nu = para.est$nu
  
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  tau = ka = matrix(1, n)
  yhat = matrix(NA, nrow=n, ncol=p)  
  tauyhat = array(NA, dim=c(n, p))
  tauy2hat = array(NA, dim=c(p, p, n))
  cent = array(NA, dim = c(p, n))
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  tau = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
  ka = digamma((nu + p)/2) - log((nu + d)/2)
  
  # individual expected score vector
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = 1 + xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = tau[j]*t(X[idx, ]) %*% Sig.inv %*% (Y[j, ]  - X[idx, ] %*% beta)
    A2 = (tau[j] * Y[j,]) %*% t(Y[j,]) - tau[j]*Y[j,] %*% t(X[idx, ] %*% beta) - (X[idx, ] %*% beta * tau[j]) %*% t(Y[j,]) + (tau[j] * X[idx, ] %*% beta) %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  EIs[(1+xp+q1), ] = (-1/nu) - .5 * ka
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}

# MSLCR
IM.MSLCR = function(para.est, EST, Yc, X, cen, ITER=20, per.iter=2)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc)
  p = ncol(Yc)
  cho = seq(1, ITER, per.iter)
  no.cho = length(cho)
  vechS = vech.posi(p)
  # parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p))
  Sig.inv = solve(Sigma)
  nu = para.est$nu
  
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
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  tau = ka = matrix(1, n)
  yhat = matrix(NA, nrow=n, ncol=p)  
  tauyhat = array(NA, dim=c(n, p))
  tauy2hat = array(NA, dim=c(p, p, n))
  cent = array(NA, dim = c(p, n))
  cent = t(Yc) - mu
  
  for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]              
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
    
    if(pjc == 0){
      cent.ind = cent[, ind]
      d = diag(t(cent.ind)%*%solve(Sigma)%*%cent.ind)
      tau[ind] = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
      ka[ind] = digamma((nu + p)/2) - log((nu + d)/2)
      idx0 = which(d==0)
      tau[ind[idx0]] = (p+nu)/(p+nu+2)
      tauyhat[ind, ] = tau[ind] * Yc[ind, ]
      for(s in 1: no.ind) tauy2hat[,,ind[s]] = tau[ind[s]]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
    } else{
      yc.hat = tau.yc.hat =  array(NA, dim = c(pjc, no.ind))
      yc2.hat = tau.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind))          
      if(pjc == p){
        mu.co = mu[, ind]
        Scc.o = Sigma
        
        for(s in 1: no.ind){
          yjc = c(Cj%*%Yc[ind[s], ]) 
          tau[ind[s]] = nu/(nu+2) * pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu+2, distr='MSL')/pmsl(lower = rep(-Inf, pjc), upper = yjc, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')
          EX = TMSL.moment(mu=mu.co[,s], Sigma=Scc.o, nu=nu+2, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
          yc.hat[,s] = EX$EY
          tauy2hat[,,ind[s]] = tau[ind[s]] * EX$EYY  
          dc = diag(t(cent.yc)%*%Sig.inv%*%cent.yc)
          ka[ind[s]] = digamma((nu+p)/2) - mean(log((nu+dc)/2))
        }
        tauyhat[ind,] = tau[ind]*(t(t(Cj) %*% yc.hat)) 
      }else{
        O = O.list[[j]]
        OSO = O %*% Sigma %*% t(O)
        Soo = t(O) %*% solve(OSO) %*% O
        cent.ind = cent[, ind]
        mu.co = Cj %*% mu[, ind] + Cj %*% Sigma %*% Soo %*% cent.ind
        Scc.o = Cj %*% (Ip - Sigma %*% Soo) %*% Sigma %*% t(Cj)
        d=diag(t(cent.ind)%*%Soo%*%cent.ind)
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
          d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
          ka[ind[s]] = digamma((nu+p)/2) - mean(log((nu+d[s]+d.co)/2))  
          tauy2hat[,,ind[s]] = tau[ind[s]]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(tau.yc.hat[,s]) %*% Cj) + (t(Cj) %*% tau.yc.hat[,s] %*% t(yjo) %*% O) + (t(Cj) %*% tau.yc2.hat[,,s] %*% Cj)
        }
        tauyhat[ind,] = tau[ind]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)))) + t(t(Cj) %*% tau.yc.hat)
        yhat[ind, ] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat)         # add yhat
      }}}  
  # individual expected score vector
  q1 = p*(p+1)/2
  xp <- ncol(X)
  m = 1 + xp + q1
  
  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
  }
  
  xloc = seq(0, p*n, p)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(j in 1: n){
    idx = (xloc[j]+1): xloc[(j+1)]
    EIs[1:xp, j] = t(X[idx, ]) %*% Sig.inv %*% (tauyhat[j, ] - tau[j]*X[idx, ] %*% beta)
    A2 = tauy2hat[,,j] - tauyhat[j,] %*% t(X[idx, ] %*% beta) - X[idx, ] %*% beta %*% t(tauyhat[j,]) + (tau[j] * X[idx, ] %*% beta) %*% t(X[idx, ] %*% beta)
    for(k in 1: q1){
      SSdot = Sig.inv %*% dot.Sig[,,k]
      EIs[(xp+k), j] = .5 * (tr(A2 %*% SSdot %*% Sig.inv) - tr(SSdot))
    }
  }
  EIs[(1+xp+q1), ] = (1/nu) + .5 * ka
  
  # Louis (1982) formula
  I.theta = EIs[,1] %*% t(EIs[,1])
  for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
  sd.theta = rep(NA, m)
  V.theta = solve(I.theta)
  sd.theta = sqrt(diag(V.theta))
  return(list(EST=EST, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}
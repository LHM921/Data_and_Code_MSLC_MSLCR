library(mvtnorm)
library(tmvtnorm)
library(mclust)

MVTR.EM = function(Y, X, distr=c('MVN','MVT'), init.par=NULL, tol=1e-5, max.iter=1000, per=10)
{
  require(mvtnorm)
  begin = proc.time()[1]
  n = nrow(Y)
  p = ncol(Y)
  vechS = vech.posi(p)
  y = as.vector(t(Y))

# initial values
  if(length(init.par) != 0){
   beta = init.par$beta
   Sigma = init.par$Sigma
   old.par = c(beta, Sigma[vechS])
   if(distr=='MVT'){
    nu = init.par$nu
    old.par = c(beta, Sigma[vechS], nu)
  }} else{
   beta = solve(t(X)%*%X) %*% t(X) %*% yc
   Sigma = cov(Y)#*(n-1)/n
   old.par = c(beta, Sigma[vechS])
   if(distr=='MVT'){
    nu = runif(1,2,50)
    old.par = c(beta, Sigma[vechS], nu)
  }}
  Sig.inv = solve(Sigma)

# old observed-data log-likelihood
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%Sig.inv%*%cent)
  if(distr=='MVN') w.den = -.5*p*log(2*pi) -.5*log(det(Sigma)) -.5*d
  if(distr=='MVT') w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)
    
    loglik.old = iter.lnL = sum(w.den)
    iter = 0
    if(distr == 'MVN') iter.EST = c(iter, beta, Sigma[vechS])
    if(distr == 'MVT') iter.EST = c(iter, beta, Sigma[vechS], nu)

    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the MNI with ", distr, "distribution.", "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")
    tauycent2 = array(NA, dim=c(p, p, n))
    repeat
    {
      iter = iter + 1 
#E-step
      if(distr=='MVN') tau = rep(1, n)
      if(distr=='MVT') tau = (nu+p)/(nu+d)

# M-step
      TSig = kronecker(diag(n), Sigma)
      TSig.inv = kronecker(diag(n), Sig.inv)
      Tga = rep(tau, each=p)
      beta = solve(t(Tga * X) %*% TSig.inv %*% X) %*% (t(Tga * X) %*% TSig.inv %*% y)
      vec.mu = X %*% beta
      mu = matrix(vec.mu, nrow=p, ncol=n)
      cent = t(Y) - mu
      for(j in 1: n) tauycent2[,,j] = tau[j] *cent[, j] %*% t(cent[, j])
      Sigma =  apply(tauycent2, 1:2, sum)/ n
      Sig.inv = solve(Sigma)
      new.par = c(beta, Sigma[vechS])
      if(distr=='MVT'){
#       nu = optim(par = nu, fn = MVTR.nu.fn, method = "L-BFGS-B", lower = 2+1e-6, upper =Inf, Y=Y, X=X, beta=beta, Sigma=Sigma)$par
       nu = nlminb(start = nu, objective = MVTR.nu.fn, lower = 2+1e-6, upper =Inf, Y=Y, X=X, beta=beta, Sigma=Sigma)$par
       new.par = c(beta, Sigma[vechS], nu)
      }

# new observed-data log-likelihood
    d = diag(t(cent)%*%Sig.inv%*%cent)
    if(distr=='MVN') w.den = -.5*p*log(2*pi) -.5*log(det(Sigma)) -.5*d
    if(distr=='MVT') w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)

    loglik.new = sum(w.den)
    diff = loglik.new - loglik.old
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vechS]))
    if(distr == 'MVT') iter.EST = rbind(iter.EST, c(iter, beta, Sigma[vechS], nu))

    # if(iter%%per == 0){
    # if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, '\t beta=', beta, "\n")
    # if(distr=='MVT') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, '\t nu=', nu, '\t beta=', beta, "\n")
    # }
    if(diff < tol | diff.par < tol | iter > max.iter) break
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
   IM = I.MVTR(para.est, Y, X, distr=distr)
   Yfit = X %*% beta
   e.hat = y - Yfit   
   cat(paste(rep("-", 50), collapse = ""), "\n")
   return(list(run.sec = run.sec, iter=iter, model.inf = model.inf, para.est = para.est, EST = EST, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST, yfit = Yfit, error = e.hat, tau=tau, Del=d))
}


# MVT log-likelihood function for nu
MVTR.nu.fn = function(par, Y, X, beta, Sigma)
{
  nu = par
  n = nrow(Y)
  p = ncol(Y) 
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)
  return(-sum(w.den))
}

# Louis' information matrix
I.MVTR = function(para.est, Y, X, distr=c('MVN','MVT'))
{
  n = nrow(Y)
  p = ncol(Y)
# parameter estimates
  beta = para.est$beta
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MVT') nu = para.est$nu 

# first two moments 
  vec.mu = X %*% beta
  mu = matrix(vec.mu, nrow=p, ncol=n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%Sig.inv%*%cent)
  if(distr=='MVN') tau = tau2 = rep(1, n)
  if(distr=='MVT'){
    tau = (nu+p)/(nu+d)
    tau2 = (2*(nu+p) + (nu+p)^2)/(nu+d)^2
  }

  TSig = kronecker(diag(n), Sigma)
  TSig.inv = kronecker(diag(n), Sig.inv)
  Tga = rep(tau, each=p)
  np = n*p
  TK = matrix(0, np, np)
  cumsum.np=seq(0, np, p)
  sumK2 = 0
  for(j in 1: n){
   np.idx = (cumsum.np[j]+1): cumsum.np[(j+1)] 
   sumK2 = sumK2 + (tau2[j]-tau[j]^2) * t(X[np.idx, ]) %*% Sig.inv %*% cent[, j] %*% t(cent[, j]) %*% Sig.inv %*% X[np.idx, ]
  }
  I.beta = t(Tga * X) %*% TSig.inv %*% X - sumK2  
  sd.beta = sqrt(diag(solve(I.beta))) 
  return(list(I.beta = I.beta, SD=sd.beta, beta.hat=rbind(c(beta), sd.beta)))
}
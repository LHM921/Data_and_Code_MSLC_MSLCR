library(mvtnorm)
library(tmvtnorm)
library(mclust)

MVT.EM = function(Y, distr=c('MVN','MVT'), init.par=NULL, tol=1e-5, max.iter=1000, per=10)
{
  require(mvtnorm)
  begin = proc.time()[1]
  n = nrow(Y)
  p = ncol(Y)
  vechS = vech.posi(p)

# initial values
  if(length(init.par) != 0){
   mu = init.par$mu
   Sigma = init.par$Sigma
   old.par = c(mu, Sigma[vechS])
   if(distr=='MVT'){
    nu = init.par$nu
    old.par = c(mu, Sigma[vechS], nu)
  }} else{
   mu = colMeans(Y)
   Sigma = cov(Y)#*(n-1)/n
   old.par = c(mu, Sigma[vechS])
   if(distr=='MVT'){
    nu = runif(1,2,50)
    old.par = c(mu, Sigma[vechS], nu)
  }}

# old observed-data log-likelihood
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  if(distr=='MVN') w.den = -.5*p*log(2*pi) -.5*log(det(Sigma)) -.5*d
  if(distr=='MVT') w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)
    
    loglik.old = iter.lnL = sum(w.den)
    iter = 0
    if(distr == 'MVN') iter.EST = c(iter, mu, Sigma[vech.posi(p)])
    if(distr == 'MVT') iter.EST = c(iter, mu, Sigma[vech.posi(p)], nu)

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
      mu = colSums(tau*Y)/sum(tau)
      cent = t(Y) - mu
      for(j in 1: n) tauycent2[,,j] = tau[j] *cent[, j] %*% t(cent[, j])
      Sigma =  apply(tauycent2, 1:2, sum)/ n
      new.par = c(mu, Sigma[vechS])
      if(distr=='MVT'){
#       nu = optim(par = nu, fn = MVT.nu.fn, method = "L-BFGS-B", lower = 2+1e-6, upper =Inf, Y=Y, mu=mu, Sigma=Sigma)$par
       nu = nlminb(start = nu, objective = MVT.nu.fn, lower = 2+1e-6, upper =Inf, Y=Y, mu=mu, Sigma=Sigma)$par
       new.par = c(mu, Sigma[vechS], nu)
      }

# new observed-data log-likelihood
    d = diag(t(cent)%*%solve(Sigma)%*%cent)
    if(distr=='MVN') w.den = -.5*p*log(2*pi) -.5*log(det(Sigma)) -.5*d
    if(distr=='MVT') w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)

    loglik.new = sum(w.den)
    diff = loglik.new - loglik.old
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)]))
    if(distr == 'MVT') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)], nu))

    # if (iter%%per == 0){
    # if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, '\t mu=', mu, "\n")
    # if(distr=='MVT') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, '\t nu=', nu, '\t mu=', mu, "\n")
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
    if(diff < 0){ 
      aic = 2 * m - 2 * loglik.old
      bic = m * log(n) - 2 * loglik.old
      model.inf = list(m = m, loglik = loglik.old, aic = aic, bic = bic)
      EST = iter.EST[iter, -1]
      mu.hat = iter.EST[iter, 2:(p+1)]
      signu.hat = iter.EST[iter, -c(1:(p+1))]
      vechS = vech.posi(p); vechS1 = vechS[,c(2,1)]
      Sigma[vechS] = Sigma[vechS1] = signu.hat[1:m1]
      if(distr == 'MVN') para.est = list(mu=mu.hat, Sigma=Sigma)
      if(distr == 'MVT') para.est = list(mu=mu.hat, Sigma=Sigma, nu=signu.hat[(m1+1)])
      cat("iter",iter, ":", "observed-data log-likelihood =", loglik.old, "\t diff =", diff, '\t mu=', mu.hat, '\t Sigma nu=', signu.hat, "\n") 
    } else{
     cat("iter",iter, ":", "observed-data log-likelihood =", loglik.new, "\t diff =", diff, '\t mu=', mu, '\t Sigma=', Sigma[vech.posi(p)], "\n")
     aic = 2 * m - 2 * loglik.new
     bic = m * log(n) - 2 * loglik.new
     model.inf = list(m = m, loglik = loglik.new, aic = aic, bic = bic)
     if(distr == 'MVN'){
       EST = c(mu, Sigma[vech.posi(p)])
       para.est = list(mu=mu, Sigma=Sigma) 
     } 
     if(distr == 'MVT'){
       EST = c(mu, Sigma[vech.posi(p)], nu)
       para.est = list(mu=mu, Sigma=Sigma, nu=nu)
     }
   }
   IM = I.MVT(para.est, Y, distr=distr)
   cat(paste(rep("-", 50), collapse = ""), "\n")
   return(list(run.sec = run.sec, iter=iter, model.inf = model.inf, para.est = para.est, EST = EST, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST))
}


# MVT log-likelihood function for nu
MVT.nu.fn = function(par, Y, mu, Sigma)
{
  nu = par
  n = nrow(Y)
  p = ncol(Y)
  
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  w.den = lgamma((nu+p)/2) - lgamma(nu/2) - .5*p*log(pi*nu)- .5*log(det(Sigma)) - .5*(nu+p)*log(1+d/nu)
  return(-sum(w.den))
}

# Louis' information matrix
I.MVT = function(para.est, Y, distr=c('MVN','MVT'))
{
  n = nrow(Y)
  p = ncol(Y)
# parameter estimates
  mu = para.est$mu
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MVT') nu = para.est$nu 

# first two moments 
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  if(distr=='MVN') tau = tau2 = rep(1, n)
  if(distr=='MVT'){
    tau = (nu+p)/(nu+d)
    tau2 = (2*(nu+p) + (nu+p)^2)/(nu+d)^2
  }
  vartau.Scent2S = array(0, dim=c(p,p,n))
  for(j in 1: n) vartau.Scent2S[,,j] = (tau2[j]-tau[j]^2) * Sig.inv %*% cent[, j] %*% t(cent[, j]) %*% Sig.inv
  I.mu = sum(tau)*Sig.inv - apply(vartau.Scent2S, 1:2, sum)  
  sd.mu = sqrt(diag(solve(I.mu))) 
  return(list(I.mu = I.mu, SD=sd.mu, mu.hat=rbind(mu, sd.mu)))
}
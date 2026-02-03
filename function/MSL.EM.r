MSL.EM = function(Y, distr=c('MVN','MSL'), init.par=NULL, tol=1e-5, max.iter=1000, per=10)
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
   if(distr=='MSL'){
    nu = init.par$nu
    old.par = c(mu, Sigma[vechS], nu)
  }} else{
   mu = colMeans(Y)
   Sigma = cov(Y)#*(n-1)/n
   old.par = c(mu, Sigma[vechS])
   if(distr=='MSL'){
    nu = runif(1,1,10)
    old.par = c(mu, Sigma[vechS], nu)
  }}
# old observed-data log-likelihood
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  if(distr=='MSL'){
     pdf = IGfn(a=((p+nu)/2), b=d/2)
     idx0 = which(d==0)
     pdf[idx0] = 2/(p+nu) 
   } 
  if(distr=='MVN') w.den = -0.5*p*log(2*pi) -0.5*log(det(Sigma)) -0.5*d
  if(distr=='MSL') w.den = log(pdf) + log(nu) - log(2) - 0.5*p*log(2*pi) -0.5*log(det(Sigma))
    
    loglik.old = iter.lnL = sum(w.den)
    iter = 0
    if(distr == 'MVN') iter.EST = c(iter, mu, Sigma[vech.posi(p)])
    if(distr == 'MSL') iter.EST = c(iter, mu, Sigma[vech.posi(p)], nu)

    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the MNI with ", distr, "distribution.", "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")
    tauycent2 = array(NA, dim=c(p, p, n))
    repeat
    {
      iter = iter + 1 
# E-step
      if(distr=='MVN') tau = rep(1, n)
      if(distr=='MSL'){
       tau = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
       idx0 = which(d==0)
       tau[idx0] = (p+nu)/(p+nu+2)
      }
# M-step
      mu = colSums(tau*Y)/sum(tau)
      cent = t(Y) - mu
      for(j in 1: n) tauycent2[,,j] = tau[j] *cent[, j] %*% t(cent[, j])
      Sigma = apply(tauycent2, 1:2, sum)/ n
      new.par = c(mu, Sigma[vechS])
      if(distr=='MSL'){
       nu = optim(par = nu, fn = MSL.nu.fn, method = "L-BFGS-B", lower = 1+1e-6, upper =Inf, Y=Y, mu=mu, Sig=Sigma)$par
       new.par = c(mu, Sigma[vechS], nu)
      }
# new observed-data log-likelihood
    d = diag(t(cent)%*%solve(Sigma)%*%cent)
    if(distr=='MSL'){
       pdf = IGfn(a=(p+nu)/2, b=d/2)
       idx0 = which(d==0)
       pdf[idx0] = 2/(p+nu) 
     } 
    if(distr=='MVN') w.den = -.5*p*log(2*pi) -.5*log(det(Sigma)) -.5*d
    if(distr=='MSL') w.den = log(pdf) + log(nu) - log(2) - 0.5*p*log(2*pi) -0.5*log(det(Sigma))

    loglik.new = sum(w.den)
    diff= loglik.new - loglik.old
    diff.par = sqrt(sum(((old.par-new.par)/old.par)^2))
    iter.lnL = c(iter.lnL, loglik.new)
    if(distr == 'MVN') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)]))
    if(distr == 'MSL') iter.EST = rbind(iter.EST, c(iter, mu, Sigma[vech.posi(p)], nu))

    # if (iter%%per == 0){
    # if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t mu=', mu, "\n")
    # if(distr=='MSL') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t nu=', nu, '\t mu=', mu, "\n")
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
      if(distr == 'MSL') para.est = list(mu=mu.hat, Sigma=Sigma, nu=signu.hat[(m1+1)])
      cat("iter",iter, ":", "observed-data log-likelihood =", loglik.old, "\t diff =", diff, '\t mu=', mu.hat, '\t Sigma nu=', signu.hat, "\n") 
    } else{
     cat("iter",iter, ":", "observed-data log-likelihood =", loglik.new, "\t diff =", diff, '\t mu=', mu, '\t Sigma=', Sigma[vechS], "\n")
     aic = 2 * m - 2 * loglik.new
     bic = m * log(n) - 2 * loglik.new
     model.inf = list(m = m, loglik = loglik.new, aic = aic, bic = bic)
     if(distr == 'MVN'){
       EST = c(mu, Sigma[vechS])
       para.est = list(mu=mu, Sigma=Sigma) 
     } 
     if(distr == 'MSL'){
       EST = c(mu, Sigma[vechS], nu)
       para.est = list(mu=mu, Sigma=Sigma, nu=nu)
     }}
   IM = I.MSL(para.est, Y, distr=distr)
   cat(paste(rep("-", 50), collapse = ""), "\n")
   return(list(run.sec = run.sec, model.inf = model.inf, para.est = para.est, EST = EST, iter.lnL = iter.lnL, IM=IM, iter.EST = iter.EST))
}


# MSL log-likelihood function for nu
MSL.nu.fn = function(par, Y, mu, Sigma)
{
  nu = par
  n = nrow(Y)
  p = ncol(Y)
  
  w.den = matrix(NA, n, 1)
  pdf = w.pdf = rep(1, n)
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  pdf = IGfn(a=(p+nu)/2, b=d/2)
  idx0 = which(d==0)
  pdf[idx0] = 2/(p+nu) 
  w.den = log(pdf) + log(nu) - log(2) - 0.5*p*log(2*pi) -0.5*log(det(Sigma))
  return(-sum(w.den))
}

# Louis' information matrix
I.MSL = function(para.est, Y, distr=c('MVN','MSL'))
{
  n = nrow(Y)
  p = ncol(Y)
# parameter estimates
  mu = para.est$mu
  Sigma = para.est$Sigma
  Sig.inv = solve(Sigma)
  if(distr=='MSL') nu = para.est$nu 

# first two moments 
  cent = t(Y) - mu
  d = diag(t(cent)%*%solve(Sigma)%*%cent)
  if(distr=='MVN') tau = tau2 = rep(1, n)
  if(distr=='MSL'){
    tau = IGfn(a=((p+nu)/2)+1, b=d/2) / IGfn(a=(p+nu)/2, b=d/2) 
    tau2 = IGfn(a=((p+nu)/2)+2, b=d/2) / IGfn(a=(p+nu)/2, b=d/2)
    idx0 = which(d==0)
    tau[idx0] = (p+nu)/(p+nu+2)
    tau2[idx0] = (p+nu)/(p+nu+4)
  }
  vartau.Scent2S = array(0, dim=c(p,p,n))
  for(j in 1: n) vartau.Scent2S[,,j] = (tau2[j]-tau[j]^2) * Sig.inv %*% cent[, j] %*% t(cent[, j]) %*% Sig.inv
  I.mu = sum(tau)*Sig.inv - apply(vartau.Scent2S, 1:2, sum)  
  sd.mu = sqrt(diag(solve(I.mu))) 
  return(list(I.mu = I.mu, SD=sd.mu, mu.hat=rbind(mu, sd.mu)))
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

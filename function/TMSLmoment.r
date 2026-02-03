library(mvtnorm)
library(tmvtnorm)
library(PearsonDS)
library(cubature)
library(invgamma)
 
## generate random sample from truncated multivariate slash distribution
## 'MVN'= multivariate normal ; 'MVT' = multivariate t ; 'MSL' = multivariate slash
rtmsl = function(n, mu, Sigma, nu=NULL, distr=c('MVN','MVT','MSL'), lower=rep(-Inf,length(mu)), upper=rep(Inf, length(mu)))
{
  begin = proc.time()[1]
  distr = distr[1]
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
  Lam = diag(sqrt(diag(Sigma)))
  Lam.inv = diag(1/sqrt(diag(Sigma)))
  R = Lam.inv %*%Sigma%*% Lam.inv 
  }
  Y = matrix(NA, nrow=n, ncol=p)
  Tau = numeric(n)
  iter = 1
  while(iter <= n){
   if(distr=='MVN') tau = 1  
   if(distr=='MVT') tau = rgamma(1, shape=nu/2, rate=nu/2)
   if(distr=='MSL') tau = rbeta(1, shape1=nu/2, shape2=1)
   tau.sq = sqrt(tau)
   if(p==1){
     yi = mu + 1/tau.sq * Sd * rnorm(1, mean=0, sd=1) 
   } else  yi = mu + 1/tau.sq * Lam %*% c(rmvnorm(1, mean=rep(0,p), sigma=R))
   if((sum(yi>=lower)+sum(yi<=upper))==(2*p)){
     Y[iter, ] = yi
     Tau[iter] = tau
     iter = iter + 1
   } else next
  }
  end = proc.time()[1]
  run.sec = end - begin
  return(list(Y=Y, tau=Tau))
}

# incomplete gamma function
IGfn = function(a, b) pgamma(1, shape=a, rate=b) * gamma(a) / b^a

# PDF of multivariate normal, t, and slash distribution
dmsl = function(x, mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, distr=c('MVN','MVT','MSL'))
{
   distr = distr[1]
   p = length(mu)
   if(distr=='MVN'){
     if(p==1){ 
           den = dnorm(x, mu, sd=sqrt(Sigma))
    } else den = dmvnorm(x, mean=mu, sigma=round(Sigma,5), log=F)
   } 
   if(distr=='MVT'){
     if(p==1){ 
           den = dt((x-mu)/sqrt(Sigma), df=nu)/sqrt(Sigma)
    } else den = dmvt(x, delta=mu, sigma=Sigma, df=nu, log=F)
   } 
   if(distr=='MSL'){
     if(p==1){
       den = c(nu / (2*(2*pi)^(p/2)*sqrt(Sigma))) * IGfn(a=p/2+nu/2, b=c(((x-mu)^2/c(Sigma))/2))
       idx = which(x==mu)
       den[idx] = c(nu / ((2*pi)^(p/2)*(p+nu)*sqrt(Sigma)))
     } else{
       den = nu / (2*(2*pi)^(p/2)*sqrt(det(Sigma))) * IGfn(a=p/2+nu/2, b=c(t(x-mu)%*%solve(Sigma)%*%(x-mu))/2) 
       if(sum(x==mu)==p) den = nu / ((2*pi)^(p/2)*(p+nu)*sqrt(det(Sigma)))
   }}
   return(den)
}

# CDF of multivariate normal, t, and slash distribution
pmsl = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, distr=c('MVN','MVT','MSL'))
{
   require(mvtnorm)
   GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
   distr = distr[1]
   p = length(mu)
   if(distr=='MVN'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd)
    } else cdf = pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
   }
   if(distr=='MVT'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pt((upper-mu)/Sd, df=nu) - pt((lower-mu)/Sd, df=nu)
    } else cdf = pmvt(lower=lower, upper=upper, delta=mu, df=nu, sigma=Sigma)[1]
   }
   if(distr=='MSL'){
    if(p==1) cdf = integrate(dmsl, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$value
    else cdf = cubintegrate(dmsl, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
#    cdf = adaptIntegrate(dmsl, lowerLimit=lower, upperLimit=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
   }  
   return(cdf)  
}

# PDF of truncated multivariate normal, t, and slash distribution
dtmsl = function(x, mu, Sigma, nu=NULL, distr=c('MVN','MVT','MSL'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
   distr = distr[1]
   p = length(mu)
   if((sum(x>=a.low)+sum(x<=a.upp))!=(2*p)){
    Tden=0
   } else {
   if(distr=='MVN'){
       cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
       Tden = cdf.inv * dmsl(x, mu, Sigma, nu, distr='MVN') 
   }
   if(distr=='MVT'){
       cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
       Tden = cdf.inv * dmsl(x, mu, Sigma, nu, distr='MVT') 
   }
   if(distr=='MSL'){
     cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MSL')
     Tden = cdf.inv * dmsl(x, mu, Sigma, nu, distr='MSL') 
   }} 
   return(Tden)
}

# CDF of truncated multivariate normal, t, and slash distribution
ptmsl = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, distr=c('MVN','MVT','MSL'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
   distr = distr[1]
   p = length(mu)
   for(i in 1: p){
     if(lower[i]<a.low[i]) lower[i]=a.low[i]
     if(upper[i]>a.upp[i]) upper[i]=a.upp[i]
   }
   if(distr=='MVN'){
       cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
       if(p==1) Tcdf = cdf.inv * (pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd))
       else Tcdf = cdf.inv * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
   }
   if(distr=='MVT'){
       cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
       if(p==1) Tcdf = cdf.inv * (pt((upper-mu)/sqrt(Sigma), df=nu) - pt((lower-mu)/sqrt(Sigma), df=nu))
       else Tcdf = cdf.inv * pmvt(lower=lower, upper=upper, delta=mu, df=nu, sigma=Sigma)[1]
   }
   if(distr=='MSL'){
     cdf.inv = 1/pmsl(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MSL')
     Tcdf = cdf.inv * cubintegrate(dmsl, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
   }   
   return(Tcdf)
}

gaqh.SL = function(x, dd, a.ast2, mu2.1, R22.1)
{
  di = length(mu2.1)
  if(di==1) del.i = (x-mu2.1)^2/c(R22.1)
  else del.i = t(x-mu2.1)%*%solve(R22.1)%*%(x-mu2.1)
  tag = IGfn(a=dd, b=c(a.ast2+del.i)/2)
  return(tag)
}

# first two moments of truncated multivariate normal, t, and slash distribution
TMSL.moment = function(mu, Sigma, nu=NULL, distr=c('MVN','MVT','MSL'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{ 
  distr=distr[1]
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
   Lambda = diag(sqrt(diag(Sigma)))
   Lambda.inv = diag(1/sqrt(diag(Sigma)))
   R = Lambda.inv %*% Sigma %*% Lambda.inv
   if(det(R)<=0) stop("The R matrix must be inversible!")
  }
  if(distr=='MVN'){
    a.low = ifelse(a.low==-Inf,rep(-5e1+mu,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(5e1+mu,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     EX = Sd*(dnorm(a)-dnorm(b))/ (pnorm(b)-pnorm(a))
     EY = mu + EX
     CovY = CovX = Sd^2*(1+(a*(dnorm(a))-b*dnorm(b))/(pnorm(b)-pnorm(a))-((dnorm(a)-dnorm(b))/(pnorm(b)-pnorm(a)))^2)
     EXX = CovX + EX^2
     EYY = CovY + EY^2
    } else{
    a = c(Lambda.inv %*% (a.low - mu))
    b = c(Lambda.inv %*% (a.upp - mu))
    al0 = pmvnorm(lower = a, upper = b, sigma = round(R,4), algorithm = GB)[1]
  ### pdf & cdf
    f1a = dnorm(a)
    f1b = dnorm(b)
    f2 = matrix(NA, p, p)
    G1a = G1b = rep(NA, p)
    G2 = matrix(NA, p, p)
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
      G1a[r] = pmvnorm(lower = low, upper = upp, sigma = round(S1,4), algorithm = GB)[1]
      mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
      G1b[r] = pmvnorm(lower = low, upper = upp, sigma = round(S1,4), algorithm = GB)[1]
    }
    qa = f1a*G1a; qb = f1b*G1b
    EX = c(R %*% (qa-qb)) / al0

    H = matrix(0,p,p)
    for(r in 1:(p-1)){
      for(s in (r+1):p){
        rs = c(r,s)
        pdf.aa = dmvnorm(t(c(a[r],a[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.ab = dmvnorm(t(c(a[r],b[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.ba = dmvnorm(t(c(b[r],a[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.bb = dmvnorm(t(c(b[r],b[s])),sigma=round(R[rs,rs],5), log =F)
        if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
        if(p>2){
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          mu.aa = c(tmp%*%c(a[r],a[s]))
          mu.ab = c(tmp%*%c(a[r],b[s]))
          mu.ba = c(tmp%*%c(b[r],a[s]))
          mu.bb = c(tmp%*%c(b[r],b[s]))
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          cdf.aa = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.aa, sigma = round(R21,5), algorithm = GB)[1]
          cdf.ab = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ab, sigma = round(R21,5), algorithm = GB)[1]
          cdf.ba = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ba, sigma = round(R21,5), algorithm = GB)[1]
          cdf.bb = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.bb, sigma = round(R21,5), algorithm = GB)[1]
        }
        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
      }}
    D = matrix(0,p,p)
    diag(D) = a * qa - b * qb - diag(R%*%H)
    EXX = R + R %*% (H + D) %*% R / al0
  }}
  if(distr=='MVT'){
    if(nu <= 2) stop("The first moment exists only when the degree of freedom is larger than 2!")
    if(nu <= 4) stop("The theoretical second moment exists only when the degrees of freedom is larger than 4!")
    a.low = ifelse(a.low==-Inf,rep(-1e12,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(1e12,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     ka = gamma((nu+1)/2)/((pt(b, df=nu)-pt(a, df=nu))*gamma(nu/2)*sqrt(nu*pi))
     EX = (ka*nu)/(nu-1) * ((1+a^2/nu)^(-(nu-1)/2) - (1+b^2/nu)^(-(nu-1)/2))
     tau.sq = sqrt((nu-2)/nu)
     EXX = (nu*(nu-1))/(nu-2) * ((pt(b*tau.sq, df=nu-2)-pt(a*tau.sq, df=nu-2)) / (pt(b, df=nu)-pt(a, df=nu))) - nu
     EY = mu + Sd * EX
     EYY = Sigma * EXX + mu^2 + 2 * mu * Sd * EX
     CovX = EXX - EX^2
     CovY = EYY - EY^2 
    } else{
    aL.ast = c(Lambda.inv %*% (a.low - mu))
    aU.ast = c(Lambda.inv %*% (a.upp - mu))
    tau.t = pmvt(lower = aL.ast, upper = aU.ast, sigma = R, df = round(nu), algorithm = GB)[1]
    la1 = (nu-2)/nu; la2 = (nu-4)/nu
    da = (nu-1)/(nu+aL.ast^2); db = (nu-1)/(nu+aU.ast^2)
    f1a = sqrt(la1)*dt(sqrt(la1)*aL.ast,df=nu-2)
    f1b = sqrt(la1)*dt(sqrt(la1)*aU.ast,df=nu-2)
    G1a = G1b = rep(NA, p)
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * aL.ast[r]; low = aL.ast[-r]-mua; upp = aU.ast[-r]-mua
      G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1),pmvt(lower = low, upper = upp, sigma = S1/da[r], df = round(nu-1), algorithm = GB)[1])
      mub = temp * aU.ast[r]; low = aL.ast[-r]-mub; upp = aU.ast[-r]-mub
      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1),pmvt(lower = low, upper = upp, sigma = S1/db[r], df = round(nu-1), algorithm = GB)[1])
    }
    qa = f1a*G1a; qb = f1b*G1b
    EX = c(R %*% (qa-qb)) / tau.t / la1
    H = matrix(0,p,p)
    for(r in 1:(p-1)){
     for(s in (r+1):p){
       rs = c(r,s)
       pdf.aa = dmvt(c(aL.ast[r],aL.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.ab = dmvt(c(aL.ast[r],aU.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.ba = dmvt(c(aU.ast[r],aL.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.bb = dmvt(c(aU.ast[r],aU.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
       if(p>2){
         tmp = R[-rs,rs]%*%solve(R[rs,rs])
         mu.aa = c(tmp%*%c(aL.ast[r],aL.ast[s]))
         mu.ab = c(tmp%*%c(aL.ast[r],aU.ast[s]))
         mu.ba = c(tmp%*%c(aU.ast[r],aL.ast[s]))
         mu.bb = c(tmp%*%c(aU.ast[r],aU.ast[s]))
         daa = (nu-2)/(nu+(aL.ast[r]^2-2*R[r,s]*aL.ast[r]*aL.ast[s]+aL.ast[s]^2)/(1-R[r,s]^2))
         dab = (nu-2)/(nu+(aL.ast[r]^2-2*R[r,s]*aL.ast[r]*aU.ast[s]+aU.ast[s]^2)/(1-R[r,s]^2))
         dba = (nu-2)/(nu+(aU.ast[r]^2-2*R[r,s]*aU.ast[r]*aL.ast[s]+aL.ast[s]^2)/(1-R[r,s]^2))
         dbb = (nu-2)/(nu+(aU.ast[r]^2-2*R[r,s]*aU.ast[r]*aU.ast[s]+aU.ast[s]^2)/(1-R[r,s]^2))
         R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
         cdf.aa = ifelse(p==3,pt((aU.ast[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((aL.ast[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.aa, upper = aU.ast[-rs]-mu.aa, sigma = R21/daa, df=round(nu-2), algorithm = GB)[1])
         cdf.ab = ifelse(p==3,pt((aU.ast[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((aL.ast[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.ab, upper = aU.ast[-rs]-mu.ab, sigma = R21/dab, df=round(nu-2), algorithm = GB)[1])
         cdf.ba = ifelse(p==3,pt((aU.ast[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((aL.ast[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.ba, upper = aU.ast[-rs]-mu.ba, sigma = R21/dba, df=round(nu-2), algorithm = GB)[1])
         cdf.bb = ifelse(p==3,pt((aU.ast[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((aL.ast[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.bb, upper = aU.ast[-rs]-mu.bb, sigma = R21/dbb, df=round(nu-2), algorithm = GB)[1])
       }
       H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
     }}
     H = H / la2
     D = matrix(0,p,p)
     diag(D) = aL.ast * qa - aU.ast * qb - diag(R%*%H)
     al1 = pmvt(lower = aL.ast, upper = aU.ast, sigma = R/la1, df=round(nu-2), algorithm = GB)[1]
     EXX = (al1 * R + R %*% (H + D) %*% R) / tau.t / la1
  }}
  if(distr=='MSL'){
    if(nu <= 1) stop("The theoretical first moment exists only when the degrees of freedom is larger than 1!")
    if(nu <= 2) stop("The theoretical second moment exists only when the degrees of freedom is larger than 2!")
    mu0 = ifelse((a.low==-Inf & a.upp==Inf), rep(1e-3, p), rep(0, p))
    a.low = ifelse(a.low==-Inf,rep(-5e1+mu,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(5e1+mu,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     Fa = pmsl(lower=-Inf, upper=a, mu=0, Sigma=1, nu=nu, distr='MSL')
     Fb = pmsl(lower=-Inf, upper=b, mu=0, Sigma=1, nu=nu, distr='MSL')
     Bnu1 = nu /(2*(sqrt(2*pi)*(Fb-Fa)))    
     EX = Bnu1 * (IGfn((nu-1)/2, .5*a^2) - IGfn((nu-1)/2, .5*b^2)) 
     Fa1 = pmsl(lower=-Inf, upper=a, mu=0, Sigma=1, nu=nu-2, distr='MSL')
     Fb1 = pmsl(lower=-Inf, upper=b, mu=0, Sigma=1, nu=nu-2, distr='MSL')    
     EXX = ((Fb1-Fa1)/(Fb-Fa))*(nu/(nu-2)) + Bnu1 * (a*IGfn((nu-1)/2, .5*a^2) - b*IGfn((nu-1)/2, .5*b^2)) 
     EY = mu + Sd * EX
     EYY = Sigma * EXX + mu^2 + 2 * mu * Sd * EX
     CovX = EXX - EX^2
     CovY = EYY - EY^2 
    } else{
    aL.ast = c(Lambda.inv %*% (a.low - mu))
    aU.ast = c(Lambda.inv %*% (a.upp - mu))
    tau.s = pmsl(lower = aL.ast, upper = aU.ast, mu=mu0, Sigma=R, nu=nu, distr='MSL')
    la1 = nu/(2*(2*pi)^(p/2))
    d1 = (p+nu-2)/2
    G1a = G1b = rep(NA, p)
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * aL.ast[r]
      G1a[r] = ifelse(p==2, integrate(gaqh.SL, lower=aL.ast[-r], upper=aU.ast[-r], dd=d1, a.ast2=aL.ast[r]^2, mu2.1=mua, R22.1=S1)$value/sqrt(c(S1)), cubintegrate(gaqh.SL, lower=aL.ast[-r], upper=aU.ast[-r], dd=d1, a.ast2=aL.ast[r]^2, mu2.1=mua, R22.1=S1)$integral/sqrt(det(S1)))
      mub = temp * aU.ast[r]
      G1b[r] = ifelse(p==2, integrate(gaqh.SL, lower=aL.ast[-r], upper=aU.ast[-r], dd=d1, a.ast2=aU.ast[r]^2, mu2.1=mub, R22.1=S1)$value/sqrt(c(S1)), cubintegrate(gaqh.SL, lower=aL.ast[-r], upper=aU.ast[-r], dd=d1, a.ast2=aU.ast[r]^2, mu2.1=mub, R22.1=S1)$integral/sqrt(det(S1)))
    }
    qa = la1*G1a; qb = la1*G1b
    EX = c(R %*% (qa-qb)) / tau.s

    H = matrix(0,p,p)
    d2 = (p+nu-4)/2
    la3 = nu/(2*2*pi*sqrt(det(R)))
    for(r in 1:(p-1)){
     for(s in (r+1):p){
       la2 = nu/(sqrt(1-R[r,s]^2)*((2*pi)^(p/2)))
       rs = c(r,s)
       del.aa = c(t(c(aL.ast[r],aL.ast[s])) %*% solve(R[rs,rs]) %*% c(aL.ast[r],aL.ast[s]))
       del.ab = c(t(c(aL.ast[r],aU.ast[s])) %*% solve(R[rs,rs]) %*% c(aL.ast[r],aU.ast[s]))
       del.ba = c(t(c(aU.ast[r],aL.ast[s])) %*% solve(R[rs,rs]) %*% c(aU.ast[r],aL.ast[s]))
       del.bb = c(t(c(aU.ast[r],aU.ast[s])) %*% solve(R[rs,rs]) %*% c(aU.ast[r],aU.ast[s]))
       if(p==2){
         cdf.aa = cdf.ab = cdf.ba= cdf.bb = 1
         pdf.aa = IGfn(a=nu/2-1, b=del.aa/2) * la3
         pdf.ab = IGfn(a=nu/2-1, b=del.ab/2) * la3
         pdf.ba = IGfn(a=nu/2-1, b=del.ba/2) * la3
         pdf.bb = IGfn(a=nu/2-1, b=del.bb/2) * la3
       }
       if(p>2){
         pdf.aa = pdf.ab = pdf.ba = pdf.bb = 1
         tmp = R[-rs,rs]%*%solve(R[rs,rs])
         mu.aa = c(tmp%*%c(aL.ast[r],aL.ast[s]))
         mu.ab = c(tmp%*%c(aL.ast[r],aU.ast[s]))
         mu.ba = c(tmp%*%c(aU.ast[r],aL.ast[s]))
         mu.bb = c(tmp%*%c(aU.ast[r],aU.ast[s]))
         R21 = R[-rs,-rs] - R[-rs,rs] %*% solve(R[rs,rs]) %*% R[rs,-rs]
         cdf.aa = ifelse(p==3, integrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.aa, mu2.1=mu.aa, R22.1=R21)$value, cubintegrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.aa, mu2.1=mu.aa, R22.1=R21)$integral)*la2/sqrt(det(R21))
         cdf.ab = ifelse(p==3, integrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.ab, mu2.1=mu.ab, R22.1=R21)$value, cubintegrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.ab, mu2.1=mu.ab, R22.1=R21)$integral)*la2/sqrt(det(R21))
         cdf.ba = ifelse(p==3, integrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.ba, mu2.1=mu.ba, R22.1=R21)$value, cubintegrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.ba, mu2.1=mu.ba, R22.1=R21)$integral)*la2/sqrt(det(R21))
         cdf.bb = ifelse(p==3, integrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.bb, mu2.1=mu.bb, R22.1=R21)$value, cubintegrate(gaqh.SL, lower=aL.ast[-rs], upper=aU.ast[-rs], dd=d2, a.ast2=del.bb, mu2.1=mu.bb, R22.1=R21)$integral)*la2/sqrt(det(R21))
       }
       H[r,s] = H[s,r] = (pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb)
     }}
     K = matrix(0,p,p)
     diag(K) = aL.ast * qa - aU.ast * qb - diag(R%*%H)
     al1 = pmsl(lower = aL.ast, upper = aU.ast, mu=mu0, Sigma=R, nu=nu-2, distr='MSL')
     EXX = (nu/(nu-2)*al1*R + R %*% (H + K) %*% R) / tau.s
  }}
  if(p >= 2){
   EY = c(mu + Lambda %*% EX)
   EYY = mu%*%t(mu) + Lambda%*%EX%*%t(mu) + mu%*%t(EX)%*%Lambda + Lambda%*%EXX%*%Lambda
   CovY = EYY-(EY)%*%t(EY)
  }
  return(list(EY=EY, EYY=EYY, CovY=CovY, EX=EX, EXX=EXX))
}

#
# Functions for estimating the stock-recruitment relationship
#

SR.funcs <- function(p, rec, ssb, gamma1=0.1, d=1,type="L2", SR="HS", HSmod="HS", AutoCor=FALSE, rho.fix=FALSE,like="exact", out="obj"){
  require(VGAM)
  a <- exp(p[1])
  if (SR=="HS") b <- (max(ssb)-min(ssb))/(1+exp(-p[2]))+min(ssb) else b <- exp(p[2])
  if (AutoCor) {if (rho.fix) rho <- 0 else rho <- 2/(1+exp(-p[3]))-1} else rho <- 0
   
  n <- length(rec)
  
  if (SR=="BH") pred <- a*ssb/(1+(1/b)*ssb)
  if (SR=="RI") pred <- a*ssb*exp(-ssb/b)
  if (SR=="HS") if (HSmod=="Mesnil") pred <- a*(ssb+sqrt(b^2+(gamma1^2)/4)-sqrt((ssb-b)^2+(gamma1^2)/4)) else pred <- ifelse(ssb > b, a*b, a*ssb)  
  
  if(type=="L1") {
    if (AutoCor){      # conditional likelihood
      log.pred.ar <- numeric(n)
      log.pred.ar[1] <- log(pred[1])
      resid.ar <- log(rec[1])-log.pred.ar[1]
      sigma <- ifelse(like=="cond",0,sqrt(1-rho^2)*abs(resid.ar))
      for (i in 2:n){
       log.pred.ar[i] <- log(pred[i])+rho*resid.ar
        resid.ar <- log(rec[i])-log.pred.ar[i]
        sigma <- sigma + abs(resid.ar)
      }
      sigma <- (1/(n-(like=="cond")))*sigma
      obj <- -(sum(dlaplace(log(rec[2:n]), log.pred.ar[2:n], sigma, log=TRUE))+(like!="cond")*dlaplace(log(rec[1]), log.pred.ar[1], sigma/sqrt(1-rho^2), log=TRUE))
    } else{
      sigma <- mean(abs(log(rec)-log(pred)))      # Laplace distribution
      obj <- -sum(dlaplace(log(rec), log(pred), sigma, log=TRUE))
    }
  } else {
    if (AutoCor){      # conditional likelihood
      log.pred.ar <- numeric(n)
      log.pred.ar[1] <- log(pred[1])      
      resid.ar <- log(rec[1])-log.pred.ar[1]
      sigma <- ifelse(like=="cond",0,(1-rho^2)*resid.ar^2)     
      for (i in 2:n){
       log.pred.ar[i] <- log(pred[i])+rho*resid.ar
        resid.ar <- log(rec[i])-log.pred.ar[i]
        sigma <- sigma + resid.ar^2
      }
      sigma <- (1/(n-(like=="cond")))*sigma
      obj <- -(sum(dnorm(log(rec[2:n]), log.pred.ar[2:n], sqrt(sigma), log=TRUE))+(like!="cond")*dnorm(log(rec[1]), log.pred.ar[1], sqrt(sigma/(1-rho^2)), log=TRUE))
    } else{
      sigma <- mean((log(rec)-log(pred))^2)
      obj <- -sum(dnorm(log(rec), log(pred), sqrt(sigma), log=TRUE))
    }
  }
  
  sigma <- ifelse(type=="L2",sqrt(sigma),sigma*sqrt(2))
  
  if (out=="obj") res <- obj
  if (out=="par") {res <- c(a,b,sigma,rho); names(res) <- c("a","b","sigma","rho")}
  
  res
}

pred.SR <- function(res,ssb){
  SR <- res$SR
  SR.par <- res$pars
  HSmod <- res$HSmod
  gamma1 <- res$gamma
  
  a <- SR.par[1]
  b <- SR.par[2]
  
  if (SR=="HS"){
      if (HSmod=="Mesnil") pred <- a*(ssb+sqrt(b^2+(gamma1^2)/4)-sqrt((ssb-b)^2+(gamma1^2)/4)) else pred <- ifelse(ssb > b, a*b, a*ssb)
  }
  if (SR=="BH") pred <- a*ssb/(1+ssb/b)
  if (SR=="RI") pred <- a*ssb*exp(-ssb/b)
  return(pred)
}

estSR <- function(SRdata, gamma1=0.1, type="L1",SR="HS", HSmod="HS",
  p0=NULL,like="exact",
  Length=20,
  rho.range=seq(-0.99,0.99,by=10),
  rho.fix=FALSE,
  AutoCor=FALSE,
  zero.cut=TRUE
  ){
  R <- SRdata$R
  SSB <- SRdata$SSB
  
  if (zero.cut){
    ok.dat <- which(R!=0 & SSB!=0)
    R <- R[ok.dat]
    SSB <- SSB[ok.dat]
  }

  a.R <- range(R/SSB)          
  a.range <- seq(a.R[1],a.R[2],len=Length)      

  b.R <- range(SSB)      
  b.range <- seq(b.R[1],b.R[2],len=Length)      
  
  if (SR=="HS") b.range <- pmin(pmax(1/((max(SSB)-min(SSB))/(b.range-min(SSB))-1),0.0000001),10^50)
  if (is.null(p0)){
    if (AutoCor & !(rho.fix)){
      rho.range <- pmin(pmax(1/(2/(1+rho.range)-1),0.000001),1000000)
      parms <- expand.grid(a.range,b.range,rho.range)
    } else parms <- expand.grid(a.range,b.range)
      res <- matrix(apply(log(parms),1,SR.funcs,rec=as.numeric(R),ssb=as.numeric(SSB),gamma1=gamma1, type=type, SR=SR, HSmod=HSmod,AutoCor=AutoCor,rho.fix=rho.fix,like=like,out="obj"),nrow=length(a.range),ncol=length(b.range))
      p <- as.numeric(parms[which.min(res),])
  } else p <- p0
  
  res <- optim(log(p), SR.funcs, rec=as.numeric(R), ssb=as.numeric(SSB), gamma1=gamma1, type=type, SR=SR, HSmod=HSmod,AutoCor=AutoCor,rho.fix=rho.fix,like=like,out="obj")
  res <- optim(res$par, SR.funcs, rec=as.numeric(R), ssb=as.numeric(SSB), gamma1=gamma1,type=type, SR=SR, HSmod=HSmod, AutoCor=AutoCor, rho.fix=rho.fix,like=like,out="obj", method="BFGS")
  
  SR.par <- SR.funcs(p=res$par, rec=as.numeric(R), ssb=as.numeric(SSB), gamma1=gamma1, type=type, SR=SR, HSmod=HSmod, AutoCor=AutoCor, rho.fix=rho.fix, like=like, out="par")
  
  res$like <- like
  res$SR <- SR
  res$HSmod <- ifelse(SR=="HS", HSmod, NA)
  res$gamma <- ifelse(HSmod=="Mesnil" & SR=="HS", gamma1, NA)
  res$type <- type
  res$pars <- SR.par
  res$AutoCor <- AutoCor
 
  pred.R <- pred.SR(res,as.numeric(SSB))
  
  res$pred <- pred.R
  res$resid <- log(R)-log(pred.R)
  
  res$logLik <- -res$value
  np <- length(SR.par)-1*rho.fix-(1-1*(AutoCor))
  n <- length(pred.R)
  res$AIC <- 2*res$value+2*np
  res$AICc <- res$AIC+2*(np^2+np)/(n-np-1)
  
  return(res)
}

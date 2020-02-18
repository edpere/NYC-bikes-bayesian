#CPO and LPML
LPML_fun<- function(sp.model){
  N=length(sp.model$Y)
  mean.draws.data=(sp.model$p.beta.recover.samples)%*%t(sp.model$X)+t(sp.model$p.w.recover.samples)
  cpo=seq(1,N)
  for (i in 1:N){
    cpo[i]=1/mean(1/dnorm(sp.model$Y[i],mean.draws.data[,i], sp.model$p.theta.recover.samples[,2])) #1=sigma 2=tau 3=phi
  }
  LPML=sum(log(cpo))
  return(LPML)
}

WAIC <- function(sp.model){
  N=length(sp.model$Y) # n.of stations
  sample.Y=(sp.model$p.beta.recover.samples)%*%t(sp.model$X)+t(sp.model$p.w.recover.samples)
  # Y = beta*X + w (using samples from the model for beta and w)
  lppd=seq(1,N)
  p.waic=seq(1,N)
  for (i in 1:N){
    lppd[i]=log(mean(dnorm(sp.model$Y[i],sample.Y[,i], sp.model$p.theta.recover.samples[,2])))
    #log of sample mean of likelihood (Y|sigma.sq,tau.sq,phi sampled from a Normal(beta*X+w,tau.sq))
    p.waic[i]=var(log(dnorm(sp.model$Y[i],sample.Y[,i], sp.model$p.theta.recover.samples[,2])))
    #sample variance of log likelihood
  }
  #sum over all the stations
  LPPD=sum(lppd)
  PWAIC=sum(p.waic)
  return(-2*(LPPD-PWAIC))
}

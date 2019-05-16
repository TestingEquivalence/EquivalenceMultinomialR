source("distance.R")
source("asymptotic_test.R")

randomPoint<-function(d){
  x=runif(d,0,1)
  x=x/sum(x)
  return(x)
}
  
randomExteriorPoint<-function(eps, q,b){
  repeat{
  x = randomPoint(length(q))
  t=smooth_tvd(x,q,b)
  if (t>eps) return(x)
  }
}

linComb<-function(x,y,a){
  return((1-a)*x+a*y) 
}

linearBoundaryPoint<-function(p,q,eps,b){
  aim<-function(a){
    lc=linComb(p,q,a)
    tvd=smooth_tvd(lc,q,b)
    return(tvd-eps)
  }
  
  aMin=uniroot(aim, c(0,1))
  return(linComb(p,q,aMin$root))
}

randomLinBoundaryPoint<-function(q,eps,b){
  p=randomExteriorPoint(eps,q,b)
  lc=linearBoundaryPoint(p,q,eps,b)
  return(lc)
}

protoBstTest<-function(p,q,n,b,eps,exteriorPoints,nSimulation){
  #calculate test statistic
  t=smooth_tvd(p,q,b)
  rp=p
  
  #estimate closest boundary point
  df<-function(bp){
    smooth_tvd(bp,p,b)
  }
  
  if (t<eps){
    bps=lapply(exteriorPoints,linearBoundaryPoint,q=q,eps=eps,b=b)
    dst=lapply(bps, df)
    pos=which.min(dst)
    rp=bps[pos]
    rp=unlist(rp)
  }
  
  #simulate bootstrap sample
  i=c(1:nSimulation)
  f<-function(k){
    v=rmultinom(n=1,size=n,prob=rp)
    v=v/sum(v)
    return(smooth_tvd(v,q,b))
  }
  sample=lapply(i,f)
  
  #bootstrap test
  pValue=sum(sample<t)/nSimulation
  return(pValue)
}

#' The bootstrap test is based on the re-sampling method called bootstrap. 
#' The bootstrap test is more precise and reliable than the asymptotic test. 
#' However, it should be used carefully because the test is approximate 
#' and may be anti conservative. 
#' In order to obtain a conservative test reducing of alpha
#' (usually halving) or slight shrinkage of the tolerance parameter epsilon
#' may be appropriate. We prefer the slight shrinkage of the tolerance parameter 
#' because it is more effective and the significance level remains unchanged.
#' 
#' \code{bootstrap_test} bootstrap test for equivalence between
#' the observed counting frequencies and fully specified multinomial distribution
#' @param p observed counting frequencies
#' @param q probability vector of the fully specified multinomial distribution
#' @param n number of observations
#' @param b smoothing parameter for the total variation distance
#' @param alpha significance level
#' @param nSimulation number of bootstrap samples, default 10000
#' @param nExteriorPoints number of random directions to search for a boundary point,
#' default is length(q)*200
#' @return tests returns the minimum tolerance parameter epsilon,
#' for which the equivalence between the theoretical and observed distribution can be shown
#' at the significance level alpha


bootstrap_test<-function(p,q,n,b, alpha, 
                         nSimulation=10000, 
                         nExteriorPoints=0){
  
  #number of search directions and seed
  if (nExteriorPoints==0) nExteriorPoints=length(q)*200
  set.seed(10071977)
  
  #find start value for min eps
  #use for this purpose the asymptotic test with 
  #small safety margin
  eps=asymptotic_test(p,q,n,b,alpha)*1.1
  
  #calculate exterior points
  f<-function(x){
    randomExteriorPoint(eps,q,b)
  }
  
  i=c(1:nExteriorPoints)
  exteriorPoints=lapply(i, f)
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(p,q,n,b,eps=x,exteriorPoints, nSimulation)
    pval-alpha
  }
  
  #check boundary values
  #check lower bound
  lb=ff(0)
  if (lb<0) return(0)
  
  #check upper bound
  ub=ff(eps)
  if (ub>0) return(NA)
  
  res=uniroot(ff,c(0,eps))
  return(res$root)
}
  



source("distance.R")

asympt_stdev<-function(p,q,b){
  vec = derivative_smooth_tvd(p,q,b)
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}



#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' Therefore the asymptotic test need some sufficiently large number of the observations.
#' It should be used carefully because the test is approximate 
#' and may be anti conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate. 
#' \code{asymptotic_test} asymptitic test for equivalence between
#' the observed counting frequencies and fully specified multinomial distribution
#' @param p observed counting frequencies
#' @param q probability vector of the fully specified multinomial distribution
#' @param n number of observations
#' @param b smoothing parameter for the total variation distance
#' @param alpha significance level
#' @return tests returns the minimum tolerance parameter epsilon,
#' for which the equivalence between the theoretical and observed distribution can be shown
#' at the significance level alpha

 
asymptotic_test<-function(p,q,n,b,alpha){
  
  vol=asympt_stdev(p,q,b)
  vol=vol / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t = smooth_tvd(p,q,b)
  min_eps = t + qt * vol
}
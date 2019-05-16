smooth_abs<-function(x,b){
  v=x*x+b*b
  v=sqrt(v)
  return(v)
}

derivative_smooth_abs<-function(x,b){
  v=smooth_abs(x,b)
  return(x/v)
}

smooth_tvd<-function(x,y,b){
  v=smooth_abs(x-y,b)
  s=sum(v)
  return(s/2)
}

derivative_smooth_tvd<-function(x,y,b){
  return(derivative_smooth_abs(x-y,b)/2)
}
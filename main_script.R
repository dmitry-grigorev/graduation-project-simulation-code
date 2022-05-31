library(matrixStats)
library(ggplot2)

K<-function(x){pmax(0, 0.75*(1-x^2))}
K.vec<-Vectorize(K)

g<-Vectorize(function(y){y})

f<-function(x) {if (x<0) 
  x - 0 
  else 
    x + 0}
f.vec<-Vectorize(f)

KDE<-function(X, x, h, is.x.boundary = TRUE)
{
  if(is.x.boundary)
  {
    mean((K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h )))/h
  }
  else
  {
    mean((K.vec( (X-x)/h )))/h
  }
}

KE_mean<-function(X,Y, x, h, is.x.boundary = TRUE){
  if(is.x.boundary)
  {
    sum(Y*(K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h )))/sum((K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h )))
  }
  else
  {
    sum(Y*K.vec( (X-x)/h ))/sum((K.vec( (X-x)/h )))
  }
}

KE_mean.derivx<-function(X,Y,x,y,h,eps, is.x.boundary = TRUE)
{
  (KE_mean(X, Y, x+eps, h, is.x.boundary) - KE_mean(X, Y, x, h, is.x.boundary))/eps
}

KE_median<-function(X,Y,x,h, is.x.boundary = TRUE)
{
  #isotone package function weighted.median finds median with given weights
  if(is.x.boundary)
  {
  weightedMedian(Y, (K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h )))
  }
  else
  {
  weightedMedian(Y, K.vec((X-x)/h))
  }
}

KE_cond.pmf<-function(X,Y,x,y,h, is.x.boundary = TRUE)
{
  hy<-1.06*sd(Y)*length(Y)^(-1/3)
  if(is.x.boundary)
  {
    sum((K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h ))*K.vec((Y-y)/hy))/sum((K.vec( (g(X)-x)/h ) + K.vec( (g(X)+x)/h )))/hy
  }
  else
  {
  sum(K.vec((X-x)/h)*K.vec((Y-y)/hy))/sum(K.vec((X-x)/h))/hy
  }
}

KE_cond.cdf<-function(X,Y,x,y,h, is.x.boundary = TRUE)
{
  KE_mean(X, as.numeric(Y<y), x, h, is.x.boundary)
}

KE_cond.cdf.derivx<-function(X,Y,x,y,h,eps, is.x.boundary = TRUE)
{
  (KE_cond.cdf(X, Y, x+eps, y, h, is.x.boundary) - KE_cond.cdf(X, Y, x, y, h, is.x.boundary))/eps
}

KE_var<-function(X, Y, x, y, h, N)
{
  N/(KE_cond.pmf(X,Y,x,y,h, TRUE))^2/KDE(X, x, h, TRUE)
}

stat.correction<-function(X,Y,x,y, h, eps, M)
{
  -2*M*KE_cond.cdf.derivx(X,Y,x,y,h,eps, TRUE)/KE_cond.pmf(X, Y, x, y, h, TRUE)
}

stat.correction.mean<-function(X,Y,x,y,h,eps,M)
{
  2*M*KE_mean.derivx(X, Y, x, y, h, eps, TRUE)
}
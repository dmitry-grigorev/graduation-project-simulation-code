library(poisson)
library(ggplot2)

intensity<-function(x){ pmin(x/2, 1) }
target<-50
rate<-1

K.epan<-Vectorize(function(x){ifelse(abs(x)<=1, 0.75*(1-x^2), 0)})

intensity.discontinuity.KE<-function(times, x, n, h)
{
  sum(K.epan( (times - x)/h )*sign(times-x))/n/h
}

source("C:\\Users\\Acer\\Documents\\courseworksimulation\\main_script.R")
library(permute)

set.seed(230022)

n<-1000

.inf<-999 # для случаев деления на ноль

permutation.test<-function(t, T.permuted, N, sign.level)
{
  k1<-floor(sign.level*N/2)
  k2<-N - k1
  
  M2<-sum(T.permuted[k2] < T.permuted)
  M1<-sum(T.permuted[k1] > T.permuted)
  M0<-{
    if(T.permuted[k1] != T.permuted[k2]){
      sum(T.permuted[k1] == T.permuted) + sum(T.permuted[k2] == T.permuted)}
    else {sum(T.permuted[k1] == T.permuted)}
  }
  a<-(sign.level*N - M2 - M1)/M0
  
  if(t > T.permuted[k2] | t < T.permuted[k1])
  { return(1) }
  else if(t == T.permuted[k2] | t == T.permuted[k1]){
    return(a)
  } 
  else
    return(0)
}

repl.prog<-0
median.permtest.simulation<-function(n, n.perm, q, d1, d2, x)
{
  T_values_perm<-matrix(c(0), nrow = 3, ncol = n.perm)
  S_values_perm<-matrix(c(0), nrow = 3, ncol = n.perm)
  t.test_values<-numeric(3)
  t_median<-numeric(3)
  st_median<-numeric(3)
  
  iters<- n.perm*3
  pb <- winProgressBar(title = "progress bar", min = 0,
  max = iters, width = 300)  
  
  indexes<-1:n
  quant<- -qnorm(0.025)
  
  x1labs<-sample(c(1,2), n, prob = c(q, 1-q), replace = TRUE)
  
  x1<-2*runif(n)-1
  x1[x1<0] <- -x1[x1<0]

  y1<- f.vec(x1) + rnorm(n)*ifelse(x1labs==1, sqrt(d1), sqrt(d2))
  y1[x1labs == 1]<- -y1[x1labs == 1]
  
  data<-data.frame(X = x1, Y = y1, W = x1labs)
  n1<-sum(data$W == 1)
  n2<-n-n1
  frac1<-n/n1
  frac2<-n/n2
  h<-c(0.1, 0.3, 0.5)
  
  for(ih in 1:3)
  {
    M1_median<-KE_median( data$X[data$W == 1], data$Y[data$W == 1], x ,h[ih], TRUE)
    M2_median<-KE_median( data$X[data$W == 2], data$Y[data$W == 2], x ,h[ih], TRUE) 
    VAR<-frac1*KE_var(data$X[data$W == 1], data$Y[data$W == 1], x, 0 , h[ih], 0.3) +
      frac2*KE_var(data$X[data$W == 2], data$Y[data$W == 2], x, 0 , h[ih], 0.3)
    t_median[ih]<-sqrt(n*h[ih])*(M2_median - 0 - h[ih]*(3/8/d2) - (M1_median + 0- h[ih]*(3/8/d1)))
    if(is.na(t_median[ih]))
      t_median[ih]<-.inf*sample(c(-1,1),1)
    st_median[ih]<-t_median[ih]/sqrt(VAR)
    if(is.nan(st_median[ih]))
      st_median[ih]<-.inf*sample(c(-1,1),1)
    t.test_values[ih]<-as.numeric( abs(st_median[ih])> quant)
  }
  
  
  progress.variable<-1
  for(m_i in 1:n.perm){
    permuted.indexes<-shuffle(indexes)
    
    W<-data$W[permuted.indexes]
    for(ih in 1:3)
    {
    M1_median<-KE_median( data$X[W == 1], data$Y[W == 1], x ,h[ih], TRUE)
    M2_median<-KE_median( data$X[W == 2], data$Y[W == 2], x ,h[ih], TRUE)
    VAR<-frac1*KE_var(data$X[W == 1], data$Y[W == 1], x, 0 ,h[ih], 0.3) +
      frac2*KE_var(data$X[W == 2], data$Y[W == 2], x, 0 ,h[ih], 0.3)
    t_median.perm<-sqrt(n*h[ih])*(M2_median - 0 - h[ih]*(3/8/d2) - (M1_median + 0- h[ih]*(3/8/d1)))
    if(is.na(t_median.perm))
      t_median.perm<-.inf*sample(c(-1,1),1)
    T_values_perm[ih, m_i]<-t_median.perm
    S_values_perm[ih, m_i]<-t_median.perm/sqrt(VAR)
    if(is.nan(S_values_perm[ih, m_i]))
      S_values_perm[ih, m_i]<-.inf*sample(c(-1,1),1)
    setWinProgressBar(pb, progress.variable, title=paste( 
     round(progress.variable/iters * 100, 0), "% of total done"))
    
    progress.variable<-progress.variable+1
    }
    
  }
  
  NS.result<-numeric(3)
  S.result<-numeric(3)
  for(ih in 1:3)
  {
    NS.result[ih]<-permutation.test(t_median[ih], sort(T_values_perm[ih,]), n.perm, 0.05)
    S.result[ih]<-permutation.test(st_median[ih], sort(S_values_perm[ih,]), n.perm, 0.05)
  }
  close(pb)
  repl.prog<-repl.prog+1
  print(repl.prog)
  list(NS.res.0.1 = NS.result[1], S.res.0.1 = S.result[1], t.test.res.0.1 = t.test_values[1],
       NS.res.0.3 = NS.result[2], S.res.0.3 = S.result[2], t.test.res.0.3 = t.test_values[2],
       NS.res.0.5 = NS.result[3], S.res.0.5 = S.result[3], t.test.res.0.5 = t.test_values[3])
}

res1<-replicate(500, {r<-median.permtest.simulation(200, 500, 0.5, 1, 1, 0)})
res2<-replicate(500, {r<-median.permtest.simulation(200, 500, 0.5, 5, 1, 0)})
res3<-replicate(500, {r<-median.permtest.simulation(200, 500, 0.3, 1, 1, 0)})
res4<-replicate(500, {r<-median.permtest.simulation(200, 500, 0.3, 5, 1, 0)})


results<-list(res1, res2, res3, res4)

for(k in 1:4)
{
  res<-results[[k]]
  for(i in 1:9)
  {
    print(mean(unlist(res[i,])))
  }
}

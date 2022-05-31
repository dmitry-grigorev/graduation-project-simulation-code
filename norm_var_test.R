source("C:\\Users\\Acer\\Documents\\courseworksimulation\\main_script.R")
library(ggplot2)

set.seed(230022)

start.n<-2500
stepsize<-2500
n<-10000
eps<-1e-6

ssize<-c(2500, 5000, 10000)

robustness.normnoise.simulation<-function(n, m, q, d, x){
  k<-4*length(ssize)
  T_values_mean<-matrix(c(0), ncol = m, nrow = k)
  T_values_median<-matrix(c(0), ncol = m, nrow = k)
  
  iters<- k * m
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = iters, width = 300)  # Создаём progress bar
  
  test.size.mean<-numeric(k)
  test.size.median<-numeric(k)  
  
  progress.variable<-1
  for(m_i in 1:m){
    x1labs<-sample(c(1,2), n, prob = c(0.5, 0.5), replace = TRUE)
    y1labs<-sample(c(0,1), n, prob = c(1-q, q), replace = TRUE)
    
    odds<-sum(y1labs == 1)
    
    x1<-rnorm(n, mean = 0)
    x1[x1<0] <- -x1[x1<0]
    
    error<-rnorm(n)
    
    y1<- f.vec(x1)
    y1[x1labs == 1]<- -y1[x1labs == 1]
    j<-1
    for(dev in c(1, 2, 5, 10))
    {
    .y1<- y1 + error*ifelse(y1labs == 0, 1, sqrt(dev))
    
    data<-data.frame(X = x1, Y = .y1, W = x1labs)
    for(n_i in ssize)
    {
      .data<-data[1:n_i,]

      h_opt<-1.06*sd(.data$X)*length(.data$X)^(-1/3)
      h_opt_neg<-h_opt
      h_opt_pos<-h_opt
      
      M1_mean<-KE_mean( .data$X[.data$W == 1], .data$Y[.data$W == 1], x ,h_opt_neg)
      M2_mean<-KE_mean( .data$X[.data$W == 2], .data$Y[.data$W == 2], x ,h_opt_pos)
      
      M1_median<-KE_median( .data$X[.data$W == 1], .data$Y[.data$W == 1], x ,h_opt_neg)
      M2_median<-KE_median( .data$X[.data$W == 2], .data$Y[.data$W == 2], x ,h_opt_pos)
      
      D1_mean.corr<-stat.correction.mean(.data$X[.data$W == 1], .data$Y[.data$W == 1], x, -0, h_opt_neg, eps = eps, M = K(0)/4)
      D2_mean.corr<-stat.correction.mean(.data$X[.data$W == 2], .data$Y[.data$W == 2], x, 0, h_opt_pos, eps = eps, M = K(0)/4)
      
      D1_median.corr<-stat.correction(.data$X[.data$W == 1], .data$Y[.data$W == 1], x, -0, h_opt_neg, eps = eps, M = K(0)/4)
      D2_median.corr<-stat.correction(.data$X[.data$W == 2], .data$Y[.data$W == 2], x, 0, h_opt_pos, eps = eps, M = K(0)/4)
      
      t_mean<-sqrt(n_i*h_opt)*(M2_mean - 0 - h_opt_pos*D2_mean.corr - (M1_mean + 0 - h_opt_neg*D1_mean.corr))
      t_median<-sqrt(n_i*h_opt)*(M2_median - 0 - h_opt_pos*D2_median.corr - (M1_median + 0 - h_opt_neg*D1_median.corr))

      setWinProgressBar(pb, progress.variable, title=paste( 
        round(progress.variable/iters * 100, 0), "% of total done"))
      
      progress.variable<-progress.variable+1
      
      T_values_mean[j, m_i]<-t_mean
      T_values_median[j, m_i]<-ifelse(is.infinite(t_median), NA, t_median)
      j<-j+1
    }
    }
  }
  close(pb)
  

  list(var_mean = apply(T_values_mean, 1, var), var_median = apply(T_values_median, 1, var, na.rm = TRUE))
}

results.nonoise<-robustness.normnoise.simulation(10000, 1000, 0, 1, 0)


noise<-c(0.01, 0.05, 0.1, 0.2, 0.5)
vars<-c(1, 2, 5, 10)

for(a in noise){
    res<-robustness.normnoise.simulation(10000, 1000, a, 1, 0)
    print(c("Level: ", a))
    print(res$var_mean)
    print(res$var_median)
}


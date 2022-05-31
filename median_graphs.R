source("C:\\Users\\Acer\\Documents\\courseworksimulation\\main_script.R")
library(ggplot2)

set.seed(23002)

start.n<-2500
stepsize<-2500
n<-10000
eps<-1e-6

variance.test.size.simulation<-function(n, m, q, d1, d2, x){
  T_values_median.nest<-matrix(c(0), ncol = m, nrow = n)
  T_values_median.west<-matrix(c(0), ncol = m, nrow = n)
  
  iters<-(n-start.n+stepsize)/stepsize * m
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = iters, width = 300)  # Создаём progress bar
  
  test.median.nest2500<-numeric(m)
  test.median.west2500<-numeric(m)
  test.median.nest5000<-numeric(m)
  test.median.west5000<-numeric(m)
  test.median.nest10000<-numeric(m)
  test.median.west10000<-numeric(m)

  
  progress.variable<-1
  for(m_i in 1:m){
    x1labs<-sample(c(1,2), n, prob = c(q, 1-q), replace = TRUE)
    
    x1<-rnorm(n, mean = 0)
    x1[x1<0] <- -x1[x1<0]
    #x1[x1labs == 1]<- -x1[x1labs == 1]
    y1<- f.vec(x1) + rnorm(n)*ifelse(x1labs==1, sqrt(d1), sqrt(d2))
    y1[x1labs == 1]<- -y1[x1labs == 1]
    
    data<-data.frame(X = x1, Y = y1, W = x1labs)
    for(n_i in seq.int(start.n, n, stepsize))
    {
      .data<-data[1:n_i,]

      h_opt<-1.06*sd(.data$X)*length(.data$X)^(-1/3)
      h_opt_neg<-h_opt
      h_opt_pos<-h_opt
      
      
      M1_median<-KE_median( .data$X[.data$W == 1], .data$Y[.data$W == 1], x ,h_opt_neg, TRUE)
      M2_median<-KE_median( .data$X[.data$W == 2], .data$Y[.data$W == 2], x ,h_opt_pos, TRUE)
      
      D1_median.corr<-stat.correction(.data$X[.data$W == 1], .data$Y[.data$W == 1], x, -0, h_opt_neg, eps = eps, M = K(0)/4)
      D2_median.corr<-stat.correction(.data$X[.data$W == 2], .data$Y[.data$W == 2], x, 0, h_opt_pos, eps = eps, M = K(0)/4)

      t_mean <- 0
      t_median.nest<-sqrt(n_i*h_opt)*(M2_median - 0 - h_opt_pos*(3/8/d2) - (M1_median + 0- h_opt_neg*(3/8/d1)))
      t_median.west<-sqrt(n_i*h_opt)*(M2_median - 0 - h_opt_pos*D2_median.corr - (M1_median + 0 - h_opt_neg*D1_median.corr))
      T_values_median.nest[n_i,m_i]<-t_median.nest
      T_values_median.west[n_i, m_i]<-t_median.west
      setWinProgressBar(pb, progress.variable, title=paste( 
        round(progress.variable/iters * 100, 0), "% of total done"))
      
      progress.variable<-progress.variable+1
      
      if(n_i == 2500)
      { test.median.nest2500[m_i]<-t_median.nest; test.median.west2500[m_i]<-ifelse(is.infinite(t_median.west), NA, t_median.west) }
      
      if(n_i == 5000)
      { test.median.nest5000[m_i]<-t_median.nest; test.median.west5000[m_i]<-ifelse(is.infinite(t_median.west), NA, t_median.west) }
      
      if(n_i == 10000)
      { test.median.nest10000[m_i]<-t_median.nest; test.median.west10000[m_i]<-ifelse(is.infinite(t_median.west), NA, t_median.west) }

    }
  }
  close(pb)
  
  T_values_median.west[is.infinite(T_values_median.west)]<-NaN
  
  list(var_median.nest = apply(T_values_median.nest, 1, var), 
       var_median.west = apply(T_values_median.west, 1, var, na.rm = TRUE),
       test.median.nest2500 = test.median.nest2500, test.median.west2500 = test.median.west2500, 
       test.median.nest5000 = test.median.nest5000, test.median.west5000 = test.median.west5000,
       test.median.nest10000 = test.median.nest10000, test.median.west10000 = test.median.west10000)
}

num.seq<-seq.int(start.n, n, stepsize)

var1<-c(1, 5)
q1<-c(0.5, 0.25)

results11<-variance.test.size.simulation(n, 1000, q1[1], var1[1], 1, 0)
results12<-variance.test.size.simulation(n, 1000, q1[2], var1[1], 1, 0)
results21<-variance.test.size.simulation(n, 1000, q1[1], var1[2], 1, 0)
results22<-variance.test.size.simulation(n, 1000, q1[2], var1[2], 1, 0)

results<-list(results11, results12, results21, results22)

median.limitvar<-function(q, v){(v/q+1/(1-q))*0.3*pi*sqrt(2*pi)}

test.size.2500<-numeric(4)
test.size.5000<-numeric(4)
test.size.10000<-numeric(4)
for(i in 1:2)
{
  for(j in 1:2)
  {
    test.size.2500[2*(i-1)+j]<-mean(abs(results[[2*(i-1)+j]]$test.median.nest2500) > -qnorm(0.025, sd = sqrt(median.limitvar(q1[j], var1[i]))))
    test.size.5000[2*(i-1)+j]<-mean(abs(results[[2*(i-1)+j]]$test.median.nest5000) > -qnorm(0.025, sd = sqrt(median.limitvar(q1[j], var1[i]))))
    test.size.10000[2*(i-1)+j]<-mean(abs(results[[2*(i-1)+j]]$test.median.nest10000) > -qnorm(0.025, sd = sqrt(median.limitvar(q1[j], var1[i]))))

      }
}

plot.var<-function(i, j)
{
  ggplot()+geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.nest[num.seq]))+
    geom_abline(slope = 0, intercept = median.limitvar(q1[j], var1[i]), colour = "blue")+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.nest[num.seq]*999/qchisq(0.975, df = 999), colour = "green"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.nest[num.seq]*999/qchisq(0.025, df = 999), colour = "green"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.nest[num.seq]*999/qchisq(0.9, df = 999), colour = "red"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.nest[num.seq]*999/qchisq(0.1, df = 999), colour = "red"))+
    xlab("Объём выборки") + ylab("Величина дисперсии")+theme(legend.position="none")
  
}

plot.var.west<-function(i, j)
{
  ggplot()+geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.west[num.seq]))+
    geom_abline(slope = 0, intercept = median.limitvar(q1[j], var1[i]), colour = "blue")+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.west[num.seq]*999/qchisq(0.975, df = 999), colour = "green"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.west[num.seq]*999/qchisq(0.025, df = 999), colour = "green"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.west[num.seq]*999/qchisq(0.9, df = 999), colour = "red"))+
    geom_line(aes(x = num.seq, y = results[[2*(i-1)+j]]$var_median.west[num.seq]*999/qchisq(0.1, df = 999), colour = "red"))+
    xlab("Объём выборки") + ylab("Величина дисперсии")+theme(legend.position="none")
  
}

plot.var(1,1)
plot.var(1,2)
plot.var(2,1)
plot.var(2,2)

plot.var.west(1,1)
plot.var.west(1,2)
plot.var.west(2,1)
plot.var.west(2,2)

plot.hist.west<-function(i, j)
{
  xval<-results[[2*(i-1)+j]]$test.median.west10000[is.finite(results[[2*(i-1)+j]]$test.median.west10000)]
  ggplot()+geom_histogram(aes(x = xval, y=..density..), bins = 25, fill="#69b3a2", color="#e9ecef", alpha=0.9, closed = "left")+ 
  geom_function(fun = dnorm, args = list(mean = 0, sd = sqrt(median.limitvar(q1[j], var1[i]))))+
  xlab("Значения статистики")+ylab("Плотность")+geom_vline(xintercept = 0, col = "red")+
  geom_vline(xintercept = 0, col = "blue")#+ scale_x_continuous(breaks=seq(-3, 4, 1))
}

plot.hist.west(1,1)
plot.hist.west(1,2)
plot.hist.west(2,1)
plot.hist.west(2,2)



source("C:\\Users\\Acer\\Documents\\courseworksimulation\\pois_script.R")


set.seed(230025)

start.n<-2500
stepsize<-2500
n<-20000
eps<-1e-6

variance.test.size.simulation<-function(n, m, true.disc, x)
{
  T_values_pois.nest<-matrix(c(0), ncol = m, nrow = n/stepsize)
  
  iters<- n/stepsize * m
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = iters, width = 300)  
  
  test.pois.nest2500<-numeric(m)
  test.pois.nest5000<-numeric(m)
  test.pois.nest10000<-numeric(m)
  test.pois.nest20000<-numeric(m)
  
  progress.variable<-1
  for(m_i in 1:m){
    n.index1<-numeric(n/stepsize)
    n.index2<-numeric(n/stepsize)
    obs1<-c(0)
    obs2<-c(0)
    
    for(i in 1:8)
    {
      .obs1<-as.vector(replicate(stepsize/2, nhpp.event.times(rate, target, intensity)))
      .obs2<-as.vector(replicate(stepsize - stepsize/2, nhpp.event.times(rate, target, intensity)))
      .obs1<-.obs1[.obs1>0 & .obs1<2]
      .obs2<-.obs2[.obs2>0 & .obs2<2]
      n.index1[i]<-length(.obs1)
      n.index2[i]<-length(.obs2)
      obs1<-c(obs1, .obs1)
      obs2<-c(obs2, .obs2)
    }
    
    n.index1<-cumsum(n.index1)
    n.index2<-cumsum(n.index2)
    
    obs1<-obs1[-1]
    obs2<-obs2[-1]
    
    for(i in 1:(n/stepsize))
    {
      h<-(i*stepsize)^(-1/3)
      
      times1<-obs1[1:(n.index1[i])]
      times2<-obs2[1:(n.index2[i])]
      
      disc1<-intensity.discontinuity.KE(times1, x, (i*stepsize)/2, h)
      disc2<-intensity.discontinuity.KE(times2, x, i*stepsize - (i*stepsize)/2, h)
      
      test.value<-sqrt(i*stepsize*h)*(disc1 + disc2 + true.disc - h*3/8)
      
      T_values_pois.nest[i, m_i]<-test.value
      
      setWinProgressBar(pb, progress.variable, title=paste( 
        round(progress.variable/iters * 100, 0), "% of total done"))
      
      progress.variable<-progress.variable+1
      
      if(i == 1)
      {
        test.pois.nest2500[m_i]<-test.value
      }
      if(i == 2)
      {
        test.pois.nest5000[m_i]<-test.value
      }
      if(i == 4)
      {
        test.pois.nest10000[m_i]<-test.value
      }
      if(i == 8)
      {
        test.pois.nest20000[m_i]<-test.value
      }
      
    }
  }
  close(pb)
  
  print(T_values_pois.nest)
  
  list(var_pois = apply(T_values_pois.nest, 1, var), 
       test.pois.2500 = test.pois.nest2500,
       test.pois.5000 = test.pois.nest5000,
       test.pois.10000 = test.pois.nest10000,
       test.pois.20000 = test.pois.nest20000
       )
}
num.seq<-seq.int(2500, 20000, 2500)

pois.test<-variance.test.size.simulation(n, 2000, 0, 1)

ggplot()+geom_line(aes(x = num.seq, y = pois.test$var_pois))+
  geom_abline(slope = 0, intercept = 1.2, colour = "blue")+
  geom_line(aes(x = num.seq, y = pois.test$var_pois*999/qchisq(0.975, df = 999), colour = "green"))+
  geom_line(aes(x = num.seq, y = pois.test$var_pois*999/qchisq(0.025, df = 999), colour = "green"))+
  geom_line(aes(x = num.seq, y = pois.test$var_pois*999/qchisq(0.9, df = 999), colour = "red"))+
  geom_line(aes(x = num.seq, y = pois.test$var_pois*999/qchisq(0.1, df = 999), colour = "red"))+
  xlab("Объём выборки") + ylab("Величина дисперсии")+theme(legend.position="none")

ggplot()+geom_histogram(aes(x = pois.test$test.pois.20000, y=..density..), 
                        binwidth = 0.3, fill="#69b3a2", color="#e9ecef", alpha=0.9, closed = "left", boundary = -3.5)+ 
  geom_function(fun = dnorm, args = list(mean = 0, sd = sqrt(1.2)))+xlab("Значения статистики")+
  ylab("Плотность")+geom_vline(xintercept = 0, col = "red")+
  geom_vline(xintercept = 0, col = "blue")+ scale_x_continuous(breaks=seq(-3, 4, 1))

ggplot()+geom_histogram(aes(x = pois.test$test.pois.20000+3/8, y=..density..), 
                        binwidth = 0.3, fill="#69b3a2", color="#e9ecef", alpha=0.9, closed = "left", boundary = -3.5+3/8)+ 
  geom_function(fun = dnorm, args = list(mean = 3/8, sd = sqrt(1.2)))+xlab("Значения статистики")+
  ylab("Плотность")+geom_vline(xintercept = 0, col = "red")+
  geom_vline(xintercept = 3/8, col = "blue")+ scale_x_continuous(breaks=seq(-3, 4, 1))

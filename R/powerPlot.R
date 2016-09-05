powerPlot<-function(power.object){
  
  require(ggplot2)
  
  error.rates<-seq(0.001, 1, by=0.005)
  model.power<-vector(length=length(error.rates))
  nsim.m<-length(power.object$imppvals)
  for(i in 1:length(error.rates)){
    model.power[i]<-(length(which(power.object$imppvals<=error.rates[i]))/nsim.m)*100
  }
  
  dat<-data.frame(err=error.rates, pow=model.power)
  p <-ggplot(dat) + geom_line(aes(err, pow), color='blue', size=1) + xlab('Error Rate') + ylab('Power') + theme_bw() + scale_y_continuous(limits=c(0,100)) + geom_vline(xintercept=0.05, linetype=2)+ geom_vline(xintercept=0.01, linetype=2)
  return(p)
}


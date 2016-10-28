
plotVariance<-function(model, simdat){
  ow<-options("warn")
  options(warn=1)
  
  dat<-model$data
  fits<-fitted(model)
  cut.fit<-cut(fits, breaks=quantile(fits, probs=seq(0,1,length=20)))
  means.mod<-tapply(fits, cut.fit, mean)
  vars.mod<-tapply(residuals(model, type='pearson'), cut.fit, var)
  
  disp<-vector(length=100)
  means.sim=vars.sim=matrix(NA, nrow=length(means.mod), ncol=100)
  for(i in 1:100){
    dat$response<-simdat[,i]
    sim<-update(model,  .~ ., data=dat)
    fitsim<-fitted(sim)
    means.sim[,i]<-tapply(fitsim, cut.fit, mean)
    vars.sim[,i]<-tapply(residuals(sim, type='pearson'), cut.fit, var)
    disp[i]<-summary(sim)$dispersion
  }
  
  test<-data.frame(FittedValues= as.vector(na.omit(means.sim)), 
                   Variance=as.vector(na.omit(vars.sim)), 
                   gp=as.factor(rep(1:100, each=19)))
  test_real<-data.frame(FittedValues= means.mod, 
                        Variance=vars.mod, 
                        gp=rep(1, nrow(means.mod)))
  library(ggplot2)
  library(gridExtra)
  
  p <- ggplot(data=test, aes(x=FittedValues, y=Variance,group=gp)) + 
    geom_path(alpha=0.1)
  p <- p + geom_path(data = test_real, aes(x=FittedValues, y=Variance,
                                           colour='red'), size=1) + theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=
            element_blank()) + xlab('Mean Fitted Values') + ylab("Residual Variance") + guides(colour=FALSE) + geom_rug(data=test, sides='b', alpha=1/5)
  
  if(summary(model)$dispersion!=1){
    p2<-ggplot(data.frame(disp), aes(disp)) + geom_histogram(colour='black', fill='lightgrey', bins=15) + theme_bw() + xlab("Estimated Dispersion") + ylab("Frequency") + ggtitle("") + geom_vline(xintercept = summary(model)$dispersion, colour='blue', size=2 , linetype=2)
    
    grid.arrange(p, p2)  
  }else{
    print(p)
  }
  
  options(ow)
}

plotMean<-function(simdat, model){
  ow<-options("warn")
  options(warn=1)
  
  dat<-model$data
  # hist(apply(simdat, 2, mean), main='', xlab='mean(response)', xlim=range(c(apply(simdat, 2, mean), mean(dat$response), mean(fitted(model)))))
  # abline(v=mean(dat$response), col='red', lwd=2, lty=4)
  # abline(v=mean(fitted(model)), col='blue', lwd=2, lty=3)  
  # 
  dat2<-apply(simdat, 2, mean)
  xlims=range(c(dat2, mean(dat$response), mean(fitted(model))))

  p3<-ggplot(data.frame(dat2), aes(dat2)) + geom_histogram(colour='black', fill='lightgrey', bins=15) + theme_bw() + xlab("Mean") + ylab("Frequency") + ggtitle("")  + geom_vline(xintercept = mean(fitted(model)), colour='blue', size=2, linetype=4 ) + geom_vline(xintercept=mean(dat$response), colour='red', size=2, linetype=3) + xlim(xlims) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  print(p3)
  options(ow)
}

#' Plotting function to assess the mean-variance relationship of generated data
#'
#' @param model model object of class glm, gam or gamMRSea
#' @param simdat matrix of simulated response data.  Each column is a new simulated data set
#' @param cuts The number of cut points for calculation of residual variance
#' @param quants vector of length 2 stating the quantiles for the confidence interval bands for the simulated mean-variance relationship
#'
#'
#' @author Lindesay Scott-Hayward (University of St Andrews)
#'
#' @export
#'


plotVariance<-function(model, simdat, cuts=20, quants=c(0.025, 0.975)){
  ow<-options("warn")
  options(warn=1)

  dat<-model$data
  fits<-fitted(model)
  cut.fit<-cut(fits, breaks=quantile(fits, probs=seq(0,1,length=cuts)))
  means.mod<-tapply(fits, cut.fit, mean)
  vars.mod<-tapply(residuals(model, type='response'), cut.fit, var)

  if(ncol(simdat)>100){
    nsim=100
  }else{
    nsim=ncol(simdat)
  }

  disp<-vector(length=nsim)
  means.sim=vars.sim=matrix(NA, nrow=length(means.mod), ncol=nsim)
  for(i in 1:nsim){
    dat$response<-simdat[,i]
    sim<-update(model,  .~ ., data=dat)
    fitsim<-fitted(sim)
    means.sim[,i]<-tapply(fitsim, cut.fit, mean)
    vars.sim[,i]<-tapply(residuals(sim, type='response'), cut.fit, var)
    disp[i]<-summary(sim)$dispersion
  }

  meansimfits<-apply(means.sim, 1, mean)
  varsimfits.m<-apply(vars.sim, 1, mean)
  varsimfits.quant<-data.frame(t(apply(vars.sim, 1, quantile, probs=quants)))

  newdat<-data.frame(meansimfits, varsimfits.m, varsimfits.quant)

  test<-data.frame(meansimfits= as.vector(na.omit(means.sim)),
                   varsimfits.m=as.vector(na.omit(vars.sim)))
  test_real<-data.frame(meansimfits= means.mod,
                        varsimfits.m=vars.mod,
                        gp=rep(1, nrow(means.mod)))
  library(ggplot2)
  library(gridExtra)

  p <- ggplot(data=newdat) + geom_line(aes(x=log(meansimfits), y=varsimfits.m)) + geom_line(data = test_real, aes(x=log(meansimfits), y=varsimfits.m), colour='red', size=1) + theme_bw() + xlab('Log Mean Fitted Values') + ylab("Residual Variance") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=
            element_blank()) + geom_ribbon(aes(x=log(meansimfits), y=varsimfits.m, ymin=X2.5., ymax=X97.5.), colour='grey', fill='grey', alpha=0.3)+ guides(colour=FALSE) + geom_line(data=test_real, aes(x=log(meansimfits),y=summary(model)$dispersion*meansimfits), colour='blue') +  geom_rug(data=test_real, aes(x=log(meansimfits)), sides='b', alpha=1/5)

  # p<-ggplot(data=newdat) + geom_line(aes(x=meansimfits, y=varsimfits.m)) + geom_ribbon(aes(x=meansimfits, y=varsimfits.m, ymin=X2.5., ymax=X97.5.), colour='grey', fill='grey', alpha=0.3)
  # p <- p + geom_line(data = test_real, aes(x=meansimfits, y=varsimfits.m), colour='red', size=1) + theme_bw() +
  #   theme(panel.grid.major=element_blank(), panel.grid.minor=
  #           element_blank()) + xlab('Mean Fitted Values') + ylab("Residual Variance") + guides(colour=FALSE) + geom_abline(slope = summary(model)$dispersion, intercept = 0, colour='blue')
  # p<-p + geom_rug(data=test, aes(x=meansimfits), sides='b', alpha=1/5)

  if(summary(model)$dispersion!=1){
    p2<-ggplot(data.frame(disp), aes(disp)) + geom_histogram(colour='black', fill='lightgrey', bins=15) + theme_bw() + xlab("Estimated Dispersion") + ylab("Frequency") + ggtitle("") + geom_vline(xintercept = summary(model)$dispersion, colour='blue', size=2 , linetype=2)

    grid.arrange(p, p2)
  }else{
    print(p)
  }

  options(ow)
}


#' Plotting function to assess the mean of generated data
#'
#' @param model model object of class glm, gam or gamMRSea
#' @param simdat matrix of simulated response data.  Each column is a new simulated data set
#' @param bins number of bins to be plotted on the histogram
#' @param quants vector of length 2 stating the quantiles for the confidence interval bands for the simulated data
#'
#'
#' @author Lindesay Scott-Hayward (University of St Andrews)
#'
#' @export
#'
plotMean<-function(model, simdat, bins=15, quants=c(0.025, 0.975)){
  ow<-options("warn")
  options(warn=1)

  dat<-model$data

  dat2<-apply(simdat, 2, mean)
  cis<-quantile(dat2, probs=quants)
  xlims=range(c(dat2, mean(dat$response), mean(fitted(model))))

  p3<-ggplot(data.frame(dat2), aes(dat2)) + geom_histogram(colour='black', fill='lightgrey', bins=bins) + theme_bw() + xlab("Mean") + ylab("Frequency") + ggtitle("")  + geom_vline(xintercept = mean(fitted(model)), colour='blue', size=2, linetype=4 ) + geom_vline(xintercept=mean(dat$response), colour='red', size=2, linetype=3) + xlim(xlims) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())  + geom_vline(xintercept=cis[1], colour='grey', size=2, linetype=2) + geom_vline(xintercept=cis[2], colour='grey', size=2, linetype=2)
  print(p3)
  options(ow)
}


#' Function to plot the confidence intervals from the change parameter
#'
#' @param power.object A power analysis object of class gamMRSea.power
#' @param truebeta.response The true coeficient of the change covariate on the scale of the response. E.g. 0.8 is a 20\% decline in overall numbers.
#' @param alternativebeta An alternative coefficient which you wish to test if it is a plausible value for the truth (i.e. sits within the confidence interval for the coefficient). Defaults to no change.
#' @param plot (default = TRUE).  Logical stating whether the results should be plotted.
#'
#' @examples
#'
#' # Load a power object that was created with a 20\% sitewide decline
#' # 50 simulations were run
#' data(nysted.power.oc)
#'
#' plot.coverage(nysted.power.oc, truebeta.response=0.8)
#'
#' @author LAS Scott-Hayward, University of St Andrews
#'
#' @export
#'
plot.coverage<-function(power.object, truebeta.response, alternativebeta=NULL, plot=TRUE){
  invlink<-power.object$link$linkinv
  link<-power.object$link$link
  responseCIs<-invlink(power.object$betacis)
  if(is.null(alternativebeta)){
    alternativebeta<-invlink(0)
  }
  if(plot){
    plot(1, responseCIs[1,1], ylim=range(responseCIs), xlim=c(0,dim(responseCIs)[1]), pch='.', xlab='Simulation Number', ylab='Response Scale Beta Coefficient')
    segments(x0 = 1:dim(responseCIs)[1], y0 = responseCIs[,1], x1 = 1:dim(responseCIs)[1], y1 = responseCIs[,2])
    abline(h=truebeta.response, col='grey', lty=2)
    abline(h=alternativebeta, col='darkgrey', lty=3)
  }
  # make a table
  model.coverage<-impact.coverage(link(truebeta.response), power.object$betacis)*100
  alt.model.coverage<-impact.coverage(link(alternativebeta), power.object$betacis)*100

  print(data.frame(Beta=c('True', 'Alternative'), 'Pct CIs inc Beta'=c(model.coverage, alt.model.coverage)))

}



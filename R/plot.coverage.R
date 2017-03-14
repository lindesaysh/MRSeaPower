#' Function to plot the confidence intervals from the change parameter
#'
#' @param power.object A power analysis object of class gamMRSea.power
#'
#' @examples
#'
#' # Load a power object that was created with a 20\% sitewide decline
#' # 50 simulations were run
#' data(nysted.power.oc)
#'
#' plot.coverage(nysted.power.oc)
#'
#' @author LAS Scott-Hayward, University of St Andrews
#'
#' @export
#'
plot.coverage<-function(power.object){
  plot(1, nysted.power.oc$betacis[1,1], ylim=range(nysted.power.oc$betacis), xlim=c(0,50), pch='.', xlab='Simulation Number', ylab='Beta Coefficient Value')
  segments(x0 = 1:dim(nysted.power.oc$betacis)[1], y0 = nysted.power.oc$betacis[,1], x1 = 1:dim(nysted.power.oc$betacis)[1], y1 = nysted.power.oc$betacis[,2])
  abline(h=log(0.8), col='grey', lty=2)
}



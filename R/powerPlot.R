#' Function to show how the power to detect change varies with the error rate.
#'
#' @param power.object An object of class gamMRSea.power
#'
#' @return
#' Plot showing how the power to detect change varies with increasing error rate. The first grey dashed line is at 1% and the second at 5%, traditionally values used as p-value cutoffs. The blue dashed lines indicate the error rate required to get a power of 80%. This value is given in the title.
#'
#' @examples
#'
#' # Load two power objects that were created with a 20% sitewide decline
#' # 50 simulations were run
#' data(nysted.power.oc)
#' data(nysted.power.oc.null)
#'
#' powerPlot(nysted.power.oc)
#'
#' @author LAS Scott-Hayward
#'
#' @export


powerPlot<-function(power.object){

  require(ggplot2)

  error.rates<-seq(0.001, 1, by=0.005)
  model.power<-vector(length=length(error.rates))
  nsim.m<-length(power.object$imppvals)
  for(i in 1:length(error.rates)){
    model.power[i]<-(length(which(power.object$imppvals<=error.rates[i]))/nsim.m)*100
  }

  dat<-data.frame(err=error.rates, pow=model.power)

  id80<-which(dat$pow>79.5)[1]
  if(length(id80)>=1){
    dat80<-as.vector(dat[id80[1],])
  }else{
    dat80<-as.vector(dat[1,])
  }

  maintitle = paste('For power = 80%, error rate = ', dat80[1,1], sep='')

  p <-ggplot(dat) + geom_line(aes(err, pow), color='blue', size=1) + xlab('Error Rate') + ylab('Power') + theme_bw() + scale_y_continuous(limits=c(0,100)) + geom_vline(xintercept=0.05, linetype=2, colour='grey')+ geom_vline(xintercept=0.01, linetype=2, colour='grey') + geom_hline(yintercept=dat80[1,2], linetype=4, colour='blue') + geom_vline(xintercept=dat80[1,1], linetype=4, colour='blue') + ggtitle(maintitle)
  return(p)
}


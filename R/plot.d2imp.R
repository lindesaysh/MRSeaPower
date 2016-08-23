#' Function to plot the impact effect with distance to windfarm
#'
#' @param
#' @param
#'
#'
#'
#' @author Lindesay Scott-Hayward
#'

plot.d2imp<-function(power.object){

  data<-power.object$d2imp.plotdata

  p<-ggplot(data.frame(data), aes(x=MeanDist, y=bootMeanDiff, ymin=LowerCI, ymax=UpperCI, colour='red', fill='red')) + geom_line() + geom_ribbon(alpha=0.3) + theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15)) + labs(x='Distance from Origin', y='Difference (post-pre)') + theme(legend.position='none') + geom_hline(yintercept = 0)
  p

  return(p)
}

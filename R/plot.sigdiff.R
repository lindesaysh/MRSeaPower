#' Function to plot the output from powerSim
#'
#' @param
#'
#' @example
#'

plot.sigdiff<-function(powerout, coordinates, tailed='two', error.rate=0.05, family=FALSE){
require(ggplot2)
  nsim=length(powerout$rawrob)
  if(family==FALSE){
    pmatrixfull=matrix(unlist(powerout$significant.differences$inividual), ncol=2*nsim)
    t<-'Individual'
  }else{
    pmatrixfull=matrix(unlist(powerout$significant.differences$family), ncol=2*nsim)
    t<-'Family'
  }

  pmatrix<-switch(tailed,
                  two=pmatrixfull[,-seq(1, 2*nsim, by=2)],
                  one=pmatrixfull[,seq(1, 2*nsim, by=2)]
  )


  sigmat<-ifelse(pmatrix<=error.rate, 1, 0)
  propsig<-apply(sigmat, 1, sum)/nsim

  if(is.null(nrow(coordinates))) stop("Error: coordinates does not contain two columns")
  if(dim(coordinates)[2]>2) warning('Warning: more than two coordinate columns, first two used')

  plotdata<-data.frame(coordinates, proportion=propsig)
  p<-ggplot(plotdata)
  p<-p + geom_tile(aes(x=coordinates[,1], y=coordinates[,2], fill=proportion), height=500, width=500) + theme_bw() + coord_equal() + xlab(names(coordinates)[1]) + ylab(names(coordinates)[2]) + ggtitle(paste(t, ' ', error.rate*100,'% Error Rate', sep='')) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ scale_fill_continuous(name='Proportion')
  p

  return(p)
}

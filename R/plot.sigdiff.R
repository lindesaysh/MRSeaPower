#' Function to plot the output from powerSim
#'
#' @param powerout power analysis object of class 'gamMRSea.power'
#' @param coordinates nx2 dataframe or matrix of data coordinates
#' @param tailed takes on the value 'one' or 'two' depending on the type of test required
#' @param error.rate number stating the error rate required, default is 0.05 (5\% error rate)
#' @param adjustment character stating the adjustment, if any, for multiple testing.  The default value is 'none'.  Alternatives are 'sidak' and 'bonferroni'.
#'
#' @author Lindesay Scott-Hayward (University of St Andrews)
#'
#' @export
#'
plot.sigdiff<-function(powerout, coordinates, tailed='two', error.rate=0.05, adjustment='none'){

  require(ggplot2)

  nsim=length(powerout$imppvals)
  # individual or family
  if(adjustment=='none'){
    pmatrixfull=matrix(unlist(powerout$significant.differences$individual), ncol=2*nsim, nrow=dim(coordinates)[1])
    t<-'Individual'
    adjustment='none'
    maintitle<-paste(t, ' ', error.rate*100,'% Error Rate, Adj = "',adjustment, '"', sep='')
  }else{

    t<-'Family'
    er.orig<-error.rate
    if(adjustment=='sidak' | adjustment=='bonferroni'){
      pmatrixfull=matrix(unlist(powerout$significant.differences$individual), ncol=2*nsim, nrow=dim(coordinates)[1])
      error.rate<-switch(adjustment,
                         sidak=round(1-(1-error.rate)^(1/nrow(coordinates)),4),
                         none=error.rate,
                         bonferroni=error.rate/nrow(coordinates)
      )
      maintitle<-paste(t, ' ER = ', er.orig*100,'%; Invididual ER = ', error.rate*100,'% , Adj = "',adjustment, '"', sep='')
    }

  }

  # one or two tailed test
  pmatrix<-switch(tailed,
                  two=pmatrixfull[,-seq(1, 2*nsim, by=2)],
                  one=pmatrixfull[,seq(1, 2*nsim, by=2)]
  )

  # # calculate significance
  # if(adjustment=='bonferroni'){
  #   pmatrix<-pmatrix*(nrow(coordinates))
  # }

  sigmat<-ifelse(pmatrix<=error.rate, 1, 0)
  propsig<-(apply(sigmat, 1, sum)/nsim)*100

  if(is.null(nrow(coordinates))) stop("Error: coordinates does not contain two columns")
  if(dim(coordinates)[2]>2) warning('Warning: more than two coordinate columns, first two used')

  height<-unique(coordinates[,2])[2]-unique(coordinates[,2])[1]
  width<-unique(coordinates[,1])[2]-unique(coordinates[,1])[1]

  plotdata<-data.frame(coordinates, percentage=propsig)
  p<-ggplot(plotdata)
  p<-p + geom_tile(aes(x=coordinates[,1], y=coordinates[,2], fill=percentage), height=height, width=width) + theme_bw() + coord_equal() + xlab(names(coordinates)[1]) + ylab(names(coordinates)[2]) + ggtitle(maintitle) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ scale_fill_continuous(name='%')
  p

  return(p)
}

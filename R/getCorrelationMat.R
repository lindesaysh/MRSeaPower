#' Function to get the correlation for each panel defined in the data
#'
#' @param panel vector identifying the panel to which each row of data belongs.  Data are considered correlated within a panel and independent between panels.
#' @param data vector of response data.
#' @param dots (default = FALSE). Logical flag stating whether to show progress dots printed in workspace.
#'
#' data(nystedA_slim)
#' nysted$panelid<-as.numeric(nysted$unique.transect.label)
#'
#' corrs<-getCorrelationMat(panel = nysted$panelid, data=nysted$response, dots = FALSE)
#'
#' @author LAS Scott-Hayward, University of St Andrews
#'
#' @export

getCorrelationMat<-function (panel, data, dots=FALSE)
{
  blocktab <- table(panel)
  overallacf<-acf(data, lag.max = max(blocktab), plot=F)$acf
  acfmat <- matrix(NA, length(unique(panel)), max(blocktab))

  for (i in 1:length(unique(panel))) {
    if(dots){if((iter/100)%%1 == 0){cat(iter, '\n')}else{cat('.')}}
    corr <- as.vector(acf(data[which(panel == unique(panel)[i])], plot = F, lag.max = max(blocktab))$acf)
    if(length(which(is.na(corr)))>0){
      #corr<-c(1, rep(0.999999, (length(corr)-1)))
      #corr<-seq(1, 0.001, length.out=(length(corr)))
      corr<-overallacf[1:length(corr)]
    }
    acfmat[i, 1:length(corr)] <- corr
  }
  #return(list(acfmat = acfmat, blocktab = blocktab))
  return(acfmat)
}

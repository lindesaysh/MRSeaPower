#'
#' Function to impose a re-distribution impact effect
#'
#' @param model A glm or gamMRSea model object.
#' @param data Data frame used to fit \code{model}
#' @param changecoef.link Coefficient(s) on the scale of the link function to determine the change in the population.  If of length 1, then the change occurs in the impact zone, whilst maintaining overall population size.  If of length 2, then there is a redistribution of the population with an an additional overall decrease.
#' @param panels Character vector denoting the column of \code{data} containing the panel structure.
#' @param imppoly A data frame containing the coordinates of a polygon defining the region of impact. The variable names must match the coordinate system in the data.
#'
#'@return Returns a data frame twice the size of the original input with additional columns for panel id, eventphase and truth.  The truth column represents the input data for the first half (pre-event; eventphase==0) and the second half has the change imposed (post-event; eventphase==1).
#'
#' @author Lindesay Scott-Hayward, University of St Andrews
#'
#' @export
#'


# imppoly
# $x.pos
# [1] 510210.0 510944.6 511353.8 511166.5 510244.5
#
# $y.pos
# [1] 6556624 6556997 6555979 6554563 6556177

genRedistData<-function(model, data, changecoef.link, panels=NULL, imppoly){

  if(is.null(panels)){
    panels<-1:nrow(data)
  }else{
    panels<-data[,panels]
  }

  if(length(changecoef.link)>1){
    overallchange=TRUE
    changecoef.link.overall=changecoef.link[2]
    changecoef.link=changecoef.link[1]
  }else(overallchange=FALSE)

  impactcells<-ifelse(inout(data[,names(imppoly)], imppoly,quiet = T), 0, 1)

  dat2<-rbind(data.frame(data, eventphase=0, panels=panels, impcells=0), data.frame(data,eventphase=1, panels=max(panels)+panels, impcells=impactcells))

  # make index column to keep track of data order
  dat2$index<-1:nrow(dat2)

  modmat<-cbind(rbind(model.matrix(model), model.matrix(model)), model.matrix(rnorm(nrow(dat2),0,1) ~ eventphase, data=dat2)[,-1])
  id1<-which(dat2$eventphase==1)
  id0<-which(dat2$eventphase==0)
  initcoeff<-coef(model)
  initcoeff<-as.vector(c(initcoeff, changecoef.link))
  impcoefID<-length(initcoeff)
  # find sum of animals pre impact
  ysum= sum((modmat[id0,]%*%initcoeff))
  xsub_sum<- sum((modmat[id1,] %*% initcoeff))
  beta_x <- (ysum - xsub_sum)/sum((dat2$impcells[id1]))
  #
  b1<- c(initcoeff, beta_x)

  if(overallchange){
    xx2<- cbind(modmat, dat2$impcells, dat2$eventphase)
    b1<- c(b1, changecoef.link.overall)
    dat2$truth<-as.numeric(exp((xx2%*%b1)))
  }else{
    # new model matrix with inout column added
    xx2<- cbind(modmat, dat2$impcells)
    dat2$truth<-as.numeric(exp((xx2%*%b1)))
  }
  return(dat2)
}

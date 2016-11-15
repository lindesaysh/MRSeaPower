#'
#' Function to impose an change impact effect
#'
#' @param pct.change percentage change to occur. e.g. 50\% gives 50\% site wide decline.  If a redistribution effect (one of \code{eventsite.bnd} or \code{noneventcells} must be specified) then the change refers only to the event site cells.  For an additional site-wide change when using a re-distribution, \code{pct.change} is a vector with the first number the event site re-distribution and the second, the site wide change.
#' @param model A glm or gamMRSea model object.
#' @param data Data frame used to fit \code{model}
#' @param panels Character vector denoting the column of \code{data} containing the panel structure.
#' @param eventsite.bnd A data frame containing the coordinates of a polygon defining the region of the defined event The variable names must match the coordinate system in the data.
#' @param noneventcells logical (0/1) vector of length the same as the data indicating which cells the event did not occur (1).
#'
#' @return Returns a data frame twice the size of the original input with additional columns for panel id, eventphase and truth.  The truth column represents the input data for the first half (pre-event; eventphase==0) and the second half has the change imposed (post-event; eventphase==1).
#'
#' @author Lindesay Scott-Hayward, University of St Andrews
#'
#' @export
#'


genChangeData<-function(pct.change, model, data, panels=NULL, eventsite.bnd=NULL, noneventcells=NULL){

  if(is.null(eventsite.bnd) & is.null(noneventcells)){type='oc'
  }else{type='re'}

  if(is.null(panels)){
    panels<-1:nrow(data)
  }else{
    panels<-data[,panels]
  }

  if(type=='re' & is.null(noneventcells)){
    require(splancs)
    noneventcells<-ifelse(inout(data[,names(eventsite.bnd)], eventsite.bnd,quiet = T), 0, 1)
  }
  
  if(length(pct.change)>1){
    overallchange.re=TRUE
    pct.change.oc=pct.change[2]
    pct.change=pct.change[1]
  }else(overallchange.re=FALSE)

  dat2<-rbind(data.frame(data, eventphase=0, panels=panels), data.frame(data,eventphase=1, panels=max(panels)+panels))

  fitbefore<- fitted(model)
  
  fitafter<-fitbefore*(pct.change/100)
  #fitafter<- fitbefore*0.5

  if(type=='re'){
    nonaffect<- which(noneventcells==1)
    fitafter<-fitbefore
    fitafter[-nonaffect]<- fitbefore[-nonaffect]*(pct.change/100)
    reldens<-fitted(model)[nonaffect]/sum(fitted(model)[nonaffect])
    
    newcells<- fitafter
    newcells[nonaffect]<-newcells[nonaffect]+(sum(fitbefore[-nonaffect])-sum(fitafter[-nonaffect]))*reldens

    if(overallchange.re){
      newcells<-newcells*(pct.change.oc/100)
    }

    dat2$truth<-c(fitbefore, newcells)
  }else{
    dat2$truth<-c(fitbefore, fitafter)
  }



  return(dat2)
}



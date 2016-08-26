#' Function to make model of class \code{gamMRSea}
#'
#'
#' @param model
#' @param panelid
#' @param splineParams
#'
#'
#' @example
#'

make.gamMRSea<-function(model, panelid=NULL, splineParams=NULL, varshortnames=NULL, gamMRSea=FALSE){
  newmodel<-model
  if(class(model)[1]!='gamMRSea'){
    class(newmodel)<-c('gamMRSea', class(model))
  }
  newmodel$panels<-if(is.null(panelid) & is.null(model$panels)) 1:nrow(model$data) else model$panels

  newmodel$splineParams<-splineParams
  newmodel$varshortnames<-varshortnames
  
  if(gamMRSea){
    newmodel$call[[1]]<-quote(gamMRSea)
    newmodel$call$splineParams<-quote(splineParams)
  }
  
  return(newmodel)
}

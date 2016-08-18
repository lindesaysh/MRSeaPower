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

make.gamMRSea<-function(model, panelid=NULL, splineParams=NULL, varshortnames=NULL){
  newmodel<-model
  class(newmodel)<-c('gamMRSea', class(model))
  model$panels<-if(is.null(panelid) & is.null(model$panels)) 1:nrow(model$data) else model$panels
  model$splineParams<-splineParams
  model$varshortnames<-varshortnames
  return(model)
}

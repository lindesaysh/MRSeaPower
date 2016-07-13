#'
#' Function to impose an overall change impact effect
#' 
#' @param model A glm or gamMRSea model object.
#' @param data Data frame used to fit \code{model}
#' @param changecoef.link Coefficient on the scale of the link function to determine the overall change in the population.  
#' @param panels Character vector denoting the column of \code{data} containing the panel structure.
#' 
#' 
#' @return Returns a data frame twice the size of the original input with additional columns for panel id, eventphase and truth.  The truth column represents the input data for the first half (pre-event; eventphase==0) and the second half has the change imposed (post-event; eventphase==1).
#' 
#' @author Lindesay Scott-Hayward, University of St Andrews
#' 
#' @export
#' 


genOverallchangeData<-function(changecoef.link, model, data, panels=NULL){
  
  
  if(is.null(panels)){
    panels<-1:nrow(data)  
  }else{
    panels<-data[,panels]
  }
  
  dat2<-rbind(data.frame(data, eventphase=0, panels=panels), data.frame(data,eventphase=1, panels=max(panels)+panels))
  
 modmat<-cbind(rbind(model.matrix(model), model.matrix(model)), model.matrix(response ~ eventphase, data=dat2)[,-1])
  
  initcoeff<-coef(model)
  initcoeff<-as.vector(c(initcoeff, changecoef.link))
  impcoefID<-length(initcoeff)
  dat2$truth<-exp(modmat%*%initcoeff)
  return(dat2)
}



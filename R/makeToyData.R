#' Function to create a toy dataset for package testing
#' 
#' @param n number denoting the amount of data to generate.
#' @param changecoef.link Coefficient on the scale of the link function to determine the overall change in the population. 
#' @param length.panels vector denoting the length of the panels to be created.
#' @param binary logical (default = FALSE) denoting whether to generate binary data.  If FALSE, poisson data is generated. 
#' 
#' @return Returns a data frame of columns 
#'    \code{x} x covariate
#'    \code{evph} covariate denoting pre(0)/post(1) event
#'    \code{mu} poisson response
#'    \code{mu.bin} binary response (only if \code{binary==TRUE})
#'    \code{panels} column denoting which panel each row corresponds to.  Used to subsequently add correlation to the data. 
#'    
#' @author Lindesay Scott-Hayward, University of St Andrews
#' 
#' @export
#' 


makeToyData<- function(n, changecoef.link=0, length.panels=1, binary=FALSE){
  # ~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Truth Generation for Toy ~~
  # ~~~~~~~~~~~~~~~~~~~~~~
  b0<- 1
  b1<- 0.3
  obsorder<- 1:(n/2)
  
  x<- seq(1,10, length=(n/2))
  evph<-c(0,1)
  dat<-expand.grid(x=x, evph=evph)
  
  eta<-b0+b1*dat$x + changecoef.link*dat$evph
  dat$mu = abs(exp(eta)) 
  
  if(binary) {
    dat$mu.bin<-exp(eta)/(1+exp(eta))
    #dat$mu.bin<-ifelse(dat$mu.bin>mean(dat$mu.bin), 1, 0)
  }
  
  dat$panelid<-rep((1:ceiling(n/length.panels)), each=(length.panels))[1:n]
  return(dat)
}


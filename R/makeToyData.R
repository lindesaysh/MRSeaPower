#' Function to create a toy dataset for package testing
#'
#' @param n number denoting the amount of data to generate.
#' @param changecoef.link Coefficient on the scale of the link function to determine the overall change in the population.
#' @param length.panels vector denoting the length of the panels to be created.
#' @param binary logical (default = FALSE) denoting whether to generate binary data.  If FALSE, poisson data is generated.
#' @param b0 intercept parameter
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


makeToyData<- function(n, changecoef.link=0, length.panels=1, binary=FALSE, b0=1){
  # ~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Truth Generation for Toy ~~
  # ~~~~~~~~~~~~~~~~~~~~~~
  b0<- b0
  b1<- 0.3
  obsorder<- 1:(n/2)
  X<- seq(1,10, length=5000)
  x<- sample(X, size = (n/2), replace = TRUE)
  evph<-c(0,1)
  x<-sort(x)
  dat<-expand.grid(x=x, evph=evph)

  eta<-b0+b1*dat$x + changecoef.link*dat$evph
  dat$mu = abs(exp(eta))

  if(binary) {
    dat$mu.bin<-exp(eta)/(1+exp(eta))
    #dat$mu.bin<-ifelse(dat$mu.bin>mean(dat$mu.bin), 1, 0)
  }

  dat$panels<-rep((1:ceiling(n/length.panels)), each=(length.panels))[1:n]
  return(dat)
}


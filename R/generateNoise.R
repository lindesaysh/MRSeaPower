#' Function to generate noisy data
#'
#' The function generates a random sample from poisson, overdispersed poisson, binomial and zero inflated binomial samples.
#'
#' @param n number of simulations to generate
#' @param response vector of 'true' means to genereate from
#' @param family one of \code{poisson}, \code{binomial} or \code{zibinomial}
#' @param ... Other parameters required for the family specified
#'
#' @details
#' An additional parameter for the Poisson distribution is the dispersion parameter, specified by d=
#' The additional parameters for the Binomial distribution can be found in \link{rbinom}
#' The zibinomial family requires the \code{VGAM} library to generate zero inflated binomial data. Additional parameters can be found in the help for \link{rzibinom}.
#'
#' @example
#'
#' @export
#'

generateNoise<-function(n, response, family, ...){
  simData<-matrix(NA, nrow=length(response), ncol=n)
  for(i in 1:n){
    if(family=='poisson'){
      simData[,i]<-rpois.od(length(response), response, ...)
    }
    if(family=='binomial'){
      simData[,i]<-rbinom(n=length(response), prob=response, ...)
    }
    if(family=='zibinomial'){
      simData[,i]<-rzibinom(n=length(response), prob=response, ...)
    }
  }
  return(simData)
}

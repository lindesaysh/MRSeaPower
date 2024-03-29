% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateNoise.R
\name{generateNoise}
\alias{generateNoise}
\title{Function to generate noisy data}
\usage{
generateNoise(n, response, family, ...)
}
\arguments{
\item{n}{number of simulations to generate}

\item{response}{vector of 'true' means to genereate from}

\item{family}{one of \code{poisson}, \code{binomial} or \code{zibinomial}}

\item{...}{Other parameters required for the family specified}
}
\description{
The function generates a random sample from poisson, overdispersed poisson, binomial and zero inflated binomial samples.
}
\details{
An additional parameter for the Poisson distribution is the dispersion parameter, specified by d=
The additional parameters for the Binomial distribution can be found in \link{rbinom}
The zibinomial family requires the \code{VGAM} library to generate zero inflated binomial data. Additional parameters can be found in the help for \link{rzibinom}.
}
\examples{

data(nystedA_slim)

model<-glm(response ~ depth  + as.factor(phase),
              family='poisson', data=nysted)

simData<-generateNoise(n=500, response=fitted(model), family='poisson')

}
\author{
LAS Scott-Hayward, University of St Andrews
}

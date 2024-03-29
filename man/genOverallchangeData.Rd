% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genOverallchangeData.R
\name{genOverallchangeData}
\alias{genOverallchangeData}
\title{Function to impose an overall change impact effect}
\usage{
genOverallchangeData(changecoef.link, model, data, panels = NULL)
}
\arguments{
\item{changecoef.link}{Coefficient on the scale of the link function to determine the overall change in the population.}

\item{model}{A glm or gamMRSea model object.}

\item{data}{Data frame used to fit \code{model}}

\item{panels}{Character vector denoting the column of \code{data} containing the panel structure.}
}
\value{
Returns a data frame twice the size of the original input with additional columns for panel id, eventphase and truth.  The truth column represents the input data for the first half (pre-event; eventphase==0) and the second half has the change imposed (post-event; eventphase==1).
}
\description{
Function to impose an overall change impact effect
}
\author{
Lindesay Scott-Hayward, University of St Andrews
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genRedistData.R
\name{genRedistData}
\alias{genRedistData}
\title{Function to impose a re-distribution impact effect}
\usage{
genRedistData(model, data, changecoef.link, panels = NULL,
  eventsite.bnd = NULL, noneventcells = NULL)
}
\arguments{
\item{model}{A glm or gamMRSea model object.}

\item{data}{Data frame used to fit \code{model}}

\item{changecoef.link}{Coefficient(s) on the scale of the link function to determine the change in the population.  If of length 1, then the change occurs in the impact zone, whilst maintaining overall population size.  If of length 2, then there is a redistribution of the population with an an additional overall decrease.}

\item{panels}{Character vector denoting the column of \code{data} containing the panel structure.}

\item{eventsite.bnd}{A data frame containing the coordinates of a polygon defining the region of the defined event The variable names must match the coordinate system in the data.}

\item{noneventcells}{vector of zeros and ones of length the same as data. zeros represent event cells and ones represent non-event cells.}
}
\value{
Returns a data frame twice the size of the original input with additional columns for panel id, eventphase and truth.  The truth column represents the input data for the first half (pre-event; eventphase==0) and the second half has the change imposed (post-event; eventphase==1).
}
\description{
Function to impose a re-distribution impact effect
}
\author{
Lindesay Scott-Hayward, University of St Andrews
}

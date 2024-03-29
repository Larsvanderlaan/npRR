% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npRR.R
\name{predict.npRR}
\alias{predict.npRR}
\title{Predict relative risk values at new observations.}
\usage{
\method{predict}{npRR}(object, newV)
}
\arguments{
\item{object}{A \code{\link{npRR}} object}

\item{newV}{A vector or matrix values at which to evaluate the relative risk function. The order of columns must match the argument \code{V} passed to \code{\link{npRR}} exactly.}
}
\value{
Returns a vector or matrix of predicted RR values.
}
\description{
Given an npRR object returned by \code{\link{npRR}}, this functions returns the relative risk values for given values of `V`.
}

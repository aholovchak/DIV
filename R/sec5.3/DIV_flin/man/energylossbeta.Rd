% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/energylossbeta.R
\name{energylossbeta}
\alias{energylossbeta}
\title{Energy Loss Calculation with Beta Scaling}
\usage{
energylossbeta(x0, x, xp, beta, verbose = FALSE)
}
\arguments{
\item{x0}{A tensor representing the target values.}

\item{x}{A tensor representing the model's stochastic predictions.}

\item{xp}{A tensor representing another draw of the model's stochastic predictions.}

\item{beta}{A numeric value for scaling the energy loss.}

\item{verbose}{A boolean indicating whether to return prediction loss s1 = E(||x0-x||) and variance loss s2 = E(||x-xp||).}
}
\value{
A scalar representing the calculated energy loss.
}
\description{
This function calculates the energy loss for given tensors. The loss is calculated
as the mean of the L2 norms between \code{x0} and \code{x} and between \code{x0} and \code{xp}, each raised to the power of \code{beta},
subtracted by half the mean of the L2 norm between \code{x} and \code{xp}, also raised to the power of \code{beta}.
}
\keyword{internal}

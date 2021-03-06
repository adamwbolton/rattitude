% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{calibration}
\alias{calibration}
\title{Calibration model for}
\usage{
calibration(V, G)
}
\arguments{
\item{V}{Matrix (n x 3) of observed response where each row denotes v_i = (v_{ix}, v_{iy}, v_{iz})}

\item{G}{Matrix (n x 3) of fixed gravity vectors where each row denotes g_i = (g_{ix}, g_{iy}, g_{iz}). Each row should have unit norm.}
}
\value{
List of estimated coefficients from linear model: BetaHat (3x3) matrix of first-order effects, Beta0Hat (1x3) vector of intercepts,
Ghat (n x3) matrix of estimated gravity vectors, sd is estimated standard error of responses.
}
\description{
Calibration model for
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdaplot.R
\name{mdaplot.formatValues}
\alias{mdaplot.formatValues}
\title{Format vector with numeric values}
\usage{
mdaplot.formatValues(data, round.only = FALSE, digits = 3)
}
\arguments{
\item{data}{vector or matrix with values}

\item{round.only}{logical, do formatting or only round the values}

\item{digits}{how many significant digits take into account}
}
\value{
matrix with formatted values
}
\description{
Format vector with values, so only significant decimal numbers are left.
}
\details{
Function takes into accound difference between values and the values themselves.
}

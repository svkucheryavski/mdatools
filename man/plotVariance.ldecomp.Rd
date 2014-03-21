\name{plotVariance.ldecomp}
\alias{plotVariance.ldecomp}
\title{Explained variance plot for linear decomposition}
\usage{
\method{plotVariance}{ldecomp}(obj, type = "b", main = "Variance",
  xlab = "Components", ylab = "Explained variance, \%", show.labels = F,
  ...)
}
\arguments{
  \item{obj}{object of \code{ldecomp} class.}

  \item{type}{type of the plot}

  \item{main}{main title for the plot}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.labels}{logical, show or not labels for the
  plot objects}

  \item{...}{most of graphical parameters from
  \code{\link{mdaplot}} function can be used.}
}
\description{
Shows a plot with explained variance values vs. number of
components.
}


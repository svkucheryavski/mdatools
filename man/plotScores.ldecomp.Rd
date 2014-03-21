\name{plotScores.ldecomp}
\alias{plotScores.ldecomp}
\title{Scores plot for linear decomposition}
\usage{
\method{plotScores}{ldecomp}(obj, comp = c(1, 2), main = "Scores",
  show.labels = F, show.axes = F, ...)
}
\arguments{
  \item{obj}{object of \code{ldecomp} class.}

  \item{comp}{which components to show the plot for (can be
  one value or vector with two values).}

  \item{main}{main title for the plot}

  \item{show.labels}{logical, show or not labels for the
  plot objects}

  \item{show.axes}{logical, show or not a axes lines
  crossing origin (0,0)}

  \item{...}{most of graphical parameters from
  \code{\link{mdaplot}} function can be used.}
}
\description{
Shows a plot with scores values for data objects.
}


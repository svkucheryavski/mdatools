\name{plot.plsdares}
\alias{plot.plsdares}
\title{Overview plot for PLS-DA results}
\usage{
\method{plot}{plsdares}(x, nc = NULL, ncomp = NULL, show.labels = F, ...)
}
\arguments{
  \item{x}{PLS-DA results (object of class
  \code{plsdares})}

  \item{ncomp}{how many components to use (if NULL - user
  selected optimal value will be used)}

  \item{nc}{which class to show the summary for (if NULL,
  will be shown for all)}

  \item{show.labels}{logical, show or not labels for the
  plot objects}

  \item{...}{other arguments}
}
\description{
Shows a set of plots (x residuals, y variance,
classification performance and predictions) for PLS-DA
results.
}
\details{
See examples in help for \code{\link{pls}} function.
}


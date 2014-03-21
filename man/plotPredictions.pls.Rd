\name{plotPredictions.pls}
\alias{plotPredictions.pls}
\title{Predictions plot for PLS}
\usage{
\method{plotPredictions}{pls}(obj, ncomp = NULL, ny = 1, main = NULL,
  xlab = NULL, ylab = NULL, legend.position = "topleft", show.line = T,
  show.labels = F, show.legend = T, colmap = "default", col = NULL, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{ncomp}{how many components to use (if NULL - user
  selected optimal value will be used)}

  \item{ny}{number of response variable to make the plot
  for (if y is multivariate)}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{legend.position}{position of legend on the plot (if
  shown)}

  \item{show.line}{logical, show or not line fit for the
  plot points}

  \item{show.labels}{logical, show or not labels for the
  plot elements}

  \item{show.legend}{logical, show or not legend for the
  plot}

  \item{colmap}{a colormap to use for coloring the plot
  items}

  \item{col}{a vector with color values for plot items}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with predicted vs. reference (measured) y values
for selected components.
}
\details{
See examples in help for \code{\link{pls}} function.
}


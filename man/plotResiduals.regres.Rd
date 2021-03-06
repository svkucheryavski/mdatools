% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regres.R
\name{plotResiduals.regres}
\alias{plotResiduals.regres}
\title{Residuals plot for regression results}
\usage{
\method{plotResiduals}{regres}(
  obj,
  ny = 1,
  ncomp = obj$ncomp.selected,
  show.lines = c(NA, 0),
  show.plot = TRUE,
  ...
)
}
\arguments{
\item{obj}{regression results (object of class \code{regres})}

\item{ny}{number of predictor to show the plot for (if y is multivariate)}

\item{ncomp}{complexity of model (e.g. number of components) to show the plot for}

\item{show.lines}{allows to show the horisontal line at y = 0}

\item{show.plot}{logical, show plot or just return plot data}

\item{...}{other plot parameters (see \code{mdaplot} for details)}
}
\description{
Shows plot with Y residuals (difference between predicted and reference values) for selected
response variable and complexity (number of components).
}

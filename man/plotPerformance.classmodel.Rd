\name{plotPerformance.classmodel}
\alias{plotPerformance.classmodel}
\title{Performance plot for classification model}
\usage{
\method{plotPerformance}{classmodel}(obj, nc = NULL, param = "specificity",
  type = "h", legend = NULL, main = NULL, xlab = "Components",
  ylab = "", ylim = c(0, 1.15), ...)
}
\arguments{
  \item{obj}{classification model (object of class
  \code{plsda}, \code{simca}, etc.).}

  \item{nc}{if there are several classes, which class to
  make the plot for (NULL - summary for all classes).}

  \item{param}{which parameter to make the plot for
  (\code{'specificity'}, \code{'sensitivity'}, or
  \code{'misclassified'})}

  \item{type}{type of the plot}

  \item{legend}{vector with legend items}

  \item{main}{main title for the plot}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{ylim}{vector with two values - limits for y axis}

  \item{...}{most of the graphical parameters from
  \code{\link{mdaplotg}} function can be used.}
}
\description{
Makes a plot with sensitivity values vs. model complexity
(e.g. number of components)
}


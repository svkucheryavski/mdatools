\name{plot.regcoeffs}
\alias{plot.regcoeffs}
\title{Regression coefficients plot}
\usage{
\method{plot}{regcoeffs}(x, ncomp = 1, ny = 1, type = NULL,
  main = "Regression coefficients", xlab = "Variables",
  ylab = "Coefficients", show.line = T, ...)
}
\arguments{
  \item{x}{regression coefficients object (class
  \code{regcoeffs})}

  \item{ncomp}{number of components to return the
  coefficients for}

  \item{ny}{number of response variable to return the
  coefficients for}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.line}{logical, show or not line for 0 value}

  \item{...}{other arguments}
}
\description{
Shows plot with regression coefficient values for every
predictor variable (x)
}


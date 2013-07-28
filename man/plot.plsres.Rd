\name{plot.plsres}
\alias{plot.plsres}

\title{
Plot overview of PLS results
}

\description{
\code{plot.plsres} is used to show a set of plots for inspecting of PLS results.
}

\usage{
plot(obj, ncomp = NULL, ny = 1, show.labels = F)
}

\arguments{
  \item{obj }{PLS results (object of class \code{plsres}).}
  \item{ncomp }{number of components to show the plot for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to show the plot for.}
  \item{show.labels}{logical, show or not object labels for the data points.}
}

\details{
By default function shows X scores, Y variance, RMSE and predictions plots. For cross-validated results, where X scores values are not available, X variance is shown instead.
}
\examples{
## Examples of using plot function for PLS results

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, cv = 4)
model = selectCompNum(model, 3)

# for calibration results
plot(model$calres, ny = 2, ncomp = 2, show.labels = T)

# for cross-validation results
plot(model$cvres, ny = 2, ncomp = 2, show.labels = T)

}

\name{plotXYLoadings.pls}
\alias{plotXYLoadings.pls}

\title{
XY loadings plot for PLS model
}

\description{
\code{plotXYLoadings.pls} is used to make a plot with loading values for PLS decomposition of predictors and responses.
}

\usage{
plotXYLoadings(model, comp = c(1, 2), show.axes = F, ...)
}

\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{comp }{vector with two values - number of components to make the plot for. }  
  \item{show.axes }{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\examples{
## Examples of using XY loadings plot

library(mdatools)
data(People)

# will predict Shoesize using the other variables
xdata = people[, -4]
ydata = people[, 4]

model = pls(xdata, ydata, ncomp = 8, scale = T)

# show cumulative variance plot to see importance of components
# as well as several XY loading plots
par(mfrow = c(2, 2))
plotYCumVariance(model)
plotXYLoadings(model, show.labels = T, show.axes = T)
plotXYLoadings(model, c(2, 3), show.labels = T)
plotXYLoadings(model, show.labels = T, colmap = c('green', 'orange'))
par(mfrow = c(1, 1))
}
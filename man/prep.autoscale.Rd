\name{prep.autoscale}
\alias{prep.autoscale}

\title{
Data autoscaling
}

\description{
\code{prep.autoscale} is used to mean center and/or standardize data values.
}

\usage{
prep.autoscale(data, center = T, scale = F)
}

\arguments{
  \item{data}{a matrix with data values}
  \item{center}{either logical value (center or not) or vector with numbers used to center the values.}
  \item{scale}{either logical value (standardize or not) or vector with numbers used to scale the values.}
}
\details{
If \code{center} and/or \code{scale} parameter is a vector it should contain as many values as many variables (columns) in the data matrix.

The function assign two attributes to the preprocessed matrix: \code{prep:center} and \code{prep:scale} with values, used to center and scale the data.
}

\value{
\item{pdata }{data matrix after preprocessing.}
}

\examples{
## Examples of using autoscaling

library(mdatools)
data(People)

# make autoscaling
odata = people[, 1:8]
pdata = prep.autoscale(odata, center = T, scale = T)

# show standard deviation of each variable
par(mfrow = c(2, 1))
mdaplot(cbind(1:ncol(odata), apply(odata, 2, sd)), type = 'h', 
   ylab = 'Standard deviation', main = 'Original data', show.labels = T)
mdaplot(cbind(1:ncol(pdata), apply(pdata, 2, sd)), type = 'h', 
   ylab = 'Standard deviation', main = 'After autoscaling', show.labels = T)
par(mfrow = c(1, 1))

# show statistics and attributes
show(apply(odata, 2, sd))
show(apply(odata, 2, mean))
show(apply(pdata, 2, sd))
show(apply(pdata, 2, mean))
show(attr(pdata, 'prep:center'))
show(attr(pdata, 'prep:scale'))

}

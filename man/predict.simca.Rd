\name{predict.simca}
\alias{predict.simca}

\title{
Predict method for SIMCA model
}

\description{
\code{predict.simca} is used to classify a new data with SICMA model.
}

\usage{
predict(model, data, c.ref = NULL)
}

\arguments{
  \item{model}{an object of class \code{simca} (a SIMCA model).}
  \item{data}{a numerical matrix with data.}
  \item{c.ref}{a logical vector with class membership for the data objects.}
}

\value{
Returns object of class \code{simcares} with results of SIMCA one-class classification of the data. See \code{\link{simcares}} examples for details.
}

\seealso{
\code{\link{simca}} to make SIMCA model, \code{\link{simcares}} to store and explore SIMCA one-class classification results.
}

\examples{
## make a SIMCA model for Iris setosa class and apply it to a new data set
library(mdatools)

data = iris[, 1:4]
class = iris[, 5]

# take first 30 objects of setosa as calibration set 
se = data[1:30, ]

# use the next 40 as test set and make logical reference values
test = data[31:70, ]
testclass = class[31:70] == 'setosa'

# make SIMCA model and apply to test set
model = simca(se, 'Se')
model = selectCompNum(model, 1)
pred = predict(model, test, c.ref = testclass)

# show infromation and summary
print(pred)
summary(pred)

# show plots for predictions
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plotPredictions(pred, show.labels = T)
plotResiduals(pred)
plotPerformance(pred, show.labels = T, legend.position = 'bottomright')
layout(1, 1, 1)
}

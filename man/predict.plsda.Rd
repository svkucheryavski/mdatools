\name{predict.plsda}
\alias{predict.plsda}

\title{
Predict method for PLS-DA model
}

\description{
predict.plsda is used to apply a PLS-DA model to a new data set.
}

\usage{
predict(model, x, c.ref = NULL)
}

\arguments{
  \item{model}{a PLS-DA model (object of class \code{plsds})}
  \item{x}{a matrix with predictor values.}
  \item{c.ref}{a vector with reference class values for objects (optional).}
}

\details{
The method applies a PLS-DA model to a new data. If response values are not provided, no classification statistics will be available.
}

\value{Returns an object of \code{\link{plsdares}} class.}


\examples{
## make a PLS-DA model a apply it to a new data

library(mdatools)

# calibration set
x.cal = iris[seq(1, nrow(iris), 2), 1:4]  
c.cal = iris[seq(1, nrow(iris), 2), 5]

# new data set
x.new = iris[seq(2, nrow(iris), 2), 1:4] 
c.new = iris[seq(2, nrow(iris), 2), 5] 

model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
model = selectCompNum(model, 2)

res = predict(model, testset, testset.c)

# show summary and basic plots
summary(res)
plot(res)

# show summary and basic plots for particular class and ncomp
summary(res, nc = 3)
plot(res, nc = 3, ncomp = 1)


# show performance plots
par(mfrow = c(2, 2))
plotPerformance(res)
plotPerformance(res, nc = 3)
plotSensitivity(res)
plotSensitivity(res, nc = 3, type = 'b')
par(mfrow = c(1, 1))

# show predictions
par(mfrow = c(2, 2))
plotPredictions(res)
plotPredictions(res, nc = 3, ncomp = 3)
plotPredictions.plsres(res)
plotPredictions.plsres(res, ny = 3, ncomp = 3)
par(mfrow = c(1, 1))

}


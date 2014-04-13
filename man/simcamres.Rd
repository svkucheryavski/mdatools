\name{simcamres}
\alias{simcamres}

\title{
Results of SIMCA multiclass classification
}

\description{
\code{simcamres} is used to store results for SIMCA multiclass classification.
}

\usage{
simcamres(cres, T2, Q2, T2lim, Q2lim)
}

\arguments{
  \item{cres }{results of classification (class \code{classres}).}
  \item{T2 }{matrix with T2 values for each object and class.}
  \item{Q2 }{matrix with Q2 values for each object and class.}
  \item{T2lim }{vector with T2 statistical limits for each class.}
  \item{Q2lim }{vector with Q2 statistical limits for each class.}
}

\details{
Class \code{simcamres} inherits all properties and methods of class \code{\link{classres}}, plus 
store values necessary to visualise prediction decisions (e.g. Cooman's plot or Residuals plot).

In cotrast to \code{simcares} here only values for optimal (selected) number of components in 
each individual SIMCA models are presented.

There is no need to create a \code{simcamres} object manually, it is created automatically when make a SIMCAM model 
(see \code{\link{simcam}}) or apply the model to a new data (see \code{\link{predict.simcam}}). The object can be used 
to show summary and plots for the results.
}

\value{
Returns an object (list) of class \code{simcamres} with the same fields as \code{\link{classres}} plus extra fields for 
Q2 and T2 values and limits:

\item{c.pred}{predicted class values.}
\item{c.ref}{reference (true) class values if provided.}
\item{T2}{matrix with T2 values for each object and class.}
\item{Q2}{matrix with Q2 values for each object and class.}
\item{T2lim}{vector with T2 statistical limits for each class.}
\item{Q2lim}{vector with Q2 statistical limits for each class.}

The following fields are available only if reference values were provided.
\item{tp}{number of true positives.}
\item{fp}{nmber of false positives.}
\item{fn}{number of false negatives.}
\item{specificity}{specificity of predictions.}
\item{sensitivity}{sensitivity of predictions.}
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}


\seealso{
Methods for \code{simcamres} objects:
  \tabular{ll}{
      \code{print.simcamres} \tab shows information about the object.\cr
      \code{summary.simcamres} \tab shows statistics for results of classification.\cr
      \code{\link{plotResiduals.simcamres}} \tab makes Q2 vs. T2 residuals plot.\cr
      \code{\link{plotCooman.simcamres}} \tab makes Cooman's plot.\cr
   }
   
Methods, inherited from \code{\link{classres}} class:
  \tabular{ll}{
      \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
      \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
   }

   Check also \code{\link{simcam}}.    
}

\examples{
## make a multiclass SIMCA model for Iris data and apply to test set
library(mdatools)

# split data 
caldata = iris[seq(1, nrow(iris), 2), 1:4]
se = caldata[1:25, ]
ve = caldata[26:50, ]
vi = caldata[51:75, ]

testdata = iris[seq(2, nrow(iris), 2), 1:4]
testdata.cref = rbind(matrix('Se', 25, 1), matrix('Vi', 25, 1), matrix('Ve', 25, 1))

# create individual models
semodel = simca(se, classname = 'Se')
semodel = selectCompNum(semodel, 1)

vimodel = simca(vi, classname = 'Vi')
vimodel = selectCompNum(vimodel, 1)

vemodel = simca(ve, classname = 'Ve')
vemodel = selectCompNum(vemodel, 1)

# combine models into SIMCAM object, show statistics 
model = simcam(list(semodel, vimodel, vemodel), info = 'Iris data')
res = predict(model, testdata, testdata.cref)
summary(res)

# show predicted values
showPredictions(res)

# plot predictions
par(mfrow = c(2, 2))
plotPredictions(res)
plotPredictions(res, nc = 1)
plotPredictions(res, nc = c(1, 2))
plotPredictions(res, show.labels = TRUE)
par(mfrow = c(1, 1))

# show residuals and Cooman's plot

par(mfrow = c(2, 2))
plotCooman(res)
plotCooman(res, nc = c(1, 3))
plotResiduals(res)
plotResiduals(res, nc = 3)
par(mfrow = c(1, 1))

}


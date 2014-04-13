\name{pca}
\alias{pca}

\title{
Principal Component Analysis
}

\description{
\code{pca} is used to build and explore a principal component analysis (PCA) model.
}

\usage{
pca(x, ncomp = 15, center = T, scale = F, cv = NULL, x.test = NULL, 
   alpha = 0.05, method = 'svd', info = '')
}

\arguments{
  \item{x}{a numerical matrix with calibration data.}
  \item{ncomp}{maximum number of components to calculate.}
  \item{center}{logical, do mean centering of data or not.}
  \item{scale}{logical, do sdandardization of data or not.}
  \item{cv}{number of segments for random cross-validation (1 for full cross-validation).}
  \item{x.test}{a numerical matrix with test data.}
  \item{alpha}{significance level for calculating limit for Q2 residuals.}
  \item{method}{method to compute principal components.}
  \item{info}{a short text line with model description.}
}

\details{
So far only SVD (Singular Value Decompisition) method is available, more coming soon. 

By default \code{pca} uses number of components (\code{ncomp}) as a minimum of number of objects - 1, number of variables and default or provided value. Besides this, there is also a parameter for selecting an optimal number of components (\code{ncomp.selected}). The optimal number of components is used to build a residuals plot (with Q2 residuals vs. Hotelling T2 values), calculate confidence limits for Q2 residuals, as well as for SIMCA classification. 

If data contains missing values (NA) the \code{pca} will use an iterative algorithm to fit the values with most probable ones. The algorithm is implemented in a function \code{\link{pca.mvreplace}}. The same center and scale options will be used. You can also do this step manually before calling \code{pca} and play with extra options.

}

\value{
Returns an object of \code{pca} class with following fields:
\item{ncomp }{number of components included to the model.} 
\item{ncomp.selected }{selected (optimal) number of components.} 
\item{loadings }{matrix with loading values (nvar x ncomp).} 
\item{eigenvals }{vector with eigenvalues for all existent components.} 
\item{expvar }{vector with explained variance for each component (in percent).} 
\item{cumexpvar }{vector with cumulative explained variance for each component (in percent).} 
\item{T2lim }{statistical limit for T2 distance.} 
\item{Q2lim }{statistical limit for Q2 distance.} 
\item{info }{information about the model, provided by user when build the model.} 
\item{calres }{an object of class \code{\link{pcares}} with PCA results for a calibration data.} 
\item{testres }{an object of class \code{\link{pcares}} with PCA results for a test data, if it was provided.} 
\item{cvres }{an object of class \code{\link{pcares}} with PCA results for cross-validation, if this option was chosen.} 
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}


\seealso{
Methods for \code{pca} objects:
  \tabular{ll}{
      \code{plot.pca} \tab makes an overview of PCA model with four plots.\cr
      \code{summary.pca} \tab shows some statistics for the model.\cr
      \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
      \code{\link{predict.pca}} \tab applies PCA model to a new data.\cr
      \code{\link{plotScores.pca}} \tab shows scores plot.\cr
      \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
      \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
      \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
      \code{\link{plotResiduals.pca}} \tab shows Q2 vs. T2 residuals plot.\cr
       }
   Most of the methods for plotting data are also available for PCA results 
   (\code{\link{pcares}}) objects.    

Other methods implemented in \code{pca}:
  \tabular{ll}{
      \code{\link{pca.mvreplace}} \tab replaces missing values in a data matrix with approximated using iterative PCA decomposition.\cr
       }
}
\examples{
### Examples for PCA class

## 1. Make PCA model for People data with autoscaling
## and full cross-validation

library(mdatools)

data(people)
model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
model = selectCompNum(model, 4)
summary(model)
plot(model, show.labels = TRUE)

## 2. Add missing values, make a new model and show plots
peoplemv = people
peoplemv[2, 7] = NA
peoplemv[6, 2] = NA
peoplemv[10, 4] = NA
peoplemv[22, 12] = NA

modelmv = pca(peoplemv, scale = TRUE, info = 'Model with missing values')
modelmv = selectCompNum(modelmv, 4)
summary(modelmv)
plot(modelmv, show.labels = TRUE)

## 3. Show scores and loadings plots for the model
par(mfrow = c(2, 2))
plotScores(model, comp = c(1, 3), show.labels = TRUE)
plotScores(model, comp = 2, type = 'h', show.labels = TRUE)
plotLoadings(model, comp = c(1, 3), show.labels = TRUE)
plotLoadings(model, comp = c(1, 2), type = 'h', show.labels = TRUE)
par(mfrow = c(1, 1))

## 4. Show residuals and variance plots for the model
par(mfrow = c(2, 2))
plotVariance(model, type = 'h')
plotCumVariance(model, show.labels = TRUE, legend.position = 'bottomright')
plotResiduals(model, show.labels = TRUE)
plotResiduals(model, ncomp = 2, show.labels = TRUE)
par(mfrow = c(1, 1))

}

\keyword{ ~pca }

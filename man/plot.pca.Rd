\name{plot.pca}
\alias{plot.pca}

\title{
Plot overview of PCA model
}

\description{
Shows four plots for a \code{pca} object: scores, loadings, residuals (Q2 vs. T2) and cumulative explained variance. 
}

\usage{
plot(model, comp = c(1, 2), show.labels = F, show.legend = T)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{\code{pca} object (PCA model)}
  \item{comp}{x and y components for scores and loadings plots}
  \item{show.labels}{logical, show labels for objects on scores and residuals plots or not}
  \item{show.legend}{logical, show legend or not (legend is available if test or cross-validation used)}
}

\details{
If model was validated with test set or cross-validation the results of calibration and validation are shown on scores, residuals and explained variance plots as several groups of points (or lines) separated with colors. 
}

\author{
Sergey Kucheryavskiy (svkucheryavskiy@gmail.com)
}

\seealso{
Main plots for pca model and results: \code{\link{plotScores.pca}}, \code{\link{plotLoadings.pca}}, \code{\link{plotResiduals.pca}}, \code{\link{plotVariance.pca}}, \code{\link{plotCumVariance.pca}}, \code{\link{plotEigenvalues.pca}}.
}

\examples{
## make a PCA model for people data and show main plots
library(mdatools)

data(People)
model = pca(people, scale = T, cv = 1)
model = pca.selectCompNum(model, 5)
plot(model, c(1, 3), show.labels = T)

}

\keyword{ ~pca }

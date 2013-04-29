\name{plotResiduals.pca}
\alias{plotResiduals.pca}

\title{
Residuals plot for PCA model
}
\description{
\code{plotResiduals.pca} is used to show plot for Hotelling T2 values vs. Q2 residuals for PCA model
}
\usage{
plotResiduals(obj, ncomp = NULL, show.labels = F, show.legend = T, show.limits = T)
}
\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{comp}{a number of component to show the residuals for.} 
  \item{show.labels}{logical, show or not labels for the plot points.}
  \item{show.legend}{logical, show or not a legend if test or cross-validated results are also presented in the plot.}
  \item{show.limits}{logical, show or not significance limits for T2 and Q2 values.}
}
\details{
If number of components is not provided, the selected value (\code{ncomp.selected}) is used. If number of components is not equal to the selected value, limits are not shown even if \code{show.limits} is set to \code{TRUE}.
}

\author{
Sergey Kucheryvskiy (svkucheryavski@gmail.com)
}

\examples{
## Examples of using residuals plot

library(mdatools)
data(People)

# make a PCA model and show residuals plots
model = pca(people, scale = T)
model = selectCompNum(model, 5)

par(mfrow = c(2, 2))
plotResiduals(model)
plotResiduals(model, show.labels = T)
plotResiduals(model, 3)
plotResiduals(model, 4, show.labels = T)
par(mfrow = c(1, 1))

# make a model with full cross-validation and show residuals plots
model = pca(people, scale = T, cv = 1)
model = selectCompNum(model, 4)

par(mfrow = c(2, 2))
plotResiduals(model, show.labels = T)
plotResiduals(model, show.labels = T, show.limits = F)
plotResiduals(model, 3, show.legend = F)
plotResiduals(model, 4, show.labels = T)
par(mfrow = c(1, 1))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

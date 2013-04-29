\name{selectCompNum.pca}
\alias{selectCompNum.pca}

\title{
Set number of optimal components
}

\description{
\code{selectCompNum} is used to select number of optimal components in a model.
}

\usage{
selectCompNum(obj, ncomp)
}

\arguments{
  \item{obj}{an object of  \code{pca} class (PCA model).}
  \item{ncomp}{number of components to select.}
}

\value{
Returns the same model as provided as an argument, but where the parameter \code{ncomp.selected} is set to \code{ncomp} and all distances and limits are recalculated. If number is not set, number of all computed components will be considered as optimal.
}

\author{
Sergey Kucheryavski (svkucheryavski@gmail.com)
}


\examples{
## make PCA model and show residuals for different number of selected components

library(mdatools)
data(People)

model = pca(people, ncomp = 8, scale = T)

par(mfrow = c(2, 2))
plotResiduals(model, show.labels = T)

model = selectCompNum(model, 5)
plotResiduals(model, show.labels = T)

model = selectCompNum(model, 3)
plotResiduals(model, show.labels = T)

model = selectCompNum(model, 2)
plotResiduals(model, show.labels = T)

}

\keyword{ ~kwd1 }

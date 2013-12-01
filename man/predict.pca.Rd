\name{predict.pca}
\alias{predict.pca}

\title{
Predict method for PCA model
}

\description{
Projects data to principal component space defined by a \code{pca} object.
}

\usage{
predict(obj, data)
}

\arguments{
  \item{obj}{an object of class \code{pca} (a PCA model).}
  \item{data}{a numerical matrix with data.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns object of class \code{pcares} with results of PCA decomposition of the data.
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{pca}} to make PCA model, \code{\link{pcares}} to store and explore PCA results.
}

\examples{
library(mdatools)
data(People)
   
# take every fourth row as a new data set
tdata = people[seq(1, nrow(people), 4), ]

# take the rest as a calibration set
cdata = people[-seq(1, nrow(people), 4), ]

# make a PCA model
model = pca(cdata, scale = T)
model = selectCompNum(model, 4)
   
# apply model to the new data and show scores and residuals 
res = predict(model, tdata)
par(mfrow = c(1, 2))
plotScores(res, show.labels = T)
plotResiduals(res, show.labels = T)
par(mfrow = c(1, 1))
}

\keyword{ ~pca }

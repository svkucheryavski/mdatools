\name{pca.mvreplace}
\alias{pca.mvreplace}

\title{
Replace missing values in data
}

\description{
\code{pca.mvreplace} is used to replace missing values in a data matrix with approximated by iterative PCA decomposition.
}

\usage{
pca.mvreplace(x, center = T, scale = F, maxncomp = 7, expvarlim = 0.95, 
               covlim = 10^-6, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix with data, containing missing values.}
  \item{center}{logical, do centering of data values or not.}
  \item{scale}{logical, do standardization of data values or not.}
  \item{maxncomp}{maximum number of components in PCA model.}
  \item{expvarlim}{minimum amount of variance, explained by chosen components (used for selection of optimal number of components in PCA models).}
  \item{covlim}{convergence criterion.}
  \item{maxiter}{maximum number of iterations if convergence criterion is not met.}
}

\details{
The function uses iterative PCA modeling of the data to approximate and impute missing values.  The result is most optimal for data sets with low or moderate level of noise and with number of missing values less than 10\% for small dataset and up to 20\% for large data.
}

\value{
Returns the same matrix \code{x} where missing values are replaced with approximated.
}

\references{
Philip R.C. Nelson, Paul A. Taylor, John F. MacGregor. Missing data methods in PCA and PLS: Score calculations with incomplete observations. Chemometrics and Intelligent Laboratory Systems, 35 (1), 1996.
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## A very simple example of imputing missing values in a data with no noise

# generate a matrix with values
s = 1:6
odata = cbind(s, 2*s, 4*s)

# make a matrix with missing values
mdata = odata
mdata[5, 2] = mdata[2, 3] = NA

# replace missing values with approximated
rdata = pca.mvreplace(mdata, scale = TRUE)

# show all matrices together
show(cbind(odata, mdata, round(rdata, 2)))
}

\keyword{ ~kwd1 }
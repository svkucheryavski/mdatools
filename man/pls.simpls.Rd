\name{pls.simpls}
\alias{pls.simpls}

\title{
SIMPLS algorithm for PLS regression
}

\description{
pls.simpls is an implementation of SIMPLS algorithm for calibration of PLS regression models.}

\usage{
pls.simpls(x, y, ncomp)
}

\arguments{
  \item{x}{matrix with predictor values, centered and preprocessed.}
  \item{y}{matrix with response values, centered and preprocessed.}
  \item{ncomp}{number of components to calculate the results for.}
}

\value{
Returns a list with following fields (calibration results):
\item{coeffs }{matrix with regression coefficients for each set of components (nvar x ncomp).} 
\item{weights }{matrix with PLS weights.} 
\item{xloadings }{matrix with loading values for x decomposition (nvar x ncomp).} 
\item{yloadings }{matrix with loading values for y decomposition (nvar x ncomp).} 
\item{xscores }{matrix with score values for x decomposition (nvar x ncomp).} 
\item{yscores }{matrix with score values for y decomposition (nvar x ncomp).} 
}

\references{
1. S. de Jong, SIMPLS: an alternative approach to partial least squares regression, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251–263.
}

\author{
Implemented by Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

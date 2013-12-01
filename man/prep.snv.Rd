\name{prep.snv}
\alias{prep.snv}

\title{
Standard Normal Variate preprocessing
}

\description{
}

\usage{
prep.snv(data)
}

\arguments{
  \item{data}{matrix with data values.}
}

\details{
Standard Normal Variate (SNV) is a simple preprocessing for removing baseline offset and, partly, slope of spectra data. For each spectrum (row of \code{data} matrix) an average value is subtracted and then it is divided to standard deviation. 
}

\value{
\item{pdata }{preprocessed data values.}
}

\references{
1. R. J. Barnes, M. S. Dhanoa, S. J. Lister. Standard Normal Variate Transformation and De-trending of NearInfrared Diffuse Reflectance Spectra. Applied Spectroscopy, 43 (5), pp. 772-777 (1989)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
## Examples of using SNV preprocessing

library(mdatools)
data(Simdata)

ospectra = spectra.c
pspectra = prep.snv(ospectra)

par(mfrow = c(2, 1))
mdaplot(cbind(wavelength, t(ospectra)), type = 'l', main = 'Original spectra')
mdaplot(cbind(wavelength, t(pspectra)), type = 'l', main = 'SNV preprocessed spectra')
par(mfrow = c(2, 2))


}

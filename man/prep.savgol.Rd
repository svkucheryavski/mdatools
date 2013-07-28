\name{prep.savgol}
\alias{prep.savgol}

\title{
Savitzky-Golay filter
}

\description{
\code{prep.savgol} is used to smooth and take a derivative of data values.
}
\usage{
prep.savgol(data, width = 3, porder = 1, dorder = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix with data values.}
  \item{width}{width of the filter window.}
  \item{porder}{polynomial order.}
  \item{dorder}{derivative order (must be ≥ \code{porder})}
}

\value{
\item{pdata }{the data after preprocessing.}
}

\references{
1. Savitzky, A.; Golay, M.J.E. (1964). Smoothing and Differentiation of Data by Simplified Least Squares Procedures. Analytical Chemistry, 36 (8), 1627–1639
}
\examples{
## Examples of using Savitzky-Golay filter

library(mdatools)
data(Simdata)

# take first 20 spectra and add a noise
odata = spectra.c[1:20, ]
noise = runif(length(odata)) * 0.05;
odata = odata + noise

# apply savgol with and without derivative
pdata1 = prep.savgol(odata, 9, 1, 0)
pdata2 = prep.savgol(odata, 11, 1, 1)

# show the results
par(mfrow = c(3, 1))
mdaplot(cbind(wavelength, t(odata)), type = 'l', main = 'Noisy')
mdaplot(cbind(wavelength, t(pdata1)), type = 'l', main = 'Smoothed')
mdaplot(cbind(wavelength, t(pdata2)), type = 'l', main = 'Smoothed + 1st derivative')
par(mfrow = c(1, 1))

}

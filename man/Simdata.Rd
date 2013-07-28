\name{Simdata}
\alias{Simdata}
\docType{data}
\title{
Spectral data of polyaromatic hydrocarbons mixing
}
\description{
Simdata contains training and test set with spectra and concentration values of polyaromatic hydrocarbons mixings.
}

\usage{data(Simdata)}
\format{
   The data includes four matrices:
  \describe{
    \item{\code{spectra.c}}{a matrix (100x150) with spectral values for the training set.}
    \item{\code{spectra.t}}{a matrix (100x150) with spectral values for the test set}
    \item{\code{conc.c}}{a matrix (100x3) with concentration of components for the training set.}
    \item{\code{conc.t}}{a matrix (100x3) with concentration of components for the test set}
  }
}

\details{
This is a simulated data containing UV/Vis spectra of three component mixings - C1, C2 and C3. The spectral range is betwen 211 and 350 nm. The spectra were simulated as a linear combination of pure component spectra plus 5\% of random noise. The concentration range is (in moles): C1 [0, 1], C2 [0, 0.5], C3 [0, 0.1].

There are 100 smaples in a training set and 100 samples in a test set. The data is used for showing examples with multivariate regression methods.
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(Simdata)

pdata = list(
   train = cbind(wavelength, t(spectra.c)),
   test = cbind(wavelength, t(spectra.t))
)

par(mfrow = c(2, 2))
mdaplotg(pdata, type = 'l', legend = c('train', 'test'))
hist(conc.c[, 1], xlab = 'C1', col = 'blue', main = 'C1')
hist(conc.c[, 2], xlab = 'C2', col = 'blue', main = 'C2')
hist(conc.c[, 3], xlab = 'C3', col = 'blue', main = 'C3')
par(mfrow = c(1, 1))

}
\keyword{datasets}

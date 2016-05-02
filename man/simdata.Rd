\name{simdata}
\alias{simdata}
\docType{data}
\title{
Spectral data of polyaromatic hydrocarbons mixing
}
\description{
Simdata contains training and test set with spectra and concentration values of polyaromatic hydrocarbons mixings.
}

\usage{data(simdata)}
\format{
   The data is a list with following fields:
   \tabular{rlll}{
    \code{$spectra.c} \tab a matrix (100x150) with spectral values for the training set. \cr
    \code{$spectra.t} \tab a matrix (100x150) with spectral values for the test set. \cr
    \code{$conc.c} \tab a matrix (100x3) with concentration of components for the training set. \cr
    \code{$conc.t} \tab a matrix (100x3) with concentration of components for the test set. \cr
    \code{$wavelength} \tab a vector with spectra wavelength in nm. \cr
  }
}

\details{
This is a simulated data containing UV/Vis spectra of three component (polyaromatic hydrocarbons) mixings - C1, C2 and C3. The spectral range is betwen 210 and 360 nm. The spectra were simulated as a linear combination of pure component spectra plus 5\% of random noise. The concentration range is (in moles): C1 [0, 1], C2 [0, 0.5], C3 [0, 0.1].

There are 100 mixings in a training set and 50 mixings in a test set. The data can be used for multivariate regression examples.
}

\keyword{datasets}

\name{carbs}
\alias{carbs}
\docType{data}
\title{
Raman spectra of carbonhydrates
}
\description{
The dataset consists of Raman spectra of fructose, lactose, and ribose as well as spectra of their
mixtures.
}

\usage{data(simdata)}
\format{
   The data is a list (\code{carbs}) with the following fields:
   \tabular{rlll}{
    \code{$D} \tab a matrix (21x1401) with spectral values for the mixtures. \cr
    \code{$S} \tab a matrix (1401x3) with spectral values for the pure components. \cr
    \code{$C} \tab a matrix (21x3) with concentration of the pure components. \cr
  }
}

\details{
The dataset consists of Raman spectra of fructose, lactose, and ribose as well as spectra of their
mixtures. The original spectra were downloaded from publicly available SPECARB library [1], created
by S.B. Engelsen. The specta were truncated to the range from 200 to 1600 cm-1.

The spectra of mixtures were created by linear combinations of the original spectra:

D = CS' + E

Concentrations of the components, C, follow a simplex lattice design with four levels. Some noise calculated as a random number uniformly distributed between 0\% and 3\% of maximum initial intensity (E) was added to each spectrum of the dataset, D, individually.
}

\references{
1. Engelsen S.B., Database on Raman spectra of carbohydrates. Available at: http://www.models.life.ku.dk/~specarb/specarb.html [visited 31.05.2020]
}

\keyword{datasets}

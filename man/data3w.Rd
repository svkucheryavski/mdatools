\name{data3w}
\alias{data3w}
\docType{data}
\title{Simulated 3-way data}
\description{
The dataset consists of three simulated 3-way arrays, which can be used for learning 3-way
classification methods \code{\link{ddsimca.parafac}} and \code{\link{ddsimca.tucker}}.
}

\usage{data(data3w)}
\format{
   The dataset is a list with three 3-way arrays:
   \tabular{rlll}{
    \code{$X3w} \tab 80x50x40 array representing training set with 80 target class samples . \cr
    \code{$X3wt} \tab 20x50x40 array representing test set with 20 target class samples. \cr
    \code{$X3wa} \tab 60x50x40 array representing test set with 60 alternative class samples. \cr
  }
}

\details{
The dataset was sumulated using procedure similar to what is described in [1]. Target samples contain
two components in random concentrations. Alternative class samples contain the same two components plys
a third component in random concentrations. The component instrumental profiles are partially
overlapped Gaussian functions, defined in 40 and 50 channels in each data mode
(hence the data matrices have 50x40 dimension). Gaussian noise was then added to all signals
at a level of 10\% with respect to the mean signal.
}

\references{
1. A. P. Pagani, G. Camargo, G. A. Ibañez, A. C. Olivieri, A. L. Pomerantsev, and O. Ye. Rodionova Data-Driven Version of Multiway Soft Independent Modeling of Class Analogy (N-Way DD-SIMCA): Theory and Application. Analytical Chemistry 2024 96 (12), 4845-4853. DOI: 10.1021/acs.analchem.3c05096
}

\keyword{datasets}

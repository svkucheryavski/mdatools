\name{people}
\alias{people}
\docType{data}

\title{
People data
}

\description{
Dataset for exploratory analysis with 32 objects (male and female persons) and 12 variables.
}

\usage{data(people)}
\format{
  a matrix with 32 observations (persons) and 12 variables.
   \tabular{rlll}{
      \code{[, 1]} \tab Height in cm. \cr
      \code{[, 2]} \tab Weight in kg. \cr
      \code{[, 3]} \tab Hair length (-1 for short, +1 for long). \cr
      \code{[, 4]} \tab Shoe size (EU standard). \cr
      \code{[, 5]} \tab Age, years. \cr
      \code{[, 6]} \tab Income, euro per year. \cr
      \code{[, 7]} \tab Beer consumption, liters per year. \cr
      \code{[, 8]} \tab Wine consumption, liters per year. \cr
      \code{[, 9]} \tab Sex (-1 for male, +1 for female). \cr
      \code{[, 10]} \tab Swimming ability (index, based on 500 m swimming time). \cr
      \code{[, 11]} \tab Region (-1 for Scandinavia, +1 for Mediterranean. \cr
      \code{[, 12]} \tab IQ (European standardized test). \cr
  }
}

\details{
The data was taken from the book [1] and is in fact a small subset of a pan-European demographic survey. It includes information about 32 persons, 16 represent northern Europe (Scandinavians) and 16 are from the Mediterranean regions. In both groups there are 8 male and 8 female persons. The data includes both quantitative and qualitative variables and is particularly useful for benchmarking exploratory data analysis methods.   
}

\source{
1. K. Esbensen. Multivariate Data Analysis in Practice. Camo, 2002.
}

\keyword{datasets}

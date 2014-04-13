\name{crossval}
\alias{crossval}
\title{Generate sequence of indices for cross-validation}
\usage{
crossval(nobj, nseg = NULL)
}
\arguments{
  \item{nobj}{number of objects in a dataset}

  \item{nseg}{number of segments to split the data to}
}
\value{
matrix with object indices for each segment
}
\description{
Generates and returns sequence of object indices for each
segment in random segmented cross-validation
}


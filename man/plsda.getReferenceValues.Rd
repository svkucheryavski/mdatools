\name{plsda.getReferenceValues}
\alias{plsda.getReferenceValues}
\title{Reference values for PLS-DA}
\usage{
plsda.getReferenceValues(c, classnames = NULL)
}
\arguments{
  \item{c}{vector with class values (discrete)}

  \item{classnames}{vector with names for the classes}
}
\value{
the generated matrix with one column for each class
}
\description{
Generates matrix with reference y values (-1 and +1) for a
vector with class values
}


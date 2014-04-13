\name{pls.calculateVIPScores}
\alias{pls.calculateVIPScores}
\title{VIP scores calculation for PLS model}
\usage{
pls.calculateVIPScores(object)
}
\arguments{
  \item{object}{a PLS model (object of class \code{pls})}
}
\value{
matrix \code{nvar x ny} with VIP score values for selected
number of components
}
\description{
Calculates VIP (Variable Importance in Projection) scores
for each component and response variable in the PLS model
}


\name{mdatools}
\alias{mdatools}

\title{
Package for Multivariate Data Analysis (Chemometrics)
}

\description{
This package contains classes and functions for most common methods used in Chemometrics. For a complete list of functions, use \code{library(help = 'mdatools')}.
}

\details{
The project is hosted on GitHub (https://svkucheryavski.github.io/mdatools/), there you can also find a Bookdown user tutorial explaining most important features of the package. There is also a dedicated YouTube channel (https://www.youtube.com/channel/UCox0H4utfMq4FIu2kymuyTA) with introductory Chemometric course with examples based on mdatools functionality.

Every method is represented by two classes: a model class for keeping all parameters and information about the model, and a class for keeping and visualising results of applying the model to particular data values.

Every model class, e.g. \code{\link{pls}}, has all needed functionality implemented as class methods, including model calibration, validation (test set and cross-validation), visualisation of the calibration and validation results with various plots and summary statistics.

So far the following modelling and validation methods are implemented:
\tabular{ll}{
      \code{\link{pca}}, \code{\link{pcares}} \tab Principal Component Analysis (PCA).\cr
      \code{\link{pls}}, \code{\link{plsres}} \tab Partial Least Squares regression (PLS).\cr
      \code{\link{simca}}, \code{\link{simcares}} \tab Soft Independent Modelling of Class Analogues (SIMCA)\cr
      \code{\link{simcam}}, \code{\link{simcamres}} \tab SIMCA for multiple classes case (SIMCA)\cr
      \code{\link{plsda}}, \code{\link{plsdares}} \tab Partial Least Squares Dscriminant Analysis (PLS-DA).\cr
      \code{\link{randtest}} \tab Randomization test for PLS-regression.\cr
      \code{\link{ipls}} \tab Interval PLS variable.\cr
      \code{\link{mcrals}} \tab Multivariate Curve Resolution with Alternating Least Squares.\cr
      \code{\link{mcrpure}} \tab Multivariate Curve Resolution with Purity approach.\cr
      \code{\link{pcv}} \tab Procrustes Cross Validation.\cr
}

Methods for data preprocessing:
\tabular{ll}{
      \code{\link{prep.autoscale}} \tab data mean centering and/or standardization.\cr
      \code{\link{prep.savgol}} \tab Savitzky-Golay transformation.\cr
      \code{\link{prep.snv}} \tab Standard normal variate.\cr
      \code{\link{prep.msc}} \tab Multiplicative scatter correction.\cr
      \code{\link{prep.norm}} \tab Spectra normalization.\cr
      \code{\link{prep.alsbasecorr}} \tab Baseline correction with Asymmetric Least Squares.\cr
}

All plotting methods are based on two functions, \code{\link{mdaplot}} and \code{\link{mdaplotg}}. The functions extend the basic functionality of R plots and allow to make automatic legend and color grouping of data points or lines with colorbar legend, automatically adjust axes limits when several data groups are plotted and so on.
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

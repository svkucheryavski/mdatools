\name{classres}
\alias{classres}

\title{
Results of classification
}

\description{
\code{classres} is used to store results classification for one or multiple classes.
}

\usage{
classres()
}
\arguments{
  \item{pres}{results of PCA decomposition of data (class \code{pcares}).}
  \item{cres}{results of classification (class \code{classres}).}
}

\details{
There is no need to create a \code{classres} object manually, it is created automatically when build a classification model (e.g. using \code{\link{simca}} or \code{\link{plsda}}) or apply the model to new data. For any classification method from \code{mdatools}, a class using to represent results of classification (e.g. \code{\link{simcares}}) inherits fields and methods of \code{classres}.
}

\value{
\item{c.pred}{predicted class values (+1 or -1).}
\item{c.ref}{reference (true) class values if provided.}

The following fields are available only if reference values were provided.
\item{tp}{number of true positives.}
\item{fp}{nmber of false positives.}
\item{fn}{number of false negatives.}
\item{specificity}{specificity of predictions.}
\item{sensitivity}{sensitivity of predictions.}
}

\seealso{
Methods \code{classres} class:
  \tabular{ll}{
      \code{\link{showPredictions.classres}} \tab shows table with predicted values.\cr
      \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
      \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components values.\cr
      \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components values.\cr
      \code{\link{plotPerformance.classres}} \tab makes plot with both specificity and sensitivity values.\cr
   }
}

\examples{
see examples for \code{\link{simcares}}, \code{\link{simcamres}} and \code{\link{plsdares}}.
}

\name{regresult}
\alias{regresult}
\title{
Regression Results
}
\description{
\code{regresult} is used to represent and visualise regression results: predictions and model performance.
}
\usage{
regresult(yp, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yp}{ response values, predicted by a regression model. }
  \item{y}{ reference response values}
}
\details{
The regression results are returned by a corresponding model (e.g. \code{mlr} or \code{pls}), there is no need to use the class manually. The class is used to show summary for the results, if reference values are provided, and plot with predictions.
}
\value{
   The prediction statistics (RMSE, R2, slope and bias) are calculated only if reference values were provided.
  \item{yp}{ response values, predicted by a regression model. }
  \item{y}{ reference response values}
  \item{rmse}{ root mean squared error of predictions }
  \item{slope}{ slope of fitted line for predicted vs. reference values }
  \item{r2}{ determination coefficient for predicted vs reference values}
  \item{bias}{ bias (average error of predictions) }
}
\references{
}
\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}
\note{
}

\seealso{
\code{\link{plot.regresult}} shows plots with predicted y values (or predicted vs. reference y values if the were provided), \code{\link{summary.regresult}} shows summary statistics for predictions, if reference values were provided.
}
\examples{
}

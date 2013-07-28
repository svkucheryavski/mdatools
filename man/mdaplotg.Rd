\name{mdaplotg}
\alias{mdaplotg}

\title{
Plotting function for several sets of objects
}

\description{
\code{mdaplotg} is used to make scatter, line or bar plots or their combination for several sets of objects.
}

\usage{
mdaplotg(data, type = 'p', single.x = T, colmap = 'default', 
   pch = 16, lty = 1, lwd = 1, bwd = 0.8, ylim = NULL, xlim = NULL, 
   xlab = NA, ylab = NA, main = NULL, 
   legend = NULL, legend.position = "topright", 
   labels = NULL, show.labels = F, show.lines = F, show.grid = T, ...)
}

\arguments{
  \item{data}{a list with data values (see details below).}
  \item{type}{type of the plot ('p', 'l', 'b', 'h')}
  \item{single.x}{logical, is first column of data matrix used for x values (or every odd)?}
  \item{colmap}{a colormap to use for coloring the plot items.}
  \item{pch}{a character for markers (same as \code{plot} parameter).}
  \item{lty}{the line type (same as \code{plot} parameter).}
  \item{lwd}{the line width (thickness) (same as \code{plot} parameter).}
  \item{bwd}{a width of a bar as a percent of a maximum space available for each bar.}
  \item{xlim}{limits for the x axis (if NULL, will be calculated automatically).}
  \item{ylim}{limits for the y axis (if NULL, will be calculated automatically).}
  \item{xlab}{a title for the x axis (same as \code{plot} parameter).}
  \item{ylab}{a title for the y axis (same as \code{plot} parameter).}
  \item{main}{an overall title for the plot (same as \code{plot} parameter).}
  \item{legend}{a vector with legend elements (if NULL, no legend will be shown).}
  \item{legend.position}{position of the legend ('topleft', 'topright', 'top', 'bottomleft', 'bottomright', 'bottom').}
  \item{labels}{a matrix with text labels for data points (if NULL, row names for each group will be used).}
  \item{show.labels}{logical, show or not labels for the data objects.}
  \item{show.lines}{vector with two coordinates (x, y) to show horizontal and vertical line cross the point.}
  \item{show.grid}{logical, show or not a grid for the plot.}
  \item{\dots}{}
}
\details{
The \code{mdaplotg} function is used to make a plot with several sets of objects. Simply speaking, use it when you need a plot with legend. For example to show line plot with spectra from calibration and test set, scatter plot for height and weight values for women and men, and so on.

Most of the parameters are similar to \code{\link{mdaplot}}, the difference is described below.

The data should be organized as a list, every item is a matrix with data for one set of objects. See examples for details.

There is no color grouping option, because color is used to separate the sets. Marker symbol, line style and type, etc. can be defined as a single value (one for all sets) and as a vector with one value for each set.
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## Examples of using mdaplotg

# show scatter plot for wine and beer consumption for 
# each group of people data

library(mdatools)

# show scatter plot with height and weight values for males and females

data(People)

sex = people[, 'Sex']
fdata = cbind(people[sex == 1, 'Height'], people[sex == 1, 'Weight'])
mdata = cbind(people[sex == -1, 'Height'], people[sex == -1, 'Weight'])
colnames(fdata) = colnames(mdata) = c('Height', 'Weight')
pdata = list(fdata, mdata)

par(mfrow = c(2, 2))
mdaplotg(pdata, legend = c('female', 'male'))
mdaplotg(pdata, legend = c('female', 'male'), pch = c(17, 19), colmap = c('green', 'orange'))
mdaplotg(pdata, legend = c('female', 'male'), legend.position = 'topleft', colmap = 'gray')
mdaplotg(pdata, legend = c('female', 'male'), legend.position = 'topleft', show.labels = T)
par(mfrow = c(1, 1))

# show line plots with spectra from calibration and test set of the Simdata

data(Simdata)

cdata = cbind(wavelength, t(spectra.c))
tdata = cbind(wavelength, t(spectra.t))
pdata = list(cdata, tdata)

par(mfrow = c(2, 2))
mdaplotg(pdata, type = 'l', legend = c('cal', 'test'))
mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), xlab = 'Wavelength, nm', ylab = 'Log(1/R)')
mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), lty = c(3, 2))
mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), lwd = 2, colmap = c('green', 'orange'))
par(mfrow = c(1, 1))

}


\name{mdaplot}
\alias{mdaplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting function for a single set of objects
}

\description{
\code{mdaplot} is used to make a scatter, line or a bar plot for one set of data objects.
}

\usage{
mdaplot(data, type = 'p', single.x = T, cgroup = NULL, colmap = 'default'
   pch = 16, col = NULL, lty = 1, lwd = 1, bwd = 0.8, 
   xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
   labels = NULL, 
   show.labels = F, show.colorbar = T, show.lines = NULL, 
   show.grid = T, show.axes = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix with data values.}
  \item{type}{type of the plot ('p', 'l', 'b', 'h')}
  \item{single.x}{logical, is first column of data matrix used for x values (or every odd)?}
  \item{cgroup}{a vector with values to use for make color groups.}
  \item{colmap}{a colormap to use for coloring the plot items.}
  \item{pch}{a character for markers (same as \code{plot} parameter).}
  \item{col}{a color for markers or lines (same as \code{plot} parameter).}
  \item{lty}{the line type (same as \code{plot} parameter).}
  \item{lwd}{the line width (thickness) (same as \code{plot} parameter).}
  \item{bwd}{a width of a bar as a percent of a maximum space available for each bar.}
  \item{xlim}{limits for the x axis (if NULL, will be calculated automatically).}
  \item{ylim}{limits for the y axis (if NULL, will be calculated automatically).}
  \item{xlab}{a title for the x axis (same as \code{plot} parameter).}
  \item{ylab}{a title for the y axis (same as \code{plot} parameter).}
  \item{main}{an overall title for the plot (same as \code{plot} parameter).}
  \item{labels}{a vector with text labels for data points (if NULL, row names will be used).}
  \item{show.labels}{logical, show or not labels for the data objects.}
  \item{show.colorbar}{logical, show or not colorbar legend if color grouping is on.}
  \item{show.lines}{vector with two coordinates (x, y) to show horizontal and vertical line cross the point.}
  \item{show.grid}{logical, show or not a grid for the plot.}
  \item{show.axes}{logical, make a normal plot or show only elements (markers, lines, bars) without axes.}
  \item{\dots}{}
}

\details{
Most of the parameters are similar to what are used with \code{plot} function. The differences are described below.

The function makes a plot of one set of values. It can be a set of points (scatter plot), bars, lines or scatter-lines. The data is organized as a matrix, in fact for scatter and bar plot only two columns are needed, but to plot of set of lines there are two ways to organise the data matrix.

If \code{single.x = T} it means first column of the matrix contains x values, as the rest columns contain y values (one column - one line, e.g. [x y1 y2 y3]). If \code{single.x = F} it means every odd column of matrix contains x values and every even column - y values of a line (e.g. [x1 y1 x2 y2 x3 y3]). This option is needed if lines have different x values.

The function allows to colorize lines and points according to values of a parameter \code{cgroup}. The parameter must be a vector with the same elements as number of objects (rows) in the data. The values are divided into up to eight intervals and for each interval a particular color from a selected color scheme is assigned. Parameter \code{show.colorbar} allows to turn off and on a color bar legend for this option.

The used color scheme is defined by the \code{colmap} parameter. The default scheme is based on color brewer (colorbrewer2.org) diverging scheme with eight colors. There is also a gray scheme (\code{colmap = 'gray'}) and user can define its own just by specifing the needed sequence of colors (e.g. \code{colmap = c('red', 'yellow', 'green')}, two colors is minimum). The scheme will then be generated automatically as a gradient among the colors.  

}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
\code{\link{mdaplotg}} - to make plot with several sets of data (groups).
}

\examples{
## Examples of using mdaplot function

library(mdatools)
data(Simdata)

# make a line plot for y = x^2 with different x values for each line

x1 = seq(-10, 10, length = 100) 
y1 = x1^2
x2 = seq(-9, 9, length = 100) + 1
y2 = x2^2 * 1.2
x3 = seq(-8, 8, length = 100) + 2
y3 = x3^2 * 1.4

mdaplot(cbind(x1, y1, x2, y2, x3, y3), type = 'l', single.x = F)


# make line plot of the spectra with coloring by concentration of first component
# using different color maps

par(mfrow = c(2, 2))
mdaplot(cbind(wavelength, t(spectra.c)), type = 'l', cgroup = conc.c[, 1])
mdaplot(cbind(wavelength, t(spectra.c)), type = 'l', cgroup = conc.c[, 1], 
   colmap = 'gray', show.colorbar = F)
mdaplot(cbind(wavelength, t(spectra.c)), type = 'l', cgroup = conc.c[, 1], 
   colmap = c('red', 'green'))
mdaplot(cbind(wavelength, t(spectra.c)), type = 'l', cgroup = conc.c[, 1], 
   colmap = c('#ffff00', '#ff00ff', '#00ffff'))
par(mfrow = c(1, 1))

# show scatter plots from spectral data with coloring and other parameters
# see how limits are adjusted if show.lines option is used

nobj = 30
pdata = cbind(spectra.c[1:nobj, 100, drop = F], spectra.c[1:nobj, 110, drop = F])
coldata = conc.c[1:nobj, 1]

par(mfrow = c(2, 2))
mdaplot(pdata, cgroup = coldata, main = 'Spectra', show.lines = c(0.1, 0.06))
mdaplot(pdata, cgroup = coldata, main = 'Spectra', show.lines = c(0.06, 0.01))
mdaplot(pdata, cgroup = coldata, xlab = '309 nm', ylab = '319 nm', 
   main = 'Spectra', show.labels = T, show.grid = F)
mdaplot(pdata, col = 'red', pch = 17, xlab = '309 nm', ylab = '319 nm',
   main = 'Spectra', show.labels = T, labels = paste('Obj', 1:nobj))
par(mfrow = c(1, 1))

}

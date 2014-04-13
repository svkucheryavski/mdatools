\name{mdaplot.showLabels}
\alias{mdaplot.showLabels}
\title{Plot labels
Shows labels for data elements (points, bars) on a plot.}
\usage{
mdaplot.showLabels(data, pos = 3, cex = 0.65, col = "darkgray",
  type = NULL)
}
\arguments{
  \item{data}{data matrix with coordinates of the points
  (x, y)}

  \item{pos}{position of the labels relative to the points}

  \item{cex}{size of the labels text}

  \item{col}{color of the labels text}

  \item{type}{type of the plot}
}
\description{
Plot labels Shows labels for data elements (points, bars)
on a plot.
}
\details{
Rownames of matrix \code{data} are used as labels. If
matrix has no rownames, row numbers will be used instead.
}


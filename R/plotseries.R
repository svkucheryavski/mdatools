#' Create plot series object based on data, plot type and parameters
#'
#' @param data
#' data to make the plot for (vector, matrix or data frame).
#' @param type
#' type of the plot.
#' @param cgroup
#' vector with values used to create a color grouping of the series instances.
#' @param col
#' color to show the series on plot with (user defined).
#' @param opacity
#' opacity of the colors (between 0 and 1).
#' @param colmap
#' colormap name to generate color/colors if they are not specified by user. See
#' \code{link{mdaplot.getColors}} for details.
#' @param labels
#' either vector with labels for the series instances or string ("names", "values", or "indices")
#' if labels should be generated automatically.
#'
#' @description
#' The `plotseries` object contains all necessary paremeters to create main plots
#' from data values, including values for x and y, correct handling of excluded rows
#' and columns, color grouping (if any), limits and labels.
#'
#' If both `col` and `cgroup` are specified, `cgroup` will be ignored.
#'
#' Labels can be either provided by user or generated automatically based on values, names or
#' indices of data rows and columns. If series is made for scatter plot `type="p"` then labels
#' are required for each row of the original dataset. Otherwise (for line, bar and errobar plot)
#' labels correspond to data columns (variables).
#'
#' The object has the following plotting methods once created:
#' \code{\link{plotScatter}}
#' \code{\link{plotLines}}
#' \code{\link{plotBars}}
#' \code{\link{plotDensity}}
#' \code{\link{plotErrorbars}}
#'
#' @export
plotseries <- function(data, type, cgroup = NULL, col = NULL, opacity = 1,
   colmap = "default", labels = NULL) {

   if (length(data) == 0) {
      stop("Seems you forgot to provide data values (vector, matrix of data frame).")
   }

   ps <- list()
   ps$type <- type
   ps$name <- attr(data, "name", exact = TRUE)
   ps$call <- match.call()
   class(ps) <- "plotseries"

   # prepare data values
   plot_data <- preparePlotData(data)
   ps <- c(ps, splitPlotData(plot_data$visible, type))

   # add values for excluded rows
   if (!is.null(plot_data$excluded)) {
      ps$excluded_cols <- plot_data$excluded_cols
      ps$excluded_rows <- plot_data$excluded_rows
      ps <- c(ps, splitExcludedData(plot_data$excluded, type))
   }

   # add labels
   ps <- c(ps, getDataLabels(ps, labels))

   # handle colors
   ps <- c(ps, getPlotColors(ps, col, opacity, cgroup, colmap))

   # compute axes limits
   ps$xlim <- range(ps$x_values)
   ps$ylim <- range(ps$y_values)
   class(ps) <- "plotseries"

   return(ps)
}

#' Take dataset and prepare them for plot
#'
#' @param data
#' dataset (vector, matrix or data frame)
#'
#' @description
#' The function checks that `data` contains correct numeric values, check for
#' mandatory attributes (row and column names, x- and y-axis values and names, etc.)
#' and add them if necessary.
#'
#' Another things is to remove hidden columns and split the rest to visible and hidden
#' values (if excluded rows are present).
#'
preparePlotData <- function(data) {

   # get data attributes to separate variable for shorten code
   attrs <- attributes(data)

   # if it is vector without dimension - make a matrix
   if (is.null(dim(data))) {
      data <- t(as.matrix(data))
   }

   # convert data frame to a matrix if needed
   if (is.data.frame(data)) {
      data <- mda.df2mat(data)
   }

   # if row names are missing - add them
   if (is.null(rownames(data))) {
      rownames(data) <- paste0("O", seq_len(nrow(data)))
   }

   # if column names are missing - add them
   if (is.null(colnames(data))) {
      colnames(data) <- paste0("X", seq_len(ncol(data)))
   }

   # if xaxis values are missing - add them
   if (is.null(attrs$xaxis.values)) {
      attrs$xaxis.values <- seq_len(ncol(data))
   }

   # if yaxis values are missing - add them
   if (is.null(attrs$yaxis.values)) {
      attrs$yaxis.values <- seq_len(nrow(data))
   }

   # handle excluded columns if any
   excluded_cols <- NULL
   if (length(attrs$exclcols) > 0) {
      excluded_cols <- mda.getexclind(attrs$exclcols, colnames(data), ncol(data))
      data <- data[, -excluded_cols, drop = FALSE]
      attrs$xaxis.values <- attrs$xaxis.values[-excluded_cols]
   }

   # check if data still has some columns
   if (is.null(ncol(data)) || ncol(data) == 0) {
      stop("No columns left after excluded hidden values.")
   }

   # handle excluded rows if any
   excluded_data <- NULL
   excluded_rows <- NULL
   if (length(attrs$exclrows > 0)) {
      excluded_rows <- mda.getexclind(attrs$exclrows, rownames(data), nrow(data))
      excluded_data <- data[excluded_rows, , drop = FALSE]
      excluded_yaxis_values <- attrs$yaxis.values[excluded_rows]

      data <- data[-excluded_rows, , drop = FALSE]
      attrs$yaxis.values <- attrs$yaxis.values[-excluded_rows]
   }

   # reassign some attribues to data
   attr(data, "xaxis.values") <- attrs$xaxis.values
   attr(data, "yaxis.values") <- attrs$yaxis.values
   attr(data, "xaxis.name") <- attrs$xaxis.name
   attr(data, "yaxis.name") <- attrs$yaxis.name
   attr(data, "name") <- attrs$name

   # add necessary attributes to excluded data
   if (!is.null(excluded_data)) {
      attr(excluded_data, "xaxis.values") <- attrs$xaxis.values
      attr(excluded_data, "yaxis.values") <- excluded_yaxis_values
   }

   return(list(
      visible = data,
      excluded = excluded_data,
      excluded_cols = excluded_cols,
      excluded_rows = excluded_rows
   ))
}

#' Split dataset to x and y values depending on plot type
#'
#' @param data
#' matrix with data values (visible or hidden)
#' @param type
#' type of plot
#'
splitPlotData <- function(data, type) {

   # shortcuts to some of parameters
   attrs <- attributes(data)


   if (type == "p" && ncol(data) == 1) {
      # if data has only one column add y values in front
      data <- cbind(attrs$yaxis.values, data)
      colnames(data)[1] <- if (is.null(attrs$yaxis.name)) "Objects" else attrs$yaxis.name
   }

   if (type %in%  c("p", "d")) {
      # scatter or density plot - take first and second columns as x and y
      x_values <- data[, 1]
      names(x_values) <- rownames(data)
      attr(x_values, "name") <- colnames(data)[1]
      y_values <- data[, 2, drop = FALSE]
      attr(y_values, "name") <- colnames(data)[2]
      return(list(x_values = x_values, y_values = y_values))
   }

   if (type == "e") {
      # errorbar plor - make data consist of three rows: (m, m - err, m + err)
      data <- rbind(
         data[1, ],
         data[1, ] - data[2, ],
         data[1, ] + if (nrow(data) > 2) data[3, ] else data[2, ]
      )
   }


   # 0.12.0: yaxis.name must not be used as axis label in line and bar plots
   if (type %in% c("b", "l", "e", "h")) {
      attrs$yaxis.name = NULL
   }

   # prepare x-axis values for other types of plots
   y_values <- data
   x_values <- attrs$xaxis.values

   if (is.null(names(x_values))) names(x_values) <- colnames(data)
   attr(x_values, "name") <- if (is.null(attrs$xaxis.name)) "Variables" else attrs$xaxis.name

   # by default take all data values as y and assign name to the x-variable
   attr(y_values, "name") <- if (is.null(attrs$yaxis.name)) "" else attrs$yaxis.name

   return(list(x_values = x_values, y_values = y_values))
}

#' Split the excluded part of data
#'
#' @param data
#' matrix with hidden data values
#' @param type
#' type of plot
#'
splitExcludedData <- function(data, type) {

   if (type == "e") {
      stop("Errorbar can not be created for data with excluded rows.")
   }

   pd <- splitPlotData(data, type)
   return(list(x_values_excluded = pd$x_values, y_values_excluded = pd$y_values))
}

#' Define colors for plot series
#'
#' @param ps
#' `plotseries` object
#' @param col
#' color specified by user (if any)
#' @param opacity
#' opacity for the color
#' @param cgroup
#' vector for color grouping (if any)
#' @param colmap
#' name or values for colormap
#'
getPlotColors <- function(ps, col, opacity, cgroup, colmap) {

   # if user specified color - use it
   if (!is.null(col)) {
      return(list(
         col = adjustcolor(col, opacity),
         colmap = colmap,
         cgroup = NULL
      ))
   }

   # if user did not specify cgroup - get one color based on colmap and opacity
   if (is.null(cgroup) || ps$type == "e") {
      return(list(
         col = mdaplot.getColors(1, colmap = colmap, opacity = opacity),
         colmap = colmap,
         cgroup = NULL
      ))
   }

   # check if cgroup is provided for all values, also from excluded rows and correct this
   if (ps$type == "h") {
      cgroup_expected_length <- length(ps$x_values)
      cgroup_excluded_values <- ps$excluded_cols
   } else {
      cgroup_expected_length <- nrow(ps$y_values)
      cgroup_excluded_values <- ps$excluded_rows
   }

   if (length(cgroup) != cgroup_expected_length && length(cgroup_excluded_values) > 0) {
      cgroup <- cgroup[-cgroup_excluded_values]
   }

   if (length(cgroup) != cgroup_expected_length) {
      stop("Parameter 'cgroup' does not match size of the dataset.")
   }

   return(list(
      col = mdaplot.getColors(cgroup = cgroup, colmap = colmap, opacity = opacity),
      colmap = colmap,
      cgroup = cgroup
   ))
}

#' Compute confidence ellipse for a set of points
#'
#' @param points
#' matrix of data frame with coordinates of the points
#' @param conf.level
#' confidence level for the ellipse
#' @param n
#' number of points in the ellipse coordinates
#'
#' @return
#' matrix with coordinates of the ellipse points (x and y)
#'
getConfidenceEllipse <- function(points, conf.level = 0.95, n = 100) {

   # compute igen vectors and values
   e <- eigen(cov(points));

   # get angle between the x-axis and the largest eigenvector
   phi <- atan2(e$vectors[[2]], e$vectors[[1]])
   if (phi < 0) phi <- phi + 2 * pi

   # compute center and radii
   chisq <- qchisq(conf.level, 2)
   center <- apply(points, 2, mean)
   radii <- sqrt(chisq * e$values)

   # generate vector of angles
   theta_grid <- seq(0, 2 * pi, length.out = n)

   # compute x and y coordinates of the ellipse
   ellipse <- cbind(radii[1] * cos(theta_grid), radii[2] * sin(theta_grid))

   # rotate the coordinates to the angle found earlier and shift to center of points
   R <- rbind(c(cos(phi), sin(phi)), c(-sin(phi), cos(phi)))
   ellipse <- sweep(ellipse %*% R, 2, center, "+")

   return(ellipse)
}

#' Compute coordinates of a closed convex hull for data points
#'
#' @param points
#' matrix of data frame with coordinates of the points
#'
#' @importFrom grDevices chull
#'
getConvexHull <- function(points) {
   ch_ind <- chull(points)
   ch_ind <- c(ch_ind, ch_ind[1])
   return(points[ch_ind, ])
}

#' Create labels from data values
#'
#' @param ps
#' `plotseries` object
#'
getLabelsAsValues <- function(ps) {
   y_values <- ps$y_values
   y_values_excluded <- ps$y_values_excluded

   if (ps$type %in% c("l", "b")) {
      labels <- apply(y_values, 2, max)
      labels_excluded <- if (!is.null(y_values_excluded)) apply(y_values_excluded, 2, max)
      return(list(labels = labels, labels_excluded = labels_excluded))
   }

   if (ps$type %in% c("h", "e")) {
      labels <- y_values[1, ]
      labels_excluded <- y_values_excluded[1, ]
      return(list(labels = labels, labels_excluded = labels_excluded))
   }

   # otherwise prepare labels as for scatter plot
   labels <- y_values
   labels_excluded <- y_values_excluded
   return(list(labels = labels, labels_excluded = labels_excluded))
}

#' Create labels as column or row indices
#'
#' @param ps
#' `plotseries` object
#'
getLabelsAsIndices <- function(ps) {

   # function which returns indices for columns or rows
   f <- function(nvalues, excluded_values) {
      indices <- 1:(nvalues + length(excluded_values))
      if (!is.null(excluded_values)) {
         indices <- indices[-excluded_values]
      }
      return(list(labels = indices, labels_excluded = excluded_values))
   }

   if (ps$type == "p") {
      # if scatter plot row indices are used and excluded rows also have indices
      return(f(nrow(ps$y_values), ps$excluded_rows))
   }

   # if non scatter plot, column indices are used and no hidden cols shown
   out <- f(length(ps$x_values), ps$excluded_cols)
   return(list(labels = out$labels, labels_excluded = NULL))
}

#' Create a vector with labels for plot series
#'
#' @param ps
#' `plotseries` object
#' @param labels
#' vector with user defined labels or type of labels to show ("values", "names", "indices")
#'
#' @description
#' For scatter plots labels correspond to rows of the data (names, values, indices, etc.). For
#' non-scatter plots labels correspond to the columns (names, indices or max value for each column)
#'
getDataLabels <- function(ps, labels = NULL) {
   y_values <- ps$y_values
   x_values <- ps$x_values

   # if labels are not specified - use names
   if (is.null(labels)) labels <- "names"

   # if user provided labels - use them
   if (length(labels) > 1) {
      excluded_cols <- ps$excluded_cols
      excluded_rows <- ps$excluded_rows

      # get values and excluded values (rows or cols depending on plot type)
      n_values <- if (ps$type == "p") nrow(y_values) else length(x_values)
      excluded_values <- if (ps$type == "p") excluded_rows else excluded_cols

      # check that labels were specified for all values
      if (length(labels) != (n_values + length(excluded_values))) {
         stop("Labels must be provided for all values (also ones which are excluded).")
      }

      if (is.null(excluded_values)) {
         labels_excluded <- NULL
         labels <- labels
         return(list(labels = labels, labels_excluded = labels_excluded))
      }

      labels_excluded <- labels[excluded_values]
      labels <- labels[-excluded_values]
      return(list(labels = labels, labels_excluded = labels_excluded))
   }

   # if labels must be values use y-values for that
   if (labels == "values") {
      return(getLabelsAsValues(ps))
   }

   # if labels must be indices use row or column indices
   if (labels == "indices") {
      return(getLabelsAsIndices(ps))
   }

   # if nothing above works - use names as labels
   labels <- names(x_values)
   labels_excluded <- if (ps$type == "p") names(ps$x_values_excluded) else NULL
   return(list(labels = labels, labels_excluded = labels_excluded))
}

#' Show labels on plot
#'
#' @param ps
#' `plotseries` object
#' @param show.excluded
#' logical, are excluded rows also shown on the plot
#' @param pos
#' position of the labels relative to the data points
#' @param cex
#' size of the labels text
#' @param col
#' color of the labels text
#' @param force.x.values
#' vector with forced x-values (or NULL)
#' @param bwd
#' bar width in case of bar plot
#'
showLabels <- function(ps, show.excluded = FALSE, pos = 3, cex = 0.65, col = "darkgray",
   force.x.values = NULL, bwd = 0.8) {

   f <- function(x, y, labels, type) {

      if (length(labels) == 0) stop("No labels available.")
      if (ps$type %in% c("h", "e")) y <- y[1, ]
      if (ps$type %in% c("l", "b")) y <- apply(y, 2, max)

      # correct x_values if they were forced by bwd
      if (is.numeric(force.x.values)) {
         x <- x - bwd / 2 + (force.x.values[1] - 0.5) * bwd / force.x.values[2]
         bwd <- bwd / force.x.values[2]
      }

      x <- as.numeric(x)
      y <- as.numeric(y)
      labels <- mdaplot.formatValues(labels)

      if (ps$type == "h") pos <- ifelse(y < 0, 1, 3)
      text(x, y, labels, cex = cex, pos = pos, col = col)
   }

   f(ps$x_values, ps$y_values, ps$labels, ps$type)

   if (show.excluded && !is.null(ps$labels_excluded)) {
      f(ps$x_values_excluded, ps$y_values_excluded, ps$labels_excluded, ps$type)
   }
}

#' Show plot series as set of points
#'
#' @param ps
#' `plotseries` object
#' @param pch
#' size of point markers
#' @param col
#' color of the points
#' @param bg
#' background color of the points if `pch=21:25`
#' @param lwd
#' line width for the error bars
#' @param cex
#' scale factor for the marker
#' @param col.excluded
#' color for excluded values (if must be shown)
#' @param pch.colinv
#' logical, should `col` and `bg` be switched if `pch=21:25` and `cgroup` is used to create colors.
#' @param show.excluded
#' logical, show or not the excluded data points
#' @param ...
#' other arguments for function `points()`.
#'
#' @export
plotScatter <- function(ps, pch = 16, col = ps$col, bg = "white", lwd = 1, cex = 1,
   col.excluded = "lightgray", pch.colinv = FALSE, show.excluded = FALSE, ...) {

   if (pch.colinv) {
      # switch colors for main and background in case pch=21:25
      bg.old <- bg
      bg <- col
      col <- bg.old
   }

   # show main set of points
   points(ps$x_values, ps$y_values, col = col, bg = bg, pch = pch, lwd = lwd, cex = cex, ...)

   # show excluded points if any
   if (show.excluded && !is.null(ps$y_values_excluded)) {
      points(ps$x_values_excluded, ps$y_values_excluded, col = col.excluded, pch = pch,
         lwd = lwd, cex = cex, ...)
   }
}

#' Show plot series as set of lines
#'
#' @param ps
#' `plotseries` object
#' @param col
#' a color for markers or lines (same as \code{plot} parameter).
#' @param lty
#' line type
#' @param lwd
#' line width
#' @param cex
#' scale factor for the marker
#' @param col.excluded
#' color for the excluded lines.
#' @param show.excluded
#' logical, show or not the excluded data points
#' @param ...
#' other arguments for function `lines()`.
#'
#' @export
plotLines <- function(ps, col = ps$col, lty = 1, lwd = 1, cex = 1,
   col.excluded = "darkgray", show.excluded = FALSE, ...) {

   # show main set of lines
   matlines(ps$x_values, t(ps$y_values), type = ps$type, col = col, lty = lty,
      lwd = lwd, cex = cex, ...)

   # show excluded rows
   if (show.excluded && !is.null(ps$y_values_excluded)) {
      matlines(ps$x_values_excluded, t(ps$y_values_excluded), type = ps$type,
         lty = lty, lwd = lwd, cex = cex, col = col.excluded, ...)
   }
}

#' Show plot series as error bars
#'
#' @param ps
#' `plotseries` object
#' @param col
#' color for the error bars
#' @param pch
#' marker symbol for the plot
#' @param lwd
#' line width for the error bars
#' @param cex
#' scale factor for the marker
#' @param ...
#' other arguments for function `points()`.
#'
#' @description
#' It is assumed that first row of dataset contains the y-coordinates of points,
#' second rows contains size of lower error bar and third - size for upper error bar. If
#' only two rows are provided it is assumed that error bars are symmetric.
#'
#' @export
plotErrorbars <- function(ps, col = ps$col, pch = 16, lwd = 1, cex = 1, ...) {

   x <- ps$x_values
   y <- ps$y_values

   dx <- diff(ps$xlim) / max(50, (length(x) * 3))

   segments(x, y[2, ], x, y[3, ], col = ps$col, lwd = lwd)
   segments(x - dx, y[2, ], x + dx, y[2, ], col = ps$col, lwd = lwd)
   segments(x - dx, y[3, ], x + dx, y[3, ], col = ps$col, lwd = lwd)
   points(x, y[1, ], col = col, pch = pch, cex = cex, ...)
}

#' Show plot series as bars
#'
#'
#' @param ps
#' `plotseries` object
#' @param col
#' colors of the bars
#' @param bwd
#' width of the bars (as a ratio for max width)
#' @param border
#' color of bar edges
#' @param force.x.values
#' vector with corrected x-values for a bar plot (needed for group plots, do not change manually).
#'
#' @description
#' First row of the data matrix is taken for creating the bar series. In case of barplot
#' color grouping is made based on columns (not rows as for all other plots).
#'
#' @export
plotBars <- function(ps, col = ps$col, bwd = 0.8, border = NA, force.x.values = NA) {
   x <- ps$x_values
   y <- ps$y_values[1, ]

   if (length(x) > 1) {
      # this gives variable width for bars and does not work well
      #bwd_left <- c(x[seq(2, length(x))] - x[seq(1, length(x) - 1)])
      #bwd_right <- -c(x[seq(1, length(x) - 1)] - x[seq(2, length(x))])
      #bwd_left <- c(bwd_left[1], bwd_left) * bwd / 2
      #bwd_right <- c(bwd_right, bwd_right[length(bwd_right)]) * bwd / 2
      dx <- min(diff(x))
      bwd_right <- bwd_left <- dx * bwd / 2
   } else {
      bwd_left <- bwd_right <- bwd * x / 2
   }

   # correct x_values if they were forced by bwd
   if (is.numeric(force.x.values)) {
      x <- x - bwd / 2 + (force.x.values[1] - 0.5) * bwd / force.x.values[2]
      bwd_left <- bwd_right <- bwd / force.x.values[2] / 2
   }

   if (length(col) != length(y)) {
      col <- rep(col, length.out = length(y))
   }

   rect(x - bwd_left, 0, x + bwd_right, y, col = col, border = border)
}

#' Show plot series as density plot (using hex binning)
#'
#' @param ps
#' `plotseries` object
#' @param nbins
#' number of bins in one dimension
#' @param colmap
#' colormap name or values used to create color gradient
#'
plotDensity <- function(ps, nbins = 60, colmap = ps$colmap) {
   x <- ps$x_values
   y <- ps$y_values

   # Compute coordinates of nearest lattice center
   getNearest <- function(value, scale, round = 1) {
      div <- floor(value / (scale / 2))
      return(scale / 2 * (div + as.numeric(div %% 2 == round)))
   }

   # Compute x or y coordinates of hexagons around given centers coordinates
   getHexpoints <- function(center_coord, r, f = sin) {
      a <- seq(30, 390, by = 60) / 180 * pi
      nrow <- length(a)
      ncol <- length(center_coord)
      hex_coord <-
         matrix(r * f(a), nrow = nrow, ncol = ncol) +
         matrix(center_coord, nrow = nrow, ncol = ncol, byrow = TRUE)
      hex_coord <- rbind(hex_coord, rep(NA, ncol))
      return(as.numeric(hex_coord))
   }

   # Compute x coordinates of hexagons
   getHexpointsX <- function(x, rx) {
      return(getHexpoints(x, rx, cos))
   }

   # Compute y coordinates of hexagons
   getHexpointsY <- function(y, ry) {
      return(getHexpoints(y, ry, sin))
   }

   # compute size of bins based on range of data values and number of bins
   x_range <- diff(range(x))
   y_range <- diff(range(y))

   bin_width <- x_range / nbins
   bin_height <- y_range / nbins

   # compute coordinates of hexbin lattice based on points coordinates

   ## basic lattice
   xx1 <- getNearest(x, bin_width, 1)
   yy1 <- getNearest(y, bin_height * sqrt(3), 1)

   ## lattice shifted by 0.5 of bin width/height related to the basic
   xx2 <- getNearest(x, bin_width, 0)
   yy2 <- getNearest(y, bin_height * sqrt(3), 0)

   # find which of the two fits every point best
   d1 <- (x - xx1)^2 + (y - yy2)^2
   d2 <- (x - xx2)^2 + (y - yy2)^2

   xx <- ifelse(d1 < d2, xx1, xx2)
   yy <- ifelse(d1 < d2, yy1, yy2)

   # compute number of points around each lattice center
   density <- as.data.frame(table(xx, yy, exclude = FALSE))
   density[, 1] <- as.numeric(as.character(density[, 1]))
   density[, 2] <- as.numeric(as.character(density[, 2]))

   # remove centers with no points around
   density <- density[density[, 3] > 0, ]
   lattice_x <- density[, 1]
   lattice_y <- density[, 2]
   density <- density[, 3]

   # compute coordinates of hexagons as vectors
   rx <- bin_width / 2 / cos(30 / 180 * pi)
   ry <- bin_height / 2 / cos(30 / 180 * pi)
   hex_x <- getHexpointsX(lattice_x, rx)
   hex_y <- getHexpointsY(lattice_y, ry)

   # add polygons to the current axes
   polygon(hex_x, hex_y, border = NA, col = mdaplot.getColors(cgroup = density, colmap = colmap))
}

#' Add confidence ellipse for groups of points on scatter plot
#'
#' @param p
#' plot data returned by function `mdaplot()`.
#' @param lwd
#' thickness of line used to show the hull.
#' @param lty
#' type of line used to show the hull.
#' @param conf.level
#' confidence level to make the ellipse for (between 0 and 1).
#' @param opacity
#' of opacity is 0 ellipse is transparent otherwise semi-transparent.
#'
#' @description
#' The method shows confidence ellipse for groups of points on a scatter plot made using
#' `mdaplot()` function with `cgroup` parameter. It will work only if `cgroup` is a factor.
#'
#' @examples
#' # adds 90% confidence ellipse with semi-transparent area over two clusters of points
#'
#' library(mdatools)
#' data(people)
#' group <- factor(people[, "Sex"], labels = c("Male", "Female"))
#'
#' # first make plot and then add confidence ellipse
#' p <- mdaplot(people, type = "p", cgroup = group)
#' plotConfidenceEllipse(p, conf.level = 0.90, opacity = 0.2)
#'
#' @export
plotConfidenceEllipse <- function(p, conf.level = 0.95, lwd = 1, lty = 1, opacity = 0) {

   plotPointsShape(
      p,
      lwd = lwd,
      lty = lty,
      opacity = opacity,
      shape_function = getConfidenceEllipse,
      conf.level = conf.level
   )
}

#' Add convex hull for groups of points on scatter plot
#'
#' @param p
#' plot data returned by function `mdaplot()`.
#' @param lwd
#' thickness of line used to show the hull.
#' @param lty
#' type of line used to show the hull.
#' @param opacity
#' of opacity is larger than 0 a semi-transparent polygon is shown over points.
#'
#' @description
#' The method shows convex hull for groups of points on a scatter plot made using
#' `mdaplot()` function with `cgroup` parameter. It will work only if `cgroup` is a factor.
#'
#' @examples
#' # adds convex hull with semi-transparent area over two clusters of points
#'
#' library(mdatools)
#' data(people)
#' group <- factor(people[, "Sex"], labels = c("Male", "Female"))
#'
#' p <- mdaplot(people, type = "p", cgroup = group)
#' plotConvexHull(p)
#'
#' @importFrom grDevices chull
#' @importFrom graphics polygon
#'
#' @export
plotConvexHull <- function(p, lwd = 1, lty = 1, opacity = 0) {
   plotPointsShape(
      p,
      lwd = lwd,
      lty = lty,
      opacity = opacity,
      shape_function = getConvexHull
   )
}

#' Add confidence ellipse or convex hull for group of points
#'
#' @param p
#' plot data returned by function `mdaplot()`
#' @param lwd
#' thickness of line used to show the hull
#' @param lty
#' type of line used to show the hull
#' @param opacity
#' of opacity is larger than 0 a semi-transparent polygon is shown over points
#' @param shape_function
#' function which calculates and return coordinates of the shape
#' @param ...
#' extra parameters for shape_function
#'
#' @importFrom graphics polygon
#'
#' @export
plotPointsShape <- function(p, lwd, lty, opacity, shape_function, ...) {

   x <- p$x_values
   y <- p$y_values
   cgroup <- p$cgroup
   type <- p$type

   if (type != "p") {
      stop("Shape can be added only to scatter plots.")
   }

   if (!is.factor(cgroup)) {
      stop("Parameter 'cgroup' must be a factor if you want to show points shape.")
   }

   if (length(unique(cgroup)) != length(unique(p$col))) {
      stop("Number of colors should be the same as number of levels in parameter 'cgroup'.")
   }

   col <- split(p$col, f = p$cgroup, drop = TRUE)
   d <- split(data.frame(x, y), f = p$cgroup, drop = TRUE)

   plot_function <- function(i) {
      # compute indices for convex hull points and draw the hull
      shape <- shape_function(d[[i]], ...)
      color <- col[[i]][1]
      lines(shape[, 1], shape[, 2], col = color, lwd = lwd, lty = lty)

      # if opacity is not zero add a polygon on top
      if (opacity > 0) {
         polygon(shape[, 1], shape[, 2], border = NA, col = adjustcolor(color, alpha.f = opacity))
      }
   }

   sapply(seq_len(length(d)), plot_function)
   invisible(NULL)
}

#' Add regression line for data points
#'
#' @description
#' Shows linear fit line for data points.
#'
#' @param p
#' plot data returned by function `mdaplot()`
#' @param col
#' color of line
#' @param ...
#' other parameters available for `abline()` function
#'
plotRegressionLine <- function(p, col = p$col, ...) {
   abline(lm(p$y_values ~ p$x_values), col = col, ...)
}

#' Hotelling ellipse
#'
#' @description
#' Add Hotelling ellipse to a scatter plot
#'
#' @param p
#' plot series (e.g. from PCA scores plot)
#' @param conf.lim
#' confidence limit
#' @param col
#' color of the ellipse line
#' @param lty
#' line type (e.g. 1 for solid, 2 for dashed, etc.)
#' @param ...
#' any argument suitable for \code{lines} function
#'
#' @details
#' The method is created to be used with PCA and PLS scores plots, so it shows the statistical
#' limits computed using Hotelling T^2 distribution in form of ellipse. The function works similar
#' to \code{\link{plotConvexHull}} and \code{\link{plotConfidenceEllipse}} but does not require
#' grouping of data points. Can be used together with functions \code{\link{plotScores.pca}},
#' \code{\link{plotScores.ldecomp}}, \code{\link{plotXScores.pls}},
#' \code{\link{plotXScores.plsres}}.
#'
#' See examples for more details.
#'
#' @examples
#'
#' # create PCA model for People data
#' data(people)
#' m <- pca(people, 4, scale = TRUE)
#'
#' # make scores plot and show Hotelling ellipse with default settings
#' p <- plotScores(m, xlim = c(-8, 8), ylim = c(-8, 8))
#' plotHotellingEllipse(p)
#'
#' # make scores plot and show Hotelling ellipse with manual settings
#' p <- plotScores(m, xlim = c(-8, 8), ylim = c(-8, 8))
#' plotHotellingEllipse(p, conf.lim = 0.99, col = "red")
#'
#' # in case if you have both calibration and test set, 'plotScores()' returns
#' # plot series data for both, so you have to subset it and take the first series
#' # (calibration set) as shown below.
#' ind <- seq(1, 32, by = 4)
#' xc <- people[-ind, , drop = FALSE]
#' xt <- people[ind, , drop = FALSE]
#' m <- pca(xc, 4, scale = TRUE, x.test = xt)
#'
#' p <- plotScores(m, xlim = c(-8, 8), ylim = c(-8, 8))
#' plotHotellingEllipse(p[[1]])
#'
#' @export
plotHotellingEllipse <- function(p, conf.lim = 0.95, col = "#a0a0a0", lty = 3, ...) {

   if (!inherits(p, "plotseries")) {
      stop("Argument 'p' does not look like a plot series, try 'p[[1]]' instead.")
   }

   if (p$type != "p") {
      stop("Hotelling ellipse can be added to scatter plot only.")
   }

   t1 <- as.numeric(p$x_values)
   t2 <- as.numeric(p$y_values)

   s1 <- sd(t1)^2
   s2 <- sd(t2)^2

   nobj <- length(t1)
   T2lim <- (2 * (nobj - 1) / (nobj - 2)) * qf(conf.lim, 2, (nobj - 2))

   a <- sqrt(T2lim * s1)
   b <- sqrt(T2lim * s2)

   ellipse(a = a, b = b, col = col, lty = lty, ...)
}

#' Create ellipse on the current plot
#'
#' @param xc
#' coordinate of center (x)
#' @param yc
#' coordinate of center (y)
#' @param a
#' major axis
#' @param b
#' minor axis
#' @param col
#' color of the ellipse line
#' @param lty
#' type of the ellipse line
#' @param ...
#' any argument suitable for \code{lines} function
#'
#' @export
ellipse <- function(xc = 0, yc = 0, a, b, col = "black", lty = 1, ...) {
   t <- seq(0, 2 * pi, 0.01)
   x <- xc + a * cos(t)
   y <- yc + b * sin(t)
   lines(x, y, type = "l", col = col, lty = lty, ...)
}

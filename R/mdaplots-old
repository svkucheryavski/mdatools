#' Check color values
#'
#' @description
#' Checks if elements of argument are valid color values
#'
#' @param palette
#' vector with possibly color values (names, RGB, etc.)
#'
mdaplot.areColors <- function(palette) {
   sapply(palette, function(x) tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE))
}

#' Format vector with numeric values
#'
#' @description
#' Format vector with values, so only significant decimal numbers are left.
#'
#' @param data
#' vector or matrix with values
#' @param round.only
#' logical, do formatting or only round the values
#' @param digits
#' how many significant digits take into account
#'
#' @details
#' Function takes into accound difference between values and the values themselves.
#'
#' @return
#' matrix with formatted values
#'
mdaplot.formatValues <- function(data, round.only = F, digits = 3) {

   # if values are not numeric - return as is
   if (!is.numeric(data[1])) return(data)

   fdata <- if (round.only) round(data, digits) else prettyNum(data, digits = digits)

   if (!is.null(dim(data))) {
      dim(fdata) <- dim(data)
      dimnames(fdata) <- dimnames(data)
   }

   return(fdata)
}

#' Calculate axes limits (common method for each axis)
#'
#' @description
#' Calculates axes limits depending on data values that have to be plotted,
#' extra plot elements that have to be shown and margins.
#'
#' @param min
#' smallest value for the axis.
#' @param max
#' largest value for the axis.
#' @param show.lines
#' numeric value - horizontal or vertical line to show for given axis.
#' @param scale
#' scale factor for limit margins (7.5% by default)
#'
#' @return
#' Returns a vector with two limits for given axis.
#'
mdaplot.getAxesLimits <- function(min, max, show.lines = NULL, scale = 0.075) {

   # if min and max are equal add 5% on each side (or just 0.05 if they are zero)
   if (min == max) {
      min <- min - 0.05 * ifelse(min == 0, 1, min)
      max <- max + 0.05 * ifelse(max == 0, 1, max)
   }

   # correct limits if some lines that have to be shown are outside the data points cloud
   if (is.numeric(show.lines)) {
      max <- max(max, show.lines, na.rm = TRUE)
      min <- min(min, show.lines, na.rm = TRUE)
   }

   # calculate margins: dx and dy
   d <- (max - min) * scale

   # define limits with margins
   lim <- c(min - d, max + d)

   # return the limits
   return(lim)
}

#' Calculate limits for x-axis.
#'
#' @description
#' Calculates limits for x-axis depending on data values that have to be plotted,
#' extra plot elements that have to be shown and margins.
#'
#' @param plot_data
#' a list with plot data and related parameters returned by `mdaplot.createPlotData()`.
#' @param xlim
#' limits provided by user
#' @param show.excluded
#' logical, show or not the excluded values
#' @param show.lines
#' logical or numeric with line coordinates to be shown on the plot.
#' @param bwd
#' if limits are computed for bar plot, this is a bar width (otherwise NULL)
#'
#' @return
#' Returns a vector with two limits.
#'
mdaplot.getXAxesLimits <- function(plot.data, xlim, show.excluded = FALSE,
   show.lines = FALSE, bwd = NULL, scale = 0.075, type = NULL) {

   # if user provided limits for x - use them strictly
   if (!is.null(xlim)) return(xlim)

   values <- plot.data$x_values
   values_excluded <- plot.data$x_values_excluded

   if (show.excluded && type == "p") {
      values <- c(values, values_excluded)
   }

   # find range for x values
   min <- min(values)
   max <- max(values)

   # find if show.lines is in use
   show.lines <- if (!is.numeric(show.lines) || is.na(show.lines[1])) NULL else show.lines[1]

   # get the limits
   lim <- mdaplot.getAxesLimits(min, max, show.lines, scale)

   # correct x axis limits and bar width if it is a bar plot
   if (!is.null(bwd) && bwd > 0) {
      bwd <- if (length(values) == 1) 2 * bwd else bwd * min(diff(values))
      lim <- c(min(values) - bwd / 2, max(values) + bwd / 2)
   }

   return(lim)
}


mdaplot.getYAxesLimits <- function(plot.data, ylim, show.excluded = FALSE,
   show.lines = FALSE, show.colorbar = FALSE, show.labels = FALSE, scale = 0.075) {

   # if user provided limits for y - use them strictly
   if (!is.null(ylim)) return(ylim)

   values <- plot.data$y_values
   values_excluded <- plot.data$y_values_excluded

   if (show.excluded) {
      values <- rbind(values, values_excluded)
   }

   # find range for x values
   min <- min(values)
   max <- max(values)

   # find if show.lines is in use
   show_y_line <- is.numeric(show.lines) && length(show.lines) == 2 && !is.na(show.lines[2])
   show.lines <- if (show_y_line) show.lines[2] else NULL

   # get the limits
   lim <- mdaplot.getAxesLimits(min, max, show.lines, scale)

   # compute margin size based on limits
   d <- -diff(lim) * scale / (1 + 2 * scale)

   # add an extra margin to y limit if colorbar must be shown
   if (show.colorbar == T) lim[2] <- lim[2] + diff(lim) * 0.075 * 3

   # add extra margins if labels must be shown
   if (show.labels == T) lim <- lim + c(-d, d)

   return(lim)
}

#' Prepare colors based on palette and opacity value
#'
#' @param palette
#' vector with main colors for current pallette
#' @param ncolors
#' number of colors to generate
#' @param opacity
#' opacity for the colors (one value or individual for each color)
#'
#' @return
#' vector with colors
#'
mdaplot.prepareColors <- function(palette, ncolors, opacity) {

   # generate colors based on color ramp and palette
   colors <- colorRampPalette(palette)(ncolors)

   # no opacity - just return colors as is
   if (is.null(opacity) || all(opacity == 1)) {
      return(colors)
   }

   # repeate opacity values for each color
   if (length(opacity) == 1) {
      opacity <- rep(opacity, ncolors)
   }

   if (length(opacity) != ncolors) {
      stop('Wrong number of values for "opacity" parameter!')
   }

   # apply opacity
   for (i in 1:ncolors) {
      colors[i] <- adjustcolor(colors[i], alpha.f = opacity[i])
   }

   return(colors)
}

#' Color values for plot elements
#'
#' @description
#' Generate vector with color values for plot objects (lines, points, bars), depending
#' on number of groups for the objects.
#'
#' @param ngroups
#' number of groups.
#' @param cgroup
#' vector of values, used for color grouping of plot points or lines.
#' @param colmap
#' which colormap to use ('default', 'gray', 'old', or user defined in form c('col1', 'col2', ...)).
#' @param opacity
#' opacity for colors (between 0 and 1)
#' @param maxsplits
#' if contenuous values are used for color gruping - how many groups to create?
#'
#' @importFrom grDevices col2rgb colorRampPalette rgb adjustcolor
#'
#' @return
#' Returns vector with generated color values
#'
#' @export
mdaplot.getColors <- function(ngroups = NULL, cgroup = NULL, colmap = "default",
   opacity = 1, maxsplits = 64) {

   # if non of the main arguments defined assume only one color is needed
   if (is.null(ngroups) && is.null(cgroup)) {
      ngroups <- 1
   }

   # define palette (if colmap has more than one value - take it as palette)
   palette <- if (length(colmap) > 1) colmap else switch(colmap,
      # gray scale based
      "gray" = c(
         "#E8E8E8", "#D6D6D6", "#C4C4C4", "#B2B2B2",
         "#9A9A9A", "#808080", "#484848", "#101010"
      ),
      # jet (like old in MATLAB)
      "jet" = c(
         "#00007F", "blue", "#007FFF", "cyan",
         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
      ),
      # old (used in mdatools in versions < 0.10.0)
      "old" = c(
         "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
         "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"
      ),
      # current default
      c(
         "#2679B2", "#1C9AA8", "#379531",
         "#EED524", "#FB7F28", "#D22C2F"
      )
   )

   if (!all(mdaplot.areColors(palette))) {
      stop("Parameter 'colmap' must contains valid color values or name of palette.")
   }

   # if grayscale palette and only one color is needed reorder pallete so the black is first
   if (all(colmap == "gray") && is.null(cgroup) && ngroups == 1) {
      palette <- rev(palette)
   }

   # if cgroup is not provided just return the colors

   if (is.null(cgroup)) {
      return(mdaplot.prepareColors(palette, ngroups, opacity))
   }

   # if cgroup is factor return vector with corresponding values
   if (is.factor(cgroup)) {
      ngroups <- length(attr(cgroup, "levels"))
      return(mdaplot.prepareColors(palette, ngroups, opacity)[as.numeric(cgroup)])
   }

   # if not split it into groups
   if (is.null(ngroups)) {
      ngroups <- length(unique(cgroup))
      ngroups <- ifelse(ngroups > maxsplits, maxsplits, ngroups)
   }

   if (ngroups > 1) {
      cgroup <- cut(as.numeric(cgroup), ngroups, include.lowest = TRUE)
   }

   out_palette <- mdaplot.prepareColors(palette, ngroups, opacity)
   colors <- out_palette[as.numeric(cgroup)]
   attr(colors, "palette") <- out_palette

   return(colors)
}

#' Plot colorbar
#'
#' @description
#' Shows a colorbar if plot has color grouping of elements (points or lines).
#'
#' @param cgroup
#' a vector with values used to make color grouping of the elements
#' @param colmap
#' a colormap to be used for color generation
#' @param lab.col
#' color for legend labels
#' @param lab.cex
#' size for legend labels
#'
mdaplot.showColorbar <- function(cgroup, colmap = "default", lab.col = "darkgray", lab.cex = 0.65) {
   # get number of levels for the cgroup

   # define if colorbar should be discrete (for factors) or not
   shift <- ifelse(is.factor(cgroup), 1, 0)

   if (!is.factor(cgroup) && length(unique(cgroup)) > 12) {
         # get colors for 8 groups based on colormap
         col <- mdaplot.getColors(ngroups = 12, colmap = colmap)
         ncol <- length(unique(col))

         # split values to intervals
         cgroupl <- levels(cut(as.vector(cgroup), ncol))

         # get left and right values for the intervals
         lvals <- as.numeric(sub("\\((.+),.*", "\\1", cgroupl))
         rvals <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", cgroupl))

         # correct issue with first element
         if (min(cgroup) != lvals[1]) {
            lvals[1] <- min(cgroup)
         }

         # combine values and define matrix for labels
         vals <- c(lvals, rvals[ncol])
         labels <- matrix(0, ncol = 2, nrow = ncol + 1)
   } else {
      if (!is.factor(cgroup)) {
         cgroup <- factor(cgroup)
      }

      nlevels <- length(attr(cgroup, "levels"))

      # no splitting is needed, just use factors as labels
      col <- mdaplot.getColors(ngroups = nlevels, colmap = colmap)
      ncol <- length(unique(col))
      vals <- levels(cgroup)
      labels <- matrix(0, ncol = 2, nrow = ncol)
   }

   # use formatted values as rownames for labels matrix
   rownames(labels) <- mdaplot.formatValues(vals)

   # get size of the plotting area and calculate size for color bar elements
   lim <- par("usr")

   dx <- lim[2] - lim[1]
   dy <- lim[4] - lim[3]

   w <- (dx * 0.8) / ncol
   h <- dy * 0.015
   shift <- shift * w * 0.02 # 2 percent of segment width

   x <- lim[1] + dx * 0.1
   y <- lim[4] - (h + 0.1 * h);

   # show colorbar and define coordinates for labels
   for (i in 1:ncol) {
      rect(shift + x + w * (i - 1), y, x + w * i, y - h, col = col[i], border = NA)
      labels[i, ] <- c(x + w * (i - 1), y - h)
   }

   # add last value or shift coordinates if labels shall be centered
   if (nrow(labels) > i)
      labels[i + 1, ] <- c(x + w * i, y - h)
   else
      labels[, 1] <- labels[, 1] + w / 2

   # show labels for colorbar regions
   mdaplot.showLabels(labels[, 1], labels[, 2], labels = rownames(labels), pos = 1, col = lab.col,
      cex = lab.cex)
}

#' Plot legend
#'
#' @description
#' Shows a legend for plot elements or their groups.
#'
#' @param legend
#' vector with text elements for the legend items
#' @param col
#' vector with color values for the legend items
#' @param pch
#' vector with marker symbols for the legend items
#' @param lty
#' vector with line types for the legend items
#' @param lwd
#' vector with line width values for the legend items
#' @param cex
#' vector with cex factor for the points
#' @param bty
#' border type for the legend
#' @param position
#' legend position ("topright", "topleft', "bottomright", "bottomleft", "top", "bottom")
#' @param plot
#' logical, show legend or just calculate and return its size
#'
mdaplot.showLegend <- function(legend, col, pch = NULL, lty = NULL, lwd = NULL, cex = 1,
                              bty = "o", position = "topright", plot = TRUE) {
   # which positions need multiple columns
   onecolpos <- c("topright", "topleft", "bottomright", "bottomleft")
   multcolpos <- c("top", "bottom", "right", "left")

   if (!(position %in% c(onecolpos, multcolpos))) {
      stop("Wrong values for 'legend.position' argument!")
   }

   # compute number of columns
   ncol <- if (position %in% onecolpos) 1 else length(legend)

   # calculate inset values depending on a ration between width and height of a plot
   lim <- par("plt")

   dx <- lim[2] - lim[1]
   dy <- lim[4] - lim[3]

   inset <- c(0.02, 0.02 * (dx / dy))

   # show legend
   legend(position, legend, col = col, pch = pch, lty = lty, pt.cex = cex, lwd = lwd,
          cex = 0.85, plot = plot, inset = inset, bg = "white", box.lwd = 0.75, box.col = "gray",
          ncol = ncol)

}

#' Plot labels
#' Shows labels for data elements (points, bars) on a plot.
#'
#' @param x_values
#' a vector with x-values
#' @param y_values
#' a vector with y-values
#' @param labels
#' a vector with labels
#' @param pos
#' position of the labels relative to the points
#' @param cex
#' size of the labels text
#' @param col
#' color of the labels text
#' @param type
#' type of the plot
#'
#' @details
#' Rownames of matrix \code{data} are used as labels. If matrix has no rownames, row numbers
#' will be used instead.
#'
mdaplot.showLabels <- function(x_values, y_values, labels, pos = 3,
   cex = 0.65, col = "darkgray", type = "") {

   if (length(labels) == 0 || length(x_values) == 0 || length(y_values) == 0) {
      warning("No labels available.")
      return()
   }

   if (type %in% c("h", "e")) {
      y_values <- y_values[1, ]
   }

   if (type %in% c("l", "b")) {
      y_values <- apply(y_values, 2, max)
   }

   # prepare labels
   x_values <- as.numeric(x_values)
   y_values <- as.numeric(y_values)
   labels <- mdaplot.formatValues(labels)

   # if x or y values contains NA do nothing
   if (any(is.nan(x_values)) || any(is.nan(y_values))) {
      return()
   }

   if (!is.null(type) && type == "h") {
      # show labels properly for bars with positive and negative values
      text(x_values, y_values, labels, cex = cex, pos = ifelse(y_values < 0, 1, 3), col = col)
      return()
   }

   # default way to show the labels
   text(x_values, y_values, labels, cex = cex, pos = pos, col = col)
}

#' Plot lines
#'
#' @description
#' Shows horisontal and vertical lines on a plot.
#'
#' @param point
#' vector with two values: x coordinate for vertical point y for horizontal
#' @param lty
#' line type
#' @param lwd
#' line width
#' @param col
#' color of lines
#'
#' @details
#' If it is needed to show only one line, the other coordinate shall be set to NA.
#'
mdaplot.showLines <- function(point, lty = 2, lwd = 0.75, col = rgb(0.2, 0.2, 0.2)) {

   if (!is.na(point[2])) {
      abline(h = point[2], lty = lty, lwd = lwd, col = col)
   }

   if (!is.na(point[1])) {
      abline(v = point[1], lty = lty, lwd = lwd, col = col)
   }
}

#' Regression line for data points
#'
#' @description
#' Shows linear fit line for data points.
#'
#' @param data
#' data values
#' @param lty
#' line type
#' @param lwd
#' line width
#' @param colmap
#' color map
#' @param col
#' color of lines
#'
mdaplot.showRegressionLine <- function(data, lty = 1, lwd = 1, colmap = "default", col = NULL) {

   if (!is.list(data)) {
      data <- list(data)
   }

   ngroups <- length(data)

   if (is.null(col)) {
      col <- mdaplot.getColors(ngroups = ngroups, colmap = colmap)
   }

   if (length(col) == 1) {
      col <- rep(col, ngroups)
   }

   for (i in 1:ngroups) {
      abline(lm(data[[i]][, 2] ~ data[[i]][, 1]), lty = lty, lwd = lwd, col = col[i])
   }
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
#' @export
mdaplot.getConfidenceEllipse <- function(points, conf.level = 0.95, n = 100) {

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
#' @export
mdaplot.getConvexHull <- function(points) {
   ch_ind <- chull(points)
   ch_ind <- c(ch_ind, ch_ind[1])
   return(points[ch_ind, ])
}


#' Create dataset to show as plot
#'
#' @param data
#' dataset (vector, matrix or data frame)
#' @param type
#' type of the plot
#' @param main
#' main plot title provided by user (if any)
#' @param xlab
#' label for x-axis provided by user (if any)
#' @param ylab
#' label for y-axis provided by user (if any)
#' @param xlim
#' limits for x-axis provided by user
#' @param ylim
#' limits for y-axis provided by user
#' @param bwd
#' width of bars provided by user
#' @param show.excluded
#' logical, show or not excluded values
#' @param show.colorbar
#' logical, show or not colorbar (needed to compute limits correctly)
#' @param show.labels
#' logical, show or not labels for data points
#' @param show.lines
#' logical, show or not lines
#' @param show.axes
#' logical, show or not axes
#' @param labels
#' vector with labels
#'
#' @export
mdaplot.createPlotData <- function(data, type, main, xlab, ylab, xlim, ylim, bwd, show.excluded,
   show.colorbar, show.labels, show.lines, show.axes, labels) {

   # check plot type

   valid.types <- c("p", "l", "b", "h", "e", "d")
   if (!(type %in% valid.types)) {
      stop("Wrong plot type!")
   }

   # get attributes and prepare dataset
   data_attrs <- attributes(data)
   plot_data <- mdaplot.prepareDataForPlot(data, type)

   # check if density plot should be shown
   plot_data$density <- FALSE
   if (type == "d") {
      if (ncol(data) < 2) stop("Dataset with at least two columns required for density plot.")
      type <- "p"
      plot_data$density <- TRUE
   }

   # split plot data to x and y values
   plot_data <- mdaplot.splitPlotData(plot_data, type, xlab, ylab)

   # process excluded rows if they have to be shown
   if (show.excluded) {
      plot_data <- mdaplot.processExcludedRows(plot_data, type)
   }

   # create labels if necessary
   if (show.labels) {
      plot_data <- mdaplot.prepareDataLabels(plot_data, type, labels)
   }

   # define name for the plot series
   plot_data$name <- if (type == "h") rownames(plot_data$data)[1] else data_attrs[["name"]]
   if (!is.null(main)) plot_data$name <- main

   # compute limits for axes
   plot_data$xlim <- mdaplot.getXAxesLimits(
      plot_data,
      xlim,
      show.lines = show.lines,
      show.excluded = show.excluded,
      bwd = if (type == "h") bwd else NULL,
      type = type
   )

   plot_data$ylim <- mdaplot.getYAxesLimits(
      plot_data,
      ylim,
      show.lines = show.lines,
      show.excluded = show.excluded,
      show.colorbar = show.colorbar
   )


   return(plot_data)
}


#' Split dataset to x and y values depending on plot type
#'
#' @param plot_data
#' list with data and parameters prepared for plot by `mdaplot.prepareDataForPlot()` function
#' @param type
#' type of plot
#' @param xlab
#' label for x-axis provided by user (if any)
#' @param ylab
#' label for y-axis provided by user (if any)
#'
#' @export
mdaplot.splitPlotData <- function(plot_data, type, xlab = NULL, ylab = NULL) {

   # shortcuts to some of parameters
   data <- plot_data$data
   data_attrs <- plot_data$attrs
   excluded_cols <- plot_data$excluded_cols
   xaxis_name <- if (!is.null(xlab)) xlab else data_attrs$xaxis.name
   yaxis_name <- if (!is.null(ylab)) ylab else data_attrs$yaxis.name

   if (type == "p") {
      # scatter plot
      plot_data$x_values <- data[, 1]
      attr(plot_data$x_values, "name") <- if (is.null(xlab)) colnames(data)[1] else xlab
      plot_data$y_values <- data[, 2, drop = F]
      attr(plot_data$y_values, "name") <- if (is.null(ylab)) colnames(data)[2] else ylab
      return(plot_data)
   }

   if (type == "e") {
      # errorbar plor
      data <- rbind(
         data[1, ],
         data[1, ] - data[2, ],
         data[1, ] + if (nrow(data) > 2) data[3, ] else data[2, ]
      )
   }

   # prepare x-axis values for other types of plots
   x_values <- seq_len(ncol(data))
   if (!is.null(data_attrs$xaxis.values)) {
      x_values <- data_attrs$xaxis.values
      if (length(excluded_cols) > 0) x_values[-excluded_cols]
   }

   # assign names to x values
   if (is.null(names(x_values))) names(x_values) <- colnames(data)

   # add x values to the plot data and assign name to x-variable
   plot_data$x_values <- x_values
   attr(plot_data$x_values, "name") <- if (is.null(xaxis_name)) "Variables" else xaxis_name

   # by default take all data values as y and assign name to the x-variable
   plot_data$y_values <- data
   attr(plot_data$y_values, "name") <- if (is.null(yaxis_name)) "" else yaxis_name

   return(plot_data)
}

#' Prepare x and y values for excluded rows
#'
#' @param plot_data
#' list with data and parameters prepared for plot by `mdaplot.prepareDataForPlot()` function
#' @param type
#' type of plot
#'
#' @export
mdaplot.processExcludedRows <- function(plot_data, type) {

   # shortcuts for some parameters
   excluded_data <- plot_data$excluded_data

   # if nothing to show - return
   if (is.null(excluded_data) || (nrow(excluded_data) == 0)) {
      return(plot_data)
   }

   if (type == "p") {
      plot_data$x_values_excluded <- excluded_data[, 1]
      plot_data$y_values_excluded <- excluded_data[, 2, drop = F]
      return(plot_data)
   }

   plot_data$x_values_excluded <- plot_data$x_values
   plot_data$y_values_excluded <- excluded_data
   return(plot_data)
}

#' Check if data argument looks correct for creatins plots, correct if needed
#'
#' @param data
#' dataset provided to mdaplot function
#' @param type
#' type of plot
#'
#' @export
mdaplot.prepareDataForPlot <- function(data, type) {

   # save data arguments
   data_attrs <- attributes(data)

   # if it is vector without dimension - make a matrix
   if (is.null(dim(data))) {
      data <- if (type == "p") as.matrix(data) else t(as.matrix(data))
   }

   # at least two rows (y values and size of error bars) are needed for errorbar plot
   if (type == "e" && nrow(data) < 2) {
      stop("Errorbar plot requires dataset with at least two rows!")
   }

   if (type == "e" && length(data_attrs$exclrows) > 0) {
      stop("Errobar plot can not be made for data with excluded rows.")
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

   # handle excluded columns if any
   excluded_cols <- data_attrs$exclcols
   col_indices <- seq_len(ncol(data))
   if (length(excluded_cols) > 0) {
      excluded_cols <- mda.getexclind(excluded_cols, colnames(data), ncol(data))
      data <- data[, -excluded_cols, drop = F]
      col_indices <- col_indices[-excluded_cols]
   }

   # check if data still has some columns
   if (is.null(ncol(data))) {
      stop("No columns left when excluded hidden values.")
   }

   # define yaxis values
   if (is.null(data_attrs$yaxis.values)) {
      data_attrs$yaxis.values <- seq_len(nrow(data))
   }

   # add columns with row numbers or yaxis values to the data if necessary
   if (type == "p" && ncol(data) == 1) {
      rownames <- rownames(data)
      new_column <- data_attrs$yaxis.values
      data <- mda.cbind(new_column, data)
      colnames(data)[1] <- if (is.null(data_attrs$yaxis.name)) "Objects" else data_attrs$yaxis.name
      rownames(data) <- rownames
   }

   # handle excluded rows if necessary
   row_indices <- seq_len(nrow(data))
   excluded_data <- NULL
   excluded_rows <- data_attrs$exclrows
   if (length(excluded_rows) > 0) {
      excluded_rows <- mda.getexclind(excluded_rows, rownames(data), nrow(data))
      excluded_data <- data[excluded_rows, , drop = F]
      data <- data[-excluded_rows, , drop = F]
      row_indices <- row_indices[-excluded_rows]
   }

   return(
      list(
         data = data,
         excluded_data = excluded_data,
         excluded_rows = excluded_rows,
         excluded_cols = excluded_cols,
         row_indices = row_indices,
         col_indices = col_indices,
         attrs = data_attrs
      )
   )
}


#' Create data labels as values
#'
#' @param plot_data
#' list with data and parameters prepared for plot by `mdaplot.prepareDataForPlot()` function
#' @param type
#' type of plot
#'
mdaplot.getLabelsAsValues <- function(plot_data, type) {
   y_values <- plot_data$y_values
   y_values_excluded <- plot_data$y_values_excluded

   if (type %in% c("l", "b")) {
      plot_data$labels <- apply(y_values, 2, max)
      plot_data$labels_excluded <- if (!is.null(y_values_excluded)) apply(y_values_excluded, 2, max)
      return(plot_data)
   }

   if (type %in% c("h", "e")) {
      plot_data$labels <- y_values[1, ]
      plot_data$labels_excluded <- y_values_excluded[1, ]
      return(plot_data)
   }

   # otherwise prepare labels as for scatter plot
   plot_data$labels <- y_values
   plot_data$labels_excluded <- y_values_excluded
   return(plot_data)
}

#' Create data labels as indices
#'
#' @param plot_data
#' list with data and parameters prepared for plot by `mdaplot.prepareDataForPlot()` function
#' @param type
#' type of plot
#'
mdaplot.getLabelsAsIndices <- function(plot_data, type) {
   excluded_rows <- plot_data$excluded_rows

   if (type == "p") {
      plot_data$labels <- plot_data$row_indices
      plot_data$labels_excluded <- excluded_rows
      return(plot_data)
   }

   plot_data$labels <- plot_data$col_indices
   plot_data$labels_excluded <- NULL
   return(plot_data)
}

#' Create a vector with labels for plot data series
#'
#' @param plot_data
#' list with data and parameters prepared for plot by `mdaplot.prepareDataForPlot()` function
#' @param type
#' type of plot
#' @param labels
#' vector or type of labels to show
#'
#' @description
#' For scatter plots labels correspond to rows of the data (names, values, indices, etc.). For
#' non-scatter plots labels correspond to the columns (names, indices or max value for each column)
#'
#' @export
mdaplot.prepareDataLabels <- function(plot_data, type, labels = NULL) {
   y_values <- plot_data$y_values
   x_values <- plot_data$x_values

   # if labels are not specified - use names
   if (is.null(labels)) labels <- "names"

   # if user provided labels - use them
   if (length(labels) > 1) {
      excluded_cols <- plot_data$excluded_cols
      excluded_rows <- plot_data$excluded_rows

      # get values and excluded values (rows or cols depending on plot type)
      n_values <- if (type == "p") nrow(y_values) else length(x_values)
      excluded_values <- if (type == "p") excluded_rows else excluded_cols

      # check that labels were specified for all values
      if (length(labels) != (n_values + length(excluded_values))) {
         stop("Labels must be provided for all values (also ones which are excluded).")
      }

      if (is.null(excluded_values)) {
         plot_data$labels_excluded <- NULL
         plot_data$labels <- labels
         return(plot_data)
      }

      plot_data$labels_excluded <- labels[excluded_values]
      plot_data$labels <- labels[-excluded_values]
      return(plot_data)
   }

   # if labels must be values use y-values for that
   if (labels == "values") {
      return(mdaplot.getLabelsAsValues(plot_data, type))
   }

   # if labels must be indices use row or column indices
   if (labels == "indices") {
      return(mdaplot.getLabelsAsIndices(plot_data, type))
   }

   # if nothing above works - use names as labels
   x_values_excluded <- plot_data$x_values_excluded
   plot_data$labels <- names(x_values)
   plot_data$labels_excluded <- if (type == "p") names(x_values_excluded) else NULL
   return(plot_data)
}



#' Prepare xticks for plot
#'
#' @param xticks
#' xticks provided by user (if any)
#' @param x_values
#' x values for the plot data object
#' @param xlim
#' limits for x axis
#' @param type
#' type of the plot
mdaplot.prepareXTicks <- function(xticks, x_values, xlim, type) {

   if (!is.null(xticks)) return(xticks)
   if (type != "p" && length(x_values) == 1) return(1)
   return(axisTicks(xlim, log = FALSE))
}

#' Prepare yticks for plot
#'
#' @param yticks
#' yticks provided by user (if any)
#' @param y_values
#' y values for the plot data object
#' @param ylim
#' limits for y axis
#' @param type
#' type of the plot
mdaplot.prepareYTicks <- function(yticks, y_values, ylim, type) {

   if (!is.null(yticks)) return(yticks)
   if (type != "p" && length(y_values) == 1) return(1)

   return(axisTicks(ylim, log = FALSE))
}

#' Prepare xticklabels for plot
#'
#' @param xticklabels
#' xticklables provided by user (if any)
#' @param xticks
#' xticks (provided or computed)
#' @param excluded_cols
#' columns excluded from plot data (if any)
#'
mdaplot.prepareXTickLabels <- function(xticklabels, xticks, excluded_cols) {

   if (is.null(xticklabels)) return(TRUE)
   if (is.null(xticks)) stop("You need to specify both 'xticklabels' and 'xticks'")

   # if xticklabels were provided - remove excluded columns if any and check the length
   if (!is.null(excluded_cols)) xticklabels <- xticklabels[-excluded_cols]
   if (length(xticks) != length(xticklabels)) {
      stop('Number of elements in "xticks" and "xticklabels" should be the same')
   }

   return(xticklabels)
}

#' Prepare yticklabels for plot
#'
#' @param yticklabels
#' yticklables provided by user (if any)
#' @param yticks
#' yticks (provided or computed)
#' @param excluded_rows
#' rows excluded from plot data (if any)
#'
mdaplot.prepareYTickLabels <- function(yticklabels, yticks, excluded_rows) {

   if (is.null(yticklabels)) return(TRUE)
   if (is.null(yticks)) stop("You need to specify both 'yticklabels' and 'yticks'")

   # if yticklabels were provided - remove excluded rows if any and check the length
   if (length(yticks) != length(yticklabels)) {
      stop('Number of elements in "yticks" and "yticklabels" should be the same')
   }

   return(yticklabels)
}

mdaplot.prepareDataForGPlots <- function(data, type, groupby) {

   if (is.list(data) && !is.data.frame(data)) return(data)

   if (is.null(groupby)) {
      # take every row of matrix or data frame as separate group

      if (!all(type %in% c("h", "l", "b"))) {
         stop("Group plot with matrix or data frame can be made only for types 'h', 'l' and 'b'.")
      }

      # split data into a list of subsets for each group
      data_list <- list()
      for (i in seq_len(nrow(data))) {
         data_list[[rownames(data)[i]]] <- mda.subset(data, subset = i)
      }

      # redefine the data with list
      return(data_list)
   }

   # if groupby is provided - use it to split rows into groups

   ## check that groupby is a factor or data frame with factor columns
   ## TODO: simplify the three lines below
   is_groupby_factor <-
      (is.data.frame(groupby) && all(unlist(lapply(groupby, is.factor)))) ||
      is.factor(groupby)

   if (!is_groupby_factor) {
      stop("Parameter 'groupby' should be a factor or data frame with several factors.")
   }

   attrs <- mda.getattr(data)
   data <- as.data.frame(data)

   # in this case if labels = indices generate labels for each case
   data$row.ind <- seq_len(nrow(data))
   data_list <- split(data, groupby)

   for (i in seq_len(length(data_list))) {
      row.ind <- data_list[[i]]$row.ind
      data_list[[i]] <- subset(data_list[[i]], select = -row.ind)
      data_list[[i]] <- mda.setattr(data_list[[i]], attrs)
      attr(data_list[[i]], "exclrows") <- which(row.ind %in% attrs$exclrows)
      attr(data_list[[i]], "labels") <- row.ind
   }

   return(data_list)
}


#' Create axes plane
#'
#' @description
#' Creates an empty axes plane for given parameters
#'
#' @param xticklabels
#' labels for x ticks
#' @param yticklabels
#' labels for y ticks
#' @param xticks
#' values for x ticks
#' @param yticks
#' values for y ticks
#' @param xlim
#' vector with limits for x axis
#' @param ylim
#' vector with limits for y axis
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param xlas
#' orientation of xticklabels
#' @param ylas
#' orientation of yticklabels
#' @param show.grid
#' logical, show or not axes grid
#' @param grid.lwd
#' line thinckness (width) for the grid
#' @param grid.col
#' line color for the grid
#'
mdaplot.plotAxes <- function(xticklabels = NULL, yticklabels = NULL,
   xlim = xlim, ylim = ylim, xticks = NULL, yticks = NULL, main = NULL, xlab = NULL, ylab = NULL,
   xlas = 0, ylas = 0, show.grid = TRUE, grid.lwd = 0.5, grid.col = "lightgray") {

   # make plot without ticks
   plot(0, 0, type = "n", main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
        xaxt = "n", yaxt = "n")

   # generate x and y ticks
   if (is.null(xticks)) xticks <- axisTicks(xlim, log = FALSE)
   if (is.null(yticks)) yticks <- axisTicks(ylim, log = FALSE)

   # show x-axis
   if (is.null(xticklabels)) xticklabels <- TRUE
   axis(1, at = xticks, labels = xticklabels, las = xlas)

   # show y-axis
   if (is.null(yticklabels)) yticklabels <- TRUE
   axis(2, at = yticks, labels = yticklabels, las = ylas)

   # show grid if needed
   if (show.grid) {
      grid(lwd = grid.lwd, col = grid.col)
   }
}

#' Make density plot using hexagon binning
#'
#' @param x
#' vector with x coordinates of data points
#' @param y
#' vector with y coordinates of data points
#' @param nbins
#' number of bins to split x and y axes into
#' @param colmap
#' colormap to use for showing points density
#'
#' @export
mdaplot.plotDensity <- function(x, y, nbins = 60, colmap = "default") {

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
#' @param plot_data
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
#' p <- mdaplot(people, cgroup = group)
#' mdaplot.plotConfidenceEllipse(p, conf.level = 0.90, opacity = 0.2)
#'
#' @export
mdaplot.plotConfidenceEllipse <- function(plot_data, conf.level = 0.95,
   lwd = 1, lty = 1, opacity = 0) {

   mdaplot.plotPointsShape(
      plot_data,
      lwd = lwd,
      lty = lty,
      opacity = opacity,
      shape_function = mdaplot.getConfidenceEllipse,
      conf.level = conf.level
   )
}


#' Add convex hull for groups of points on scatter plot
#'
#' @param plot_data
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
#' p <- mdaplot(people, cgroup = group)
#' mdaplot.plotConvexHull(p)
#'
#' @importFrom grDevices chull
#' @importFrom graphics polygon
#'
#' @export
mdaplot.plotConvexHull <- function(plot_data, lwd = 1, lty = 1, opacity = 0) {
   mdaplot.plotPointsShape(
      plot_data,
      lwd = lwd,
      lty = lty,
      opacity = opacity,
      shape_function = mdaplot.getConvexHull
   )
}


#' Add confidence ellipse or convex hull for group of points
#'
#' @param plot_data
#' plot data returned by function `mdaplot()`
#' @param lwd
#' thickness of line used to show the hull
#' @param lty
#' type of line used to show the hull
#' @param conf.level
#' confidence level to make the ellipse for (between 0 and 1)
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
mdaplot.plotPointsShape <- function(plot_data, lwd, lty, opacity, shape_function, ...) {

   x <- plot_data$x_values
   y <- plot_data$y_values
   cgroup <- plot_data$cgroup
   type <- plot_data$type

   if (type != "p") {
      stop("Shape can be added only to scatter plots.")
   }

   if (!is.factor(cgroup)) {
      stop("Parameter 'cgroup' must be a factor if you want to show points shape.")
   }

   if (length(unique(cgroup)) != length(unique(plot_data$col))) {
      stop("Number of colors should be the same as number of levels in parameter 'cgroup'.")
   }

   col <- split(plot_data$col, f = plot_data$cgroup, drop = TRUE)
   d <- split(data.frame(x, y), f = cgroup, drop = TRUE)

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


#' Show bars on axes
#'
#' @description
#' Shows bars (bar plot) on predefined axes
#'
#' @param x
#' vector with x values (centers of bars)
#' @param y
#' vector with y values (height of bars)
#' @param col
#' colors of the bars
#' @param bwd
#' width of the bars (as a ratio for max width)
#' @param border
#' color of bar edges
#' @param force.x.values
#' vector with corrected x-values for a bar plot (do not specify this manually).
#'
mdaplot.plotBars <- function(x, y, col = NULL, bwd = 0.8, border = NA, force.x.values = NA) {
   # correct x_values if they were forced by bwd
   if (is.numeric(force.x.values)) {
      x <- x - bwd / 2 + (force.x.values[1] - 0.5) * bwd / force.x.values[2]
      bwd <- bwd / force.x.values[2]
   }

   if (length(bwd) == 1) {
      bwd <- matrix(bwd, ncol = length(x))
   }

   if (length(col) != length(y)) {
      col <- rep(col, length.out = length(y))
   }

   rect(x - bwd / 2, 0, x + bwd / 2, y[1, ], col = col, border = border)
}


#' Show error bars on a plot
#'
#' @description
#' Shows error bars (errorbar plot) on predefined axes
#'
#' @param x
#' vector with x values
#' @param y
#' matrix with y values (matrix with three or two rows)
#' @param col
#' color for the error bars
#' @param pch
#' marker symbol for the plot
#' @param cex
#' cex factor for the marker
#'
mdaplot.plotErrorbars <- function(x, y, col = NULL, pch = 16, cex = 1) {

   e2 <- (max(x) - min(x)) / 50
   e1 <- (max(x) - min(x)) / (length(x) * 5)
   e <- min(e1, e2)
   if (e == 0) x[1] * 0.05

   segments(x, y[2, ], x, y[3, ], col = col)
   segments(x - e, y[2, ], x + e, y[2, ], col = col)
   segments(x - e, y[3, ], x + e, y[3, ], col = col)
   points(x, y[1, ], col = col, pch = pch, cex = cex)
}

#' Show set of points on a plot
#'
#' @description
#' Show scatter plot on predefined axes
#'
#' @param x
#' vector with x values
#' @param y
#' matrix with y values
#' @param x.excluded
#' vector with x values for excluded rows
#' @param y.excluded
#' matrix with y values for excluded rows
#' @param pch.colinv
#' allows to swap values for `col` and `bg` for scatter plots with `pch` valyes from 21 to 25.
#' @param pch
#' a character for markers (same as \code{plot} parameter).
#' @param col
#' a color for markers or lines (same as \code{plot} parameter).
#' @param bg
#' background color for scatter plots wich `pch=21:25`.
#' @param col.excluded
#' color for the excluded points.
#' @param ...
#' other arguments for function `points()`.
#'
mdaplot.plotScatter <- function(x, y, x.excluded, y.excluded, pch.colinv, pch, col,
   bg, col.excluded, ...) {

   if (pch.colinv) {
      # switch colors for main and background in case pch=21:25
      bg.old <- bg
      bg <- col
      col <- bg.old
   }

   # show main set of points
   points(x, y, col = col, bg = bg, pch = pch, ...)

   # show excluded points if any
   if (!is.null(x.excluded) && !is.null(y.excluded)) {
      points(x.excluded, y.excluded, col = col.excluded, pch = pch, ...)
   }
}

#' Show set of lines on a plot
#'
#' @description
#' Show line plot on predefined axes
#'
#' @param x
#' vector with x values
#' @param y
#' matrix with y values
#' @param x.excluded
#' vector with x values for excluded rows
#' @param y.excluded
#' matrix with y values for excluded rows
#' @param type
#' type of the plot ("l" or "b").
#' @param pch
#' a character for markers (same as \code{plot} parameter).
#' @param col
#' a color for markers or lines (same as \code{plot} parameter).
#' @param col.excluded
#' color for the excluded lines.
#' @param ...
#' other arguments for function `points()`.
#'
mdaplot.plotLines <- function(x, y,  x.excluded, y.excluded, type, pch, col, col.excluded, ...) {

   # show main set of lines
   matlines(x, t(y), type = type, col = col, pch = pch, ...)

   # show excluded rows
   if (!is.null(x.excluded) && !is.null(y.excluded)) {
      matlines(x.excluded, t(y.excluded), type = type, col = col.excluded, pch = pch, ...)
   }
}

#' Plotting function for a single set of objects
#'
#' @description
#' \code{mdaplot} is used to make different kinds of plot for one set of data objects.
#'
#' @param data
#' a vector, matrix or a data.frame with data values.
#' @param plot.data
#' a list of parameters and values obtained after preprocessing of original data provided
#' by a user (if NULL it will be created automatically).
#' @param type
#' type of the plot ('p', 'l', 'b', 'h', 'e', 'i').
#' @param cgroup
#' a vector with values to use for make color groups.
#' @param colmap
#' a colormap to use for coloring the plot items.
#' @param pch
#' a character for markers (same as \code{plot} parameter).
#' @param col
#' a color for markers or lines (same as \code{plot} parameter).
#' @param bg
#' background color for scatter plots wich `pch=21:25`.
#' @param lty
#' the line type (same as \code{plot} parameter).
#' @param lwd
#' the line width (thickness) (same as \code{plot} parameter).
#' @param cex
#' the cex factor for markers (same as \code{plot} parameter).
#' @param bwd
#' a width of a bar as a percent of a maximum space available for each bar.
#' @param border
#' color for border of bars (if barplot is used)
#' @param xlim
#' limits for the x axis (if NULL, will be calculated automatically).
#' @param ylim
#' limits for the y axis (if NULL, will be calculated automatically).
#' @param xlab
#' a title for the x axis (same as \code{plot} parameter).
#' @param ylab
#' a title for the y axis (same as \code{plot} parameter).
#' @param main
#' an overall title for the plot (same as \code{plot} parameter).
#' @param labels
#' a vector with text labels for data points or one of the following: 'names', 'indices', 'values'.
#' @param show.labels
#' logical, show or not labels for the data objects.
#' @param show.colorbar
#' logical, show or not colorbar legend if color grouping is on.
#' @param show.lines
#' vector with two coordinates (x, y) to show horizontal and vertical line cross the point.
#' @param show.grid
#' logical, show or not a grid for the plot.
#' @param grid.lwd
#' line thinckness (width) for the grid.
#' @param grid.col
#' line color for the grid.
#' @param show.axes
#' logical, make a normal plot or show only elements (markers, lines, bars) without axes.
#' @param xticks
#' values for x ticks.
#' @param yticks
#' values for y ticks.
#' @param xticklabels
#' labels for x ticks.
#' @param yticklabels
#' labels for y ticks.
#' @param xlas
#' orientation of xticklabels.
#' @param ylas
#' orientation of yticklabels.
#' @param lab.col
#' color for data point labels.
#' @param lab.cex
#' size for data point labels.
#' @param show.excluded
#' logical, show or hide rows marked as excluded (attribute `exclrows`).
#' @param col.excluded
#' color for the excluded objects (rows).
#' @param nbins
#' if scatter density plot is shown, number of segments to split the plot area into.
#' (see also ?smoothScatter)
#' @param colramp
#' Colramp function for density scatter plot.
#' @param force.x.values
#' vector with corrected x-values for a bar plot (do not specify this manually).
#' @param opacity
#' opacity for plot colors (value between 0 and 1).
#' @param pch.colinv
#' allows to swap values for `col` and `bg` for scatter plots with `pch` valyes from 21 to 25.
#' @param ...
#' other plotting arguments.
#'
#' @details
#' Most of the parameters are similar to what are used with standard \code{plot} function. The
#' differences are described below.
#'
#' The function makes a plot of one set of objects. It can be a set of points (scatter plot),
#' bars, lines, scatter-lines, errorbars og an image. The data is organized as a data frame,
#' matrix or vector. For scatter and only first two columns will be used, for bar plot only
#' values from the first row. It is recommended to use \code{\link{mda.subset}} method if plot
#' should be made only for a subset of the data, especially if you have any excluded rows or
#' columns or other special attributed, described in the Bookdown tutorial.
#'
#' If data is a data frame and contains one or more factors, they will be converted to a dummy
#' variables (using function \code{\link{mda.df2mat}}) and appears at the end (last columns) if
#' line or bar plot is selected.
#'
#' The function allows to colorize lines and points according to values of a parameter
#' \code{cgroup}. The parameter must be a vector with the same elements as number of objects (rows)
#' in the data. The values are divided into up to eight intervals and for each interval a
#' particular color from a selected color scheme is assigned. Parameter \code{show.colorbar}
#' allows to turn off and on a color bar legend for this option.
#'
#' The used color scheme is defined by the \code{colmap} parameter. The default scheme is based
#' on color brewer (colorbrewer2.org) diverging scheme with eight colors. There is also a gray
#' scheme (\code{colmap = 'gray'}) and user can define its own just by specifing the needed
#' sequence of colors (e.g. \code{colmap = c('red', 'yellow', 'green')}, two colors is minimum).
#' The scheme will then be generated automatically as a gradient among the colors.
#'
#' Besides that the function allows to change tick values and corresponding tick labels for x and
#' y axis, see Bookdown tutorial for more details.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' \code{\link{mdaplotg}} - to make plots for several sets of data objects (groups of objects).
#'
#' @examples
#' # See all examples in the tutorial.
#'
#' @export
mdaplot <- function(data = NULL, plot.data = NULL, type = "p",
      pch = 16, col = NULL, bg = par("bg"), bwd = 0.8, lty = 1, lwd = 1, cex = 1, border = NA,
      cgroup = NULL, xlim = NULL, ylim = NULL, colmap = "default", labels = NULL,
      main = NULL, xlab = NULL, ylab = NULL, show.labels = F,
      show.colorbar = !is.null(cgroup), show.lines = FALSE, show.grid = TRUE, grid.lwd = 0.5,
      grid.col = "lightgray", show.axes = TRUE, xticks = NULL, yticks = NULL,
      xticklabels = NULL, yticklabels = NULL,
      xlas = 0, ylas = 0, lab.col = "darkgray", lab.cex = 0.65,
      show.excluded = FALSE, col.excluded = "#E0E0E0", nbins = 60,
      colramp = mdaplot.getColors, force.x.values = NA, opacity = 1,
      pch.colinv = FALSE, ...) {

   if (is.null(plot.data)) {
      # get the data for plot
      plot.data <- mdaplot.createPlotData(
         data = data, type = type, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         bwd = bwd, show.excluded = show.excluded, show.colorbar = show.colorbar,
         show.labels = show.labels, show.lines = show.lines, show.axes = show.axes, labels = labels
      )
   }

   # get values from the plot data list to make code cleaner
   plot.data$x_values -> x_values
   plot.data$y_values -> y_values
   plot.data$x_values_excluded -> x_values_excluded
   plot.data$y_values_excluded -> y_values_excluded
   plot.data$excluded_rows -> excluded_rows
   plot.data$excluded_cols -> excluded_cols
   plot.data$density -> density
   plot.data$labels -> labels
   plot.data$labels_excluded -> labels_excluded
   plot.data$xlim -> xlim
   plot.data$ylim -> ylim

   # show axes if needed
   if (show.axes) {

      # check and prepare xticklabels
      xticklabels <- mdaplot.prepareXTickLabels(xticklabels, xticks, excluded_cols)
      xticks <- mdaplot.prepareXTicks(xticks, x_values, xlim, type)

      # check and prepare yticklabels
      yticklabels <- mdaplot.prepareYTickLabels(yticklabels, yticks, excluded_rows)
      yticks <- mdaplot.prepareYTicks(yticks, y_values, ylim, type)

      main <- plot.data$name
      xlab <- attr(x_values, "name")
      ylab <- attr(y_values, "name")

      # make an empty plot with proper limits and axis labels
      mdaplot.plotAxes(xticklabels = xticklabels, yticklabels = yticklabels, xticks = xticks,
         yticks = yticks, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab,
         xlas = xlas, ylas = ylas, show.grid = show.grid, grid.lwd = grid.lwd, grid.col = grid.col
      )
   }

   # if some rows are excluded remove part of values from cgroup
   if (length(excluded_rows) > 0 && length(cgroup) > 1) {
      cgroup <- cgroup[-excluded_rows]
   }

   #  get proper colors
   if (density == FALSE && !is.null(cgroup) && !(type == "e")) {
      # show color groups according to cdata values
      col <- mdaplot.getColors(cgroup = cgroup, colmap = colmap, opacity = opacity)
   } else {
      # show all points with the same color
      col <- if (is.null(col)) mdaplot.getColors(1, colmap = colmap, opacity = opacity)
         else adjustcolor(col, opacity)

      # set cgroup to NULL so method will not try to show color groups or colorbar legend
      cgroup <- NULL
   }

   # make plot for the data
   if (type == "p" && density == TRUE) type <- "d"
   switch(type,
      "p" = mdaplot.plotScatter(x_values, y_values, x_values_excluded, y_values_excluded,
            pch.colinv = pch.colinv, pch = pch, col = col, bg = bg, col.excluded = col.excluded,
            cex = cex, lwd = lwd, ...),
      "d" = mdaplot.plotDensity(x_values, y_values, nbins = nbins, colmap = colmap),
      "b" = mdaplot.plotLines(x_values, y_values, x_values_excluded, y_values_excluded,
               type = "b", col = col, pch = pch, col.excluded = col.excluded,
               lty = lty, lwd = lwd, cex = cex, ...),
      "l" = mdaplot.plotLines(x_values, y_values, x_values_excluded, y_values_excluded,
               type = "l", col = col, pch = pch, col.excluded = col.excluded,
               lty = lty, lwd = lwd, cex = cex, ...),
      "h" = mdaplot.plotBars(x_values, y_values, col = col, bwd = bwd, border = border,
               force.x.values = force.x.values, ...),
      "e" = mdaplot.plotErrorbars(x_values, y_values, col = col, pch = pch,
               cex = cex, ...)
   )

   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2) {
      mdaplot.showLines(show.lines)
   }

   # show lables
   if (show.labels && !is.null(labels)) {
      mdaplot.showLabels(x_values, y_values, labels, type = type, col = lab.col, cex = lab.cex)
   }

   # show lables for excluded rows
   if (show.labels && show.excluded && length(labels_excluded) > 0) {
      mdaplot.showLabels(
         x_values_excluded, y_values_excluded, labels_excluded,
         type = type, col = lab.col, cex = lab.cex
      )
   }

   # show colorbar if needed
   if (!is.null(cgroup) && show.colorbar) {
      mdaplot.showColorbar(cgroup, colmap, lab.col = lab.col, lab.cex = lab.cex)
   }

   plot.data[["col"]] <- col
   plot.data[["cgroup"]] <- cgroup
   plot.data[["type"]] <- type
   invisible(plot.data)
}

#' Plotting function for several sets of objects
#'
#' @description

#' \code{mdaplotg} is used to make different kinds of plots or their combination for several sets

#' of objects.
#'
#' @param data
#' a matrix, data frame or a list with data values (see details below).
#' @param type
#' type of the plot ('p', 'l', 'b', 'h', 'e').
#' @param pch
#' a character for markers (same as \code{plot} parameter).
#' @param lty
#' the line type (same as \code{plot} parameter).
#' @param lwd
#' the line width (thickness) (same as \code{plot} parameter).
#' @param cex
#' the cex factor for the markers (same as \code{plot} parameter).
#' @param bwd
#' a width of a bar as a percent of a maximum space available for each bar.
#' @param legend
#' a vector with legend elements (if NULL, no legend will be shown).
#' @param xlab
#' a title for the x axis (same as \code{plot} parameter).
#' @param ylab
#' a title for the y axis (same as \code{plot} parameter).
#' @param main
#' an overall title for the plot (same as \code{plot} parameter).
#' @param labels
#' what to use as labels ('names' - row names, 'indices' - row indices, 'values' - values).
#' @param ylim
#' limits for the y axis (if NULL, will be calculated automatically).
#' @param xlim
#' limits for the x axis (if NULL, will be calculated automatically).
#' @param col
#' colors for the plot series
#' @param colmap
#' a colormap to generate colors if \code{col} is not provided
#' @param legend.position
#' position of the legend ('topleft', 'topright', 'top', 'bottomleft', 'bottomright', 'bottom').
#' @param show.legend
#' logical, show or not legend for the data objects.
#' @param show.labels
#' logical, show or not labels for the data objects.
#' @param show.lines
#' vector with two coordinates (x, y) to show horizontal and vertical line cross the point.
#' @param show.grid
#' logical, show or not a grid for the plot.
#' @param grid.lwd
#' line thinckness (width) for the grid
#' @param grid.col
#' line color for the grid
#' @param xticks
#' tick values for x axis.
#' @param xticklabels
#' labels for x ticks.
#' @param yticks
#' tick values for y axis.
#' @param yticklabels
#' labels for y ticks.
#' @param xlas
#' orientation of xticklabels
#' @param ylas
#' orientation of yticklabels
#' @param lab.col
#' color for data point labels.
#' @param lab.cex
#' size for data point labels.
#' @param show.excluded
#' logical, show or hide rows marked as excluded (attribute `exclrows`)
#' @param groupby
#' one or several factors used to create groups of data matrix rows (works if data is a matrix)
#' @param opacity
#' opacity for plot colors (value between 0 and 1)
#' @param ...
#' other plotting arguments.
#'
#' @details
#' The \code{mdaplotg} function is used to make a plot with several sets of objects. Simply
#' speaking, use it when you need a plot with legend. For example to show line plot with spectra
#' from calibration and test set, scatter plot for height and weight values for women and men, and
#' so on.
#'
#' Most of the parameters are similar to \code{\link{mdaplot}}, the difference is described below.
#'
#' The data should be organized as a list, every item is a matrix with data for one set of objects.
#' Alternatively you can provide data as a matrix and use parameter \code{groupby} to create groups.
#' See tutorial for more details.
#'
#' There is no color grouping option, because color is used to separate the sets. Marker symbol,
#' line style and type, etc. can be defined as a single value (one for all sets) and as a vector
#' with one value for each set.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @importFrom graphics abline axis grid hist lines matlines par plot points rect segments text
#' image plot.new rasterImage smoothScatter
#' @importFrom grDevices axisTicks dev.cur
#'
#' @export
mdaplotg <- function(
   data, groupby = NULL, type = "p", pch = 16,  lty = 1, lwd = 1, cex = 1,
   col = NULL, bwd = 0.8, legend = NULL, xlab = NULL, ylab = NULL, main = NULL, labels = NULL,
   ylim = NULL, xlim = NULL, colmap = "default", legend.position = "topright",
   show.legend = TRUE, show.labels = FALSE, show.lines = FALSE, show.grid = TRUE, grid.lwd = 0.5,
   grid.col = "lightgray", xticks = NULL, xticklabels = NULL, yticks = NULL, yticklabels = NULL,
   show.excluded = FALSE, lab.col = "darkgray", lab.cex = 0.65, xlas = 1,
   ylas = 1, opacity = 1, ...) {

   # split data into groups
   name <- attr(data, "name", exact = TRUE)
   data <- mdaplot.prepareDataForGPlots(data, type, groupby)
   ngroups <- length(data)

   # check if plot.new() should be called first
   if (dev.cur() == 1) plot.new()

   # check numeric parameters and multiply them if necessary
   processParam <- function(param, name, is.type, ngroups) {

      param <- if (length(param) == 1) rep(param, ngroups) else param

      if (!all(is.type(param))) {
         stop(paste0('Parameter "', name, '" mush be numeric!'))
      }

      if (length(param) != ngroups)
         stop(paste0('Parameter "', name, '" should be specified for each group or one for all!'))

      return(param)
   }

   type <- processParam(type, "type", is.character, ngroups)
   pch <- processParam(pch, "pch", is.numeric, ngroups)
   lty <- processParam(lty, "lty", is.numeric, ngroups)
   lwd <- processParam(lwd, "lwd", is.numeric, ngroups)
   cex <- processParam(cex, "cex", is.numeric, ngroups)
   opacity <- processParam(opacity, "opacity", is.numeric, ngroups)
   lab.col <- processParam(lab.col, "lab.col", mdaplot.areColors, ngroups)

   # check and define colors if necessary
   if (is.null(col)) col <- mdaplot.getColors(ngroups = ngroups, colmap = colmap)
      col <- processParam(col, "col", mdaplot.areColors, ngroups)

   # get plot data for each group
   pd <- list()
   for (i in 1:ngroups) {
      pd[[i]] <- mdaplot.createPlotData(
         data = data[[i]], type = type[i], main = main, xlab = xlab, ylab = ylab, xlim = xlim,
         ylim = ylim, bwd = bwd, show.excluded = show.excluded, show.colorbar = FALSE,
         show.labels = show.labels, labels = labels, show.lines = show.lines, show.axes = TRUE
      )
   }

   # process legend
   get_legend <- function(legend, data, pd, show.legend) {

      if (!show.legend) {
         return(NULL)
      }

      if (!is.null(legend)) {
         return(legend)
      }

      if (!is.null(names(data))) {
         return(names(data))
      }

      if (all(type == "h")) {
         return(unlist(lapply(pd, function(x) rownames(x$y_values)[1])))
      }

      return(unlist(lapply(pd, function(x) x$data_attrs$name)))
   }

   legend <- get_legend(legend, data, pd, show.legend)

   # check legend if required
   if (show.legend) {

      if (is.null(legend)) {
         stop("Can not find values for the legend items.")
      }

      if (length(legend) != ngroups) {
         stop("Number of values for 'legend' is not the same as number of plot series.")
      }
   }

   # compute x-limits
   if (is.null(xlim)) {
      xlim <- matrix(unlist(lapply(pd, function(x) x$xlim)), ncol = 2, byrow = TRUE)
      xlim <- c(min(xlim[, 1]), max(xlim[, 2]))

      # stretch limits if bar plot should be shown
      xlim <- if (any(type == "h")) xlim + c(-0.35, 0.35) else xlim
   }

   # compute y-limits
   if (is.null(ylim)) {
      ylim <- matrix(unlist(lapply(pd, function(y) y$ylim)), ncol = 2, byrow = TRUE)
      ylim <- c(min(ylim[, 1]), max(ylim[, 2]))
   }

   # add extra margins if legend must be shown
   if (show.legend == T && !is.null(legend)) {

      scale <- 0.1

      # calculate margins: dx and dy
      dx <- diff(xlim) * scale
      dy <- diff(ylim) * scale

      # table for margin factors depending on legend size and position
      legend_margins <- data.frame(
                         #  i  w    dx   j   h    dy
         "bottom"       = c(1, 0,  0.0,  1, -1, -0.2),
         "bottomleft"   = c(1, 0, -0.2,  1,  0, -0.2),
         "bottomright"  = c(2, 0,  0.2,  1,  0, -0.2),
         "top"          = c(1, 0,  0.0,  2,  1,  0.2),
         "topleft"      = c(1, 0, -0.2,  2,  0,  0.2),
         "topright"     = c(2, 0,  0.2,  2,  0,  0.2)
      )

      # get size of legend
      legend_size <- mdaplot.showLegend(legend, position = legend.position, plot = F)
      w <- legend_size$rect$w
      h <- legend_size$rect$h

      # correct limits
      mrg <- matrix(legend_margins[[legend.position]], nrow = 2, byrow = TRUE)
      xlim[mrg[1, 1]] <- xlim[mrg[1, 1]] + w * mrg[1, 2] + dx * mrg[1, 3]
      ylim[mrg[2, 1]] <- ylim[mrg[2, 1]] + h * mrg[2, 2] + dy * mrg[2, 3]
   }


   # define main title if not provided (either as "name" or as "name" attr of first dataset)
   main <- if (is.null(main)) name else main
   main <- if (is.null(main)) pd[[1]]$data_attrs[["name"]] else main

   # define labels for axes
   xlab <- if (is.null(xlab)) attr(pd[[1]]$x_values, "name", exact = TRUE) else xlab
   ylab <- if (is.null(ylab)) attr(pd[[1]]$y_values, "name", exact = TRUE) else ylab

   # make an empty plot with proper limits and axis labels
   mdaplot.plotAxes(xticklabels = xticklabels, yticklabels = yticklabels, xticks = xticks,
      yticks = yticks, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab,
      xlas = xlas, ylas = ylas, show.grid = show.grid,
      grid.lwd = grid.lwd, grid.col = grid.col
   )

   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2) {
      mdaplot.showLines(show.lines)
   }

   # count how many plots are bar plots
   nbarplots <- sum(type == "h")

   # make a plot for each group
   for (i in 1:ngroups) {

      # decide if x values should be forced as group index
      force.x.values <- if (type[i] == "h") c(i, nbarplots) else NA

      # if error bars are shown and i > 1 do not show labels
      show.labels <- if (i > 1 && type[i] == "e") FALSE else show.labels

      # use mdaplot with show.axes = FALSE to create the plot
      mdaplot(plot.data = pd[[i]], type = type[i], col = col[i], pch = pch[i], lty = lty[i],
              lwd = lwd[i], cex = cex[i], force.x.values = force.x.values, bwd = bwd,
              show.grid = F, show.labels = show.labels, opacity = opacity[i],
              lab.col = lab.col[i], lab.cex = lab.cex, show.axes = FALSE, ...
      )
   }

   # show legend if required
   if (show.legend == TRUE) {
      lty[type == "p" | type == "h"] <- 0
      pch[type == "l"] <- NA_integer_
      pch[type == "h"] <- 15

      mdaplot.showLegend(
         legend, col, pch = pch, lty = lty, lwd = lwd, cex = cex,
         position = legend.position
      )
   }
}

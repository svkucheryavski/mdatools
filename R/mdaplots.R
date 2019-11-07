#' Check color values
#'

#' @description
#' Checks if elements of argument are valid color values
#'

#' @param palette
#' vector with possibly color values (names, RGB, etc.)
#'

mdaplot.areColors = function(palette) {
   sapply(palette, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)} )
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
   if (round.only == T)
      fdata <- round(data, digits)
   else
      fdata <- prettyNum(data, digits = digits)

   if (!is.null(dim(data))){
      dim(fdata) <- dim(data)
      dimnames(fdata) <- dimnames(data)
   }

   return(fdata)
}

#' Calculate axes limits
#'

#' @description
#' Calculates axes limits depending on data values that have to be plotted,

#' extra plot elements that have to be shown and margins.

#'

#' @param plot_data
#' a list with plot data and related parameters returned by `create_plot_data()`.
#' @param show.colorbar
#' logical, show or not the colorbar on the plot.
#' @param show.lines
#' logical or numeric with line coordinates to be shown on the plot.
#' @param legend
#' vector with legend items.
#' @param show.legend
#' logical, show or not legend on the plot.
#' @param legend.position
#' position of the legend (see \code{\link{mdaplotg}} for details).
#' @param show.labels
#' logical, show or not labels for the data objects
#' @param show.excluded
#' logical, show or not excluded values
#'
#' @details
#' Data can be a list with several matrices or just one matrix. The matrices can have single.x
#' configuration, where first column is x values and the others are y values or normal

#' configuration, where every odd column is x values and every even is corresponding y values.
#'
#' @return
#' Returns a list with four limits for the x and y axes.
#'
mdaplot.getAxesLim <- function(plot_data, show.colorbar = F, show.lines = F, legend = NULL,
   show.legend = F, legend.position = "topright", show.labels = F, show.excluded = F) {

   if (is.null(plot_data)) return(NULL)

   if (is.null(show.labels)) {
      show.labels <- FALSE
   }

   x_values <- plot_data$x_values
   y_values <- plot_data$y_values

   if (show.excluded) {
      x_values <- c(x_values, plot_data$x_values_excluded)
      y_values <- rbind(y_values, plot_data$y_values_excluded)
   }

   lower <- plot_data$lower
   upper <- plot_data$upper

   # scale for margins - 7.5% of plot width or height
   scale <- 0.075

   # find range for x values
   xmin <- min(x_values)
   xmax <- max(x_values)

   # find range for y values
   ymin <- min(y_values, if (is.null(lower)) y_values else y_values - lower, na.rm = TRUE)
   ymax <- max(y_values, if (is.null(lower)) y_values else y_values + upper, na.rm = TRUE)

   # if min and max are equal add 5% on each side (or just 0.05 if they are zero)
   if (xmin == xmax) {
      xmin <- xmin - 0.05 * ifelse(xmin == 0, 1, xmin)
      xmax <- xmax + 0.05 * ifelse(xmax == 0, 1, xmax)
   }

   # if min and max are equal add 5% on each side (or just 0.05 if they are zero)
   if (ymin == ymax) {
      ymin <- ymin - 0.05 * ifelse(ymin == 0, 1, ymin)
      ymax <- ymax + 0.05 * ifelse(ymax == 0, 1, ymax)
   }

   # correct limits if some lines that have to be shown are outside the data points cloud
   if (is.numeric(show.lines)) {
      xmax <- max(xmax, show.lines[1], na.rm = TRUE)
      xmin <- min(xmin, show.lines[1], na.rm = TRUE)
      ymax <- max(ymax, show.lines[2], na.rm = TRUE)
      ymax <- max(ymax, show.lines[2], na.rm = TRUE)
   }

   # calculate margins: dx and dy
   dx <- (xmax - xmin) * scale
   dy <- (ymax - ymin) * scale

   # define limits with margins
   xlim <- c(xmin - dx, xmax + dx)
   ylim <- c(ymin - dy, ymax + dy)

   # add an extra margin to y limit if colorbar must be shown
   if (show.colorbar == T) ylim[2] <- ylim[2] + dy * 3

   # add extra margins if legend must be shown
   if (show.legend == T && !is.null(legend)) {
      legend.size <- mdaplot.showLegend(legend, position = legend.position, plot = F)
      w <- legend.size$rect$w
      h <- legend.size$rect$h

      if (legend.position == "topleft") {
         xlim[1] <- xlim[1] - dx * 0.2
         ylim[2] <- ylim[2] + dy * 0.2
      } else if (legend.position == "top") {
         ylim[2] = ylim[2] + h + 0.2 * dy
      } else if (legend.position == "topright") {
         xlim[2] <- xlim[2] + dx * 0.2
         ylim[2] <- ylim[2] + dy * 0.2
      } else if (legend.position == "bottomleft") {
         xlim[1] <- xlim[1] - dx * 0.2
         ylim[1] <- ylim[1] - dy * 0.2
      } else if (legend.position == "bottom") {
         ylim[1] = ylim[1] - h - 0.2 * dy
      } else if (legend.position == "bottomright") {
         xlim[2] <- xlim[2] + dx * 0.2
         ylim[1] <- ylim[1] - dy * 0.2
      }
   }

   # add extra margins if labels must be shown
   if (show.labels == T) {
      ylim[1] <- ylim[1] - dy
      ylim[2] <- ylim[2] + dy
   }

   # return the limits
   return(list(xlim = xlim, ylim = ylim))
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
mdaplot.getColors = function(ngroups = NULL, cgroup = NULL, colmap = 'default', opacity = 1, maxsplits = 64) {

   # if non of the main arguments defined assume only one color is needed
   if (is.null(ngroups) && is.null(cgroup)) {
      ngroups <- 1
   }

   # returns palette for given colormap
   get_palette <- function(colmap) {

      if (length(colmap) > 1) {
         return(colmap)
      }

      if (colmap == 'gray') {
         # grayscale
         return(c("#E8E8E8", "#D6D6D6", "#C4C4C4", "#B2B2B2", "#9A9A9A", "#808080", "#484848", "#101010"))
      }

      if (colmap == 'jet') {
         # jet colors
         return(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      }

      if (colmap == 'old') {
         # color brewer colormap used as default in versions before 0.10.0
         return(c("#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"))
      }

      # the current default colormap

      return(c("#2679B2", "#1C9AA8", "#379531", "#EED524", "#FB7F28", "#D22C2F"))
   }

   # returns vector with colors
   prepare_colors <- function(colfunc, ncolors, opacity) {
      colors <- colfunc(ncolors)

      if (is.null(opacity) || all(opacity == 1)) {
         return (colors)
      }

      if (length(opacity) == 1){
         opacity <- rep(opacity, ncolors)
      }

      if (length(opacity) != ncolors) {
         stop('Wrong number of values for "opacity" parameter!')
      }

      for (i in 1:ncolors) {
         colors[i] <- adjustcolor(colors[i], alpha.f = opacity[i])
      }

      return (colors)
   }

   # get palette
   palette <- get_palette(colmap)
   if (!all(mdaplot.areColors(palette))) {
      stop("Parameter 'colmap' must contains valid color values or name of palette!")
   }

   # if grayscale palette and only one color is needed reorder pallete so the black is first
   if (all(colmap == 'gray') && is.null(cgroup) && ngroups == 1) {
      palette <- rev(palette)

   }

   # define color function based on the palette
   colfunc <- colorRampPalette(palette)

   # if cgroup is not provided just return the colors

   if (is.null(cgroup)) {
      return (prepare_colors(colfunc, ngroups, opacity))
   }

   # if cgroup is factor return vector with corresponding values
   if (is.factor(cgroup)) {
      ngroups <- length(attr(cgroup, 'levels'))
      return (prepare_colors(colfunc, ngroups, opacity)[as.numeric(cgroup)])
   }

   # if not split it into groups
   if (is.null(ngroups)) {
      ngroups <- length(unique(cgroup))
      ngroups <- ifelse(ngroups > maxsplits, maxsplits, ngroups)
   }

   cgroup <- cut(as.numeric(cgroup), ngroups, include.lowest = TRUE)
   out.palette <- prepare_colors(colfunc, ngroups, opacity)
   colors <- out.palette[as.numeric(cgroup)]
   attr(colors, 'palette') <- out.palette

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
      labels[i + 1, ] = c(x + w * i, y - h)
   else
      labels[, 1] = labels[, 1] + w/2

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
#' legend position ('topright', 'topleft', 'bottomright', 'bottomleft', 'top', 'bottom')
#' @param plot
#' logical, show legend or just calculate and return its size
#'
mdaplot.showLegend = function(legend, col, pch = NULL, lty = NULL, lwd = NULL, cex = 1,

                              bty = 'o', position = 'topright', plot = T) {
   # which positions need multiple columns
   onecolpos = c('topright', 'topleft', 'bottomright', 'bottomleft')
   multcolpos = c('top', 'bottom', 'right', 'left')

   if (!(position %in% c(onecolpos, multcolpos))) {
      stop("Wrong values for 'legend.position' argument!")
   }

   # compute number of columns
   ncol <- if (position %in% onecolpos) 1 else length(legend)

   # calculate inset values depending on a ration between width and height of a plot
   lim <- par('plt')

   dx <- lim[2] - lim[1]
   dy <- lim[4] - lim[3]

   inset <- c(0.02, 0.02 * (dx/dy))

   # show legend
   legend(position, legend, col = col, pch = pch, lty = lty, pt.cex = cex, lwd = lwd,

          cex = 0.85, plot = plot, inset = inset, bg = 'white', box.lwd = 0.75, box.col = 'gray',
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
   cex = 0.65, col = "darkgray", type = NULL) {


   if (length(labels) == 0 || length(x_values) == 0 || length(y_values) == 0) {
      warning("No labels available.")
      return()
   }

   if (!is.null(dim(y_values)) && nrow(y_values) > 1 && type != "p") {
      y_values <- apply(y_values, 2, max)
   }

   # show labels
   x_values <- as.vector(x_values)
   y_values <- as.vector(y_values)

   if (is.numeric(labels)) {
      labels <- mdaplot.formatValues(labels)
   }

   if (!(any(is.nan(x_values)) || any(is.nan(y_values)))) {
      if (!is.null(type) && type == "h") {
         # show labels properly for bars with positive and negative values
         text(x_values, y_values, labels, cex = cex, pos = ifelse(y_values < 0, 1, 3), col = col)
      } else {
         text(x_values, y_values, labels, cex = cex, pos = pos, col = col)
      }
   }
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
      x <- data[[i]][, 1]
      y <- data[[i]][, 2]
      abline(lm(y ~ x), lty = lty, lwd = lwd, col = col[i])
   }
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
#'
bars <- function(x, y, col = NULL, bwd = 0.8, border = NA) {

   if (length(bwd) == 1)
      bwd <- matrix(bwd, ncol = length(x))

   if (length(col) != length(y))
      col <- rep(col, length.out = length(y))

   rect(x - bwd / 2, 0, x + bwd / 2, y, col = col, border = border)
}

#' Show error bars on a plot
#'
#' @description
#' Shows error bars (errorbar plot) on predefined axes
#'
#' @param x
#' vector with x values
#' @param lower
#' vector with lower limits for the bars
#' @param upper
#' vector with upper limits for the bars
#' @param y
#' vector with y values (bid points)
#' @param col
#' color for the error bars
#' @param pch
#' marker symbol for the plot
#' @param cex
#' cex factor for the marker
#'
errorbars <- function(x, lower, upper, y = NULL, col = NULL, pch = 16, cex = 1) {
   e2 <- (max(x) - min(x)) / 50
   e1 <- (max(x) - min(x)) / (length(x) * 5)
   e <- min(e1, e2)

   segments(x, y - lower, x, y + upper, col = col)
   segments(x - e, y - lower, x + e, y - lower, col = col)
   segments(x - e, y + upper, x + e, y + upper, col = col)

   if (!is.null(y)) {
      points(x, y, col = col, pch = pch, cex = cex)
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
get_confidence_ellipse <- function(points, conf.level = 0.95, n = 100) {

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
get_convex_hull <- function(points) {
   ch_ind <- chull(points)
   ch_ind <- c(ch_ind, ch_ind[1])
   return(points[ch_ind, ])
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
#'
#' @importFrom graphics polygon
#'
#' @export
add_points_shape <- function(plot_data, lwd, lty, conf.level, opacity, shape_function, ...) {
   
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
#' add_confidence_ellipse(p, conf.level = 0.90, opacity = 0.2)
#'
#' @export
add_confidence_ellipse <- function(plot_data, conf.level = 0.95, lwd = 1, lty = 1, opacity = 0) {
   add_points_shape(
      plot_data,
      lwd = lwd,
      lty = lty,
      opacity = opacity,
      shape_function = get_confidence_ellipse,
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
#' add_convex_hull(p, opacity = 0.2)
#'
#' @importFrom grDevices chull
#' @importFrom graphics polygon
#'
#' @export
add_convex_hull <- function(plot_data, lwd = 1, lty = 1, opacity = 0) {
   add_points_shape(
      plot_data, 
      lwd = lwd, 
      lty = lty, 
      opacity = opacity, 
      shape_function = get_convex_hull
   )
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
#' @param lim
#' vector with limits for x and y axis
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
mdaplot.plotAxes <- function(xticklabels = NULL, yticklabels = NULL, xticks = NULL, yticks = NULL,
                            lim = NULL, main = NULL, xlab = NULL, ylab = NULL, xlas = 0, ylas = 0,
                            show.grid = TRUE, grid.lwd = 0.5, grid.col = "lightgray") {

   # check xticklabels and xticks
   if (!is.null(xticklabels)) {

      if (is.null(xticks)) {
         stop('Please, specify "xticks" vector for corresponding "xticklabels".')
      }

      if (length(xticks) != length(xticklabels)) {
         stop('Number of elements in "xticks" and "xticklabels" should be the same')
      }
   }

   # check yticklabels and yticks
   if (!is.null(yticklabels)) {

      if (is.null(yticks)) {
         stop('Please, specify "yticks" vector for corresponding "yticklabels".')
      }

      if (length(yticklabels) != length(yticks)) {
         stop('Number of elements in "yticks" and "yticklabels" should be the same')
      }

   }

   # make plot without ticks
   plot(0, 0, type = "n", main = main, xlab = xlab, ylab = ylab, xlim = lim$xlim, ylim = lim$ylim,
        xaxt = "n", yaxt = "n")

   # generate x ticks
   if (is.null(xticks)) {
      xticks <- axisTicks(lim$xlim, log = FALSE)
   }

   # show x-axis
   if (is.null(xticklabels)) {
      axis(1, at = xticks, las = xlas)
   } else {
      axis(1, at = xticks, labels = xticklabels, las = xlas)
   }

   # generate y ticks
   if (is.null(yticks)) {
      yticks = axisTicks(lim$ylim, log = FALSE)
   }

   # show y-axis
   if (is.null(yticklabels)) {
      axis(2, at = yticks, las = ylas)
   } else {
      axis(2, at = yticks, labels = yticklabels, las = ylas)
   }

   # show grid if needed
   if (show.grid) {
      grid(lwd = grid.lwd, col = grid.col)
   }
}

#' Check if data argument looks correct for creatins plots, correct if needed
#'

#' @param data
#' dataset provided to mdaplot function
#' @param type
#' type of plot
#'
#' @export
prepare_data_for_plot <- function(data, type) {

   # save data arguments
   data_attrs <- attributes(data)

   # make sure the data set is provided
   if (is.null(data) || length(data) < 1) {
      stop("The provided dataset is empty!")
   }

   # if it is vector without dimension - make a matrix

   if (is.null(dim(data))) {
      data <- if (type == "p") as.matrix(data) else t(as.matrix(data))
   }

   # convert data frame to a matrix if needed
   if (is.data.frame(data)) {
      data <- mda.df2mat(data)
   }

   # if row names are missing - add them
   if (is.null(rownames(data))) {
      rownames(data) <- paste0("O", 1:nrow(data))
   }

   # if column names are missing - add them
   if (is.null(colnames(data))) {
      colnames(data) <- paste0("X", 1:ncol(data))
   }

   # handle excluded columns if any

   excluded_cols <- data_attrs$exclcols
   col_indices <- 1:ncol(data)
   if (length(excluded_cols) > 0) {
      excluded_cols <- mda.getexclind(excluded_cols, colnames(data), ncol(data))

      data <- data[, -excluded_cols, drop = F]
      col_indices <- col_indices[-excluded_cols]
   }

   # check if data still has some columns
   if (is.null(dim(data)) || ncol(data) < 1) {
      stop("No columns left when excluded hidden values.")
   }

   # add columns with row numbers or yaxis values to the data if necessary
   if (ncol(data) == 1 && type == "p") {
      rownames = rownames(data)
      new_column <- if(is.null(data_attrs$yaxis.values)) 1:nrow(data) else data_attrs$yaxis.values
      data <- mda.cbind(new_column, data)
      colnames(data)[1] <- if(is.null(data_attrs$yaxis.name)) 'Objects' else data_attrs$yaxis.name
      rownames(data) <- rownames
   }

   # handle excluded rows if necessary
   row_indices <- 1:nrow(data)
   excluded_data <- NULL
   excluded_rows <- data_attrs$exclrows
   if (length(excluded_rows) > 0) {
      excluded_rows <- mda.getexclind(excluded_rows, rownames(data), nrow(data))

      excluded_data <- data[excluded_rows, , drop = F]
      data <- data[-excluded_rows, , drop = F]
      row_indices <- row_indices[-excluded_rows]
   }

   # check that number of rows is still sufficient
   if (is.null(dim(data)) || nrow(data) < 1) {
      stop("No rows left when excluded hidden values.")
   }

   # at least two rows (y values and size of error bars) are needed for errorbar plot
   if (type == "e" && nrow(data) < 2) {
      stop("Errorbar plot requires dataset with at least two rows!")
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

#' Split dataset to x and y values depending on plot type
#'

#' @param plot_data
#' list with data and parameters prepared for plot by `prepare_data_for_plot()` function
#' @param type
#' type of plot
#'

#' @export
split_plot_data <- function(plot_data, type) {

   # shortcuts to some of parameters
   data <- plot_data$data
   data_attrs <- plot_data$attrs
   excluded_cols <- plot_data$excluded_cols
   xaxis_name <- data_attrs$xaxis.name
   yaxis_name <- data_attrs$yaxis.name

   if (type == "p") {
      # scatter plot
      plot_data$x_values <- data[, 1]
      attr(plot_data$x_values, "name") <- colnames(data)[1]
      plot_data$y_values <- data[, 2, drop = F]
      attr(plot_data$y_values, "name") <- colnames(data)[2]
      return(plot_data)
   }

   # prepare x-axis values for other types of plots
   x_values <- 1:ncol(data)
   if (!is.null(data_attrs$xaxis.values)) {
      x_values <- data_attrs$xaxis.values
      if (length(excluded_cols) > 0) x_values[-excluded_cols]
   }

   if (is.null(names(x_values))) names(x_values) <- colnames(data)

   # add name to x-values (and names)
   attr(x_values, "name") <- if (is.null(xaxis_name)) "Variables" else xaxis_name
   plot_data$x_values <- x_values

   # for bar plot we show only the first row
   if (type == "h") {
      plot_data$y_values <- data[1, , drop = F]
      attr(plot_data$y_values, "name") <- if (is.null(yaxis_name)) "" else yaxis_name
      return(plot_data)
   }

   # for line plot we use all values from data
   if (type == "l" || type == "b") {
      plot_data$y_values <- data
      attr(plot_data$y_values, "name") <- if (is.null(yaxis_name)) "" else yaxis_name

      return(plot_data)
   }

   # for errobar plot first row is y-values the other one or two rows are limits for error bars
   if (type == "e") {
      plot_data$y_values <- data[1, , drop = F]
      attr(plot_data$y_values, "name") <- rownames(data)[1]
      plot_data$lower <- data[2, ]
      plot_data$upper <- if (nrow(data) > 2) data[3, ] else plot_data$lower

      return(plot_data)
   }

   stop("Something went wrong on plot data preparation step (check plot type).")
}

#' Prepare x and y values for excluded rows

#'

#' @param plot_data
#' list with data and parameters prepared for plot by `prepare_data_for_plot()` function
#' @param type
#' type of plot
#'

#' @export
process_excluded_rows <- function(plot_data, type) {

   # shortcuts for some parameters
   excluded_data <- plot_data$excluded_data

   # if nothing to show - return
   if (is.null(excluded_data) || (nrow(excluded_data) == 0)) {
      return(plot_data)
   }

   if (type == "e" || type == "h") {
      stop("Bar plor and errorbar plot can not be made for matrix with excluded rows.")
   }

   if (type == "p") {
      plot_data$x_values_excluded <- excluded_data[, 1]
      plot_data$y_values_excluded <- excluded_data[, 2, drop = F]

      return(plot_data)
   }

   if (type == "l" || type == "b") {
      plot_data$x_values_excluded <- plot_data$x_values
      plot_data$y_values_excluded <- excluded_data

      return(plot_data)
   }
}

#' Create a vector with labels for plot data series
#'

#' @param plot_data
#' list with data and parameters prepared for plot by `prepare_data_for_plot()` function
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
prepare_plot_data_labels <- function(plot_data, type, labels) {

   # shortcuts for some parameters
   x_values <- plot_data$x_values
   y_values <- plot_data$y_values
   x_values_excluded <- plot_data$x_values_excluded
   y_values_excluded <- plot_data$y_values_excluded
   excluded_rows <- plot_data$excluded_rows
   excluded_cols <- plot_data$excluded_cols

   # if user provided labels - use them
   if (!is.null(labels) && length(labels) > 1) {

      if (type == "p") {
         # scatter plots

         if (length(labels) != (nrow(y_values) + length(excluded_rows))) {
            stop("Labels must be provided for all data rows (including excluded)")
         }

         if (is.null(excluded_rows)) {
            plot_data$labels_excluded <- NULL
            plot_data$labels <- labels

            return(plot_data)
         }

         plot_data$labels_excluded <- labels[excluded_rows]
         plot_data$labels <- labels[-excluded_rows]

         return(plot_data)
      }

      if (type != "p") {
         # non scatter plots

         if (length(labels) != (length(x_values) + length(excluded_cols))) {
            stop("Labels must be provided for all data columns (including excluded)")
         }

         if (is.null(excluded_cols)) {
            plot_data$labels_excluded <- NULL
            plot_data$labels <- labels

            return(plot_data)
         }

         plot_data$labels_excluded <- labels[excluded_cols]
         plot_data$labels <- labels[-excluded_cols]

         return(plot_data)
      }
   }

   # if labels were not provided, by default use names
   if (is.null(labels) || labels == "names") {
      plot_data$labels_excluded <- if (type == "p") names(x_values_excluded) else NULL
      plot_data$labels <- names(x_values)

      return(plot_data)
   }

   # if labels must be values use y-values for that

   if (labels == "values") {
      labels <- y_values
      labels_excluded <- y_values_excluded

      if (type == "l" || type == "b") {
         labels <- apply(labels, 2, max)
         labels_excluded <- if (!is.null(labels_excluded)) apply(labels_excluded, 2, max)
      }

      plot_data$labels <- labels
      plot_data$labels_excluded <- labels_excluded
      return(plot_data)
   }

   # if labels must be indices use row or column indices

   if (labels == "indices") {

      if (type == "p") {
         plot_data$labels <- plot_data$row_indices
         plot_data$labels_excluded <- excluded_rows
         return(plot_data)
      }

      plot_data$labels <- plot_data$col_indices
      plot_data$labels_excluded <- NULL
      return(plot_data)
   }
}

#' Create dataset to show as plot
#'

#' @param data
#' dataset (vector, matrix or data frame)
#' @param type
#' type of the plot
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
create_plot_data <- function(data, type, xlim, ylim, bwd, show.excluded,
   show.colorbar, show.labels, show.lines, show.axes, labels) {

   # check plot type

   valid.types <- c("p", "l", "b", "h", "e", "d")
   if (!(type %in% valid.types)) {
      stop("Wrong plot type!")
   }

   # get attributes and prepare dataset
   data_attrs <- attributes(data)
   plot_data <- prepare_data_for_plot(data, type)

   # check if density plot should be shown

   plot_data$density <- FALSE
   if (type == "d") {
      if (ncol(data) < 2) stop("Dataset with at least two columns required for density plot.")
      type <- "p"
      plot_data$density <- TRUE
   }

   # split plot data to x and y values
   plot_data <- split_plot_data(plot_data, type)

   # process excluded rows if they have to be shown
   if (show.excluded) {
      plot_data <- process_excluded_rows(plot_data, type)
   }

   # create labels if necessary
   if (show.labels) {
      plot_data <- prepare_plot_data_labels(plot_data, type, labels)
   }

   # define name for the plot series
   plot_data$name <- if (type == "h") rownames(plot_data$data)[1] else data_attrs[["name"]]

   # compute limits for axes
   plot_data$lim <- mdaplot.getAxesLim(
      plot_data,
      show.colorbar = show.colorbar,
      show.labels = show.labels,
      show.lines = show.lines,
      show.excluded = show.excluded
   )

   return(plot_data)
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
#' @param force.x_values
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
      pch = 16, col = NULL, bg = par("bg"), 
      lty = 1, lwd = 1, cex = 1, bwd = 0.8,
      cgroup = NULL, xlim = NULL, ylim = NULL, colmap = "default", labels = NULL,
      main = NULL, xlab = NULL, ylab = NULL, show.labels = F,
      show.colorbar = TRUE, show.lines = FALSE, show.grid = TRUE, grid.lwd = 0.5,
      grid.col = "lightgray", show.axes = TRUE, xticks = NULL, yticks = NULL,
      xticklabels = NULL, yticklabels = NULL,
      xlas = 0, ylas = 0, lab.col = "darkgray", lab.cex = 0.65,
      show.excluded = FALSE, col.excluded = "#E0E0E0", nbins = 256,
      colramp = mdaplot.getColors, force.x_values = NA, opacity = 1,
      pch.colinv = FALSE, ...) {

   if (is.null(plot.data)) {
      # get the data for plot
      plot.data <- create_plot_data(
         data, type, xlim, ylim, bwd, show.excluded, show.colorbar,

         show.labels, show.lines, show.axes, labels
      )
   }

   # get values from the plot data list to make code cleaner
   plot.data$lower -> lower
   plot.data$upper -> upper
   plot.data$x_values -> x_values
   plot.data$y_values -> y_values
   plot.data$x_values_excluded -> x_values_excluded
   plot.data$y_values_excluded -> y_values_excluded
   plot.data$excluded_rows -> excluded_rows
   plot.data$excluded_cols -> excluded_cols
   plot.data$lim -> lim
   plot.data$attrs -> data_attrs
   plot.data$density -> density
   plot.data$labels -> labels
   plot.data$labels_excluded -> labels_excluded

   # if some rows are excluded remove part of values from cgroup

   if (length(excluded_rows) > 0 && !is.null(cgroup) && length(cgroup) > 1) {
      cgroup <- cgroup[-excluded_rows]
   }

   # processs labels and ticklabels
   if (show.axes == T) {
      # if some columns are excluded and xticklabels provided for all - exclude hidden values
      if (!is.null(excluded_cols) && !is.null(xticklabels) && length(xticklabels) == ncol(data)) {
         xticklabels <- xticklabels[-excluded_cols]
      }

      # number of x-values

      nxvals <- length(x_values)

      # define title and axis labels
      if (is.null(main)) main <- plot.data$name
      if (is.null(xlab)) xlab <- attr(x_values, "name")
      if (is.null(ylab)) ylab <- attr(y_values, "name")

      # correct xticklabels and xticks if only one variable is provided for line or bar plot
      if (type != "p" && nxvals == 1) {
         xticks <- 1
         if (is.null(xticklabels)) xticklabels <- x_values
      }

      # correct x axis limits and bar width if it is a bar plot
      if (type == "h") {
         if (nxvals == 1) {
            lim$xlim <- c(x_values - bwd, x_values + bwd)
         } else {
            bwd <- bwd * min(diff(x_values))
            lim$xlim <- c(min(x_values) - bwd / 2, max(x_values) + bwd / 2)
         }
      }

      # redefine x axis limits if user specified values
      if (!is.null(xlim)) {
         lim$xlim <- xlim
      }

      # redefine y axis limits if user specified values
      if (!is.null(ylim)) {
         lim$ylim <- ylim
      }

      # if no x-values were provided to line or bar plot we generate them as integer numbers
      if ( (type %in% c("h", "l", "b")) && is.null(xticks) && is.null(data_attrs$xaxis.values)) {
         xticks <- axisTicks(lim$xlim, log = FALSE)
         xticks <- unique(round(xticks[xticks > 0 & xticks <= max(x_values)]))
      }

      # make an empty plot with proper limits and axis labels
      mdaplot.plotAxes(xticklabels, yticklabels, xticks, yticks, lim, main, xlab, ylab, xlas, ylas,
         show.grid = show.grid, grid.lwd = grid.lwd, grid.col = grid.col)
   }

   #  get proper colors
   if (density == FALSE && !is.null(cgroup) && !(type == "e")) {
      # show color groups according to cdata values
      col <- mdaplot.getColors(cgroup = cgroup, colmap = colmap, opacity = opacity)
   } else {
      # show all points with the same color
      if (!is.null(col)) {
         col <- adjustcolor(col, opacity)
      }  else {
         col <- mdaplot.getColors(1, colmap = colmap, opacity = opacity)
      }

      # set cgroup to NULL so method will not try to show color groups or colorbar legend
      cgroup <- NULL
   }

   # correct x_values if they were forced by bwd
   if (is.numeric(force.x_values)) {
      x_values <- x_values - bwd / 2 + (force.x_values[1] - 0.5) * bwd / force.x_values[2]
      bwd <- bwd / force.x_values[2]
   }

   # make plot for the data

   if (type == "p" && density == FALSE) {
      if (pch.colinv) {
         bg.old <- bg
         bg <- col
         col <- bg.old
      }
      points(x_values, y_values, col = col, bg = bg, pch = pch, lwd = lwd, cex = cex, ...)
   } else if (type == "d" && density == TRUE) {
      par(new = T)
      smoothScatter(
         x_values, y_values, nbin = nbins, cex = cex,
         colramp = colramp, axes = F, ann = F,
         ...
      )
   } else if (type == "l" || type == "b") {
      matlines(
         x_values, t(y_values), type = type, col = col,
         pch = pch, cex = cex, lty = lty, lwd = lwd,
         ...
      )
   } else if (type == "h") {
      bars(x_values, y_values, col = col, bwd = bwd)
   } else if (type == "e") {
      errorbars(x_values, lower, upper, y = y_values, col = col, cex = cex, pch = pch)
   }

   # show excluded rows
   if (!is.null(x_values_excluded) && !is.null(y_values_excluded)) {
      if (type == "p") {
         points(x_values_excluded, y_values_excluded, type = type, col = col.excluded,
            pch = pch, lwd = lwd, cex = cex, ...)
      } else if (type == "l" || type == "b") {
         matlines(x_values_excluded, t(y_values_excluded), type = type, col = col.excluded,
            pch = pch, lty = lty, lwd = lwd, cex = cex, ...)
      }
   }

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
   show.legend = T, show.labels = F, show.lines = F, show.grid = TRUE, grid.lwd = 0.5,
   grid.col = "lightgray", xticks = NULL, xticklabels = NULL, yticks = NULL, yticklabels = NULL,
   show.excluded = FALSE, lab.col = "darkgray", lab.cex = 0.65, xlas = 1,
   ylas = 1, opacity = 1, ...) {

   # name for the plot (mdain)
   name <- NULL

   # prepare list with groups of objects
   if (is.matrix(data) || is.data.frame(data)) {
      if (is.null(groupby)) {
         # take every line as a group
         name <- attr(data, "name", exact = TRUE)

         if (!all(type %in% c("h", "l", "b"))) {
            stop("Group plot with matrix or dataframe can be made only for types 'h', 'l' and 'b'.")
         }

         # split data into a list of subsets for each group
         data.list <- list()
         for (i in seq_len(nrow(data))) {
            data.list[[i]] <- mda.subset(data, subset = i)
         }

         # get rownames as legend values
         if (is.null(legend)) {
            legend <- rownames(data)
         }

         # redefine the data with list
         data <- data.list
      } else {
         # split rows into groups using data frame with factors
         # check that groupby is a factor or data frame with factor columns
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
         data <- split(data, groupby)
         for (i in seq_len(length(data))) {
            row.ind <- data[[i]]$row.ind
            data[[i]] <- subset(data[[i]], select = -row.ind)
            data[[i]] <- mda.setattr(data[[i]], attrs)
            attr(data[[i]], "exclrows") <- which(row.ind %in% attrs$exclrows)
            attr(data[[i]], "labels") <- row.ind
         }
      }
   }

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
      pd[[i]] <- create_plot_data(
         data[[i]], type[i], xlim, ylim, bwd, show.excluded,
         show.colorbar = FALSE, show.labels, labels = labels,
         show.lines, show.axes = TRUE # ! why the hell it is TRUE?
      )
   }

   # process legend
   getLegend <- function(legend, data, pd, show.legend) {

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
         return(unlist(lapply(pd, function(x) {rownames(x$y_values)[1]})))
      }

      return(unlist(lapply(pd, function(x) {x$data_attrs$name})))
   }

   legend <- getLegend(legend, data, pd, show.legend)

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
      xlim <- matrix(unlist(lapply(pd, function(x) {x$lim$xlim})), ncol = 2, byrow = TRUE)
      xlim <- c(min(xlim[, 1]), max(xlim[, 2]))
      # stretch limits if bar plot should be shown
      xlim <- if (any(type == "h")) xlim + c(-0.35, 0.35) else xlim
   }

   # compute y-limits
   if (is.null(ylim)) {
      ylim <- matrix(unlist(lapply(pd, function(y) {y$lim$ylim})), ncol = 2, byrow = TRUE)
      ylim <- c(min(ylim[, 1]), max(ylim[, 2]))
   }

   # combine limits to list
   lim <- list(xlim = xlim, ylim = ylim)

   # define main title if not provided (either as "name" or as "name" attr of first dataset)
   main <- if (is.null(main)) name else main
   main <- if (is.null(main)) pd[[1]]$data_attrs[["name"]] else main

   # define labels for axes
   xlab <- if (is.null(xlab)) attr(pd[[1]]$x_values, "name", exact = TRUE) else xlab
   ylab <- if (is.null(ylab)) attr(pd[[1]]$y_values, "name", exact = TRUE) else ylab

   # make an empty plot with proper limits and axis labels
   mdaplot.plotAxes(xticklabels, yticklabels, xticks, yticks, lim, main, xlab, ylab, xlas, ylas,
      show.grid = show.grid, grid.lwd = grid.lwd, grid.col = grid.col)

   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2) {
      mdaplot.showLines(show.lines)
   }

   # count how many plots are bar plots
   nbarplots <- sum(type == "h")

   # make a plot for each group
   for (i in 1:ngroups) {

      # decide if x values should be forced as group index
      force.x_values <- if (type[i] == 'h') c(i, nbarplots) else NA

      # if error bars are shown and i > 1 do not show labels
      show.labels <- ifelse(i > 1 && type[i] == "e", FALSE, show.labels)

      # use mdaplot with show.axes = FALSE to create the plot
      mdaplot(plot.data = pd[[i]], type = type[i], col = col[i], pch = pch[i], lty = lty[i],
              lwd = lwd[i], cex = cex[i], force.x_values = force.x_values, bwd = bwd,
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

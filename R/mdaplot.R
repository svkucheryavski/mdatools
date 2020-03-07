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
   for (i in seq_len(ncol)) {
      rect(shift + x + w * (i - 1), y, x + w * i, y - h, col = col[i], border = NA)
      labels[i, ] <- c(x + w * (i - 1), y - h)
   }

   # add last value or shift coordinates if labels shall be centered
   if (nrow(labels) > i) {
      labels[i + 1, ] <- c(x + w * i, y - h)
   } else {
      labels[, 1] <- labels[, 1] + w / 2
   }

   # show labels for colorbar regions
   text(labels[, 1], labels[, 2], labels = rownames(labels), pos = 1, col = lab.col,
      cex = lab.cex)
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


#' Color values for plot elements
#'
#' @description
#' Generate vector with color values for plot objects (lines, points, bars), depending
#' on number of groups for the objects.
#'
#' @param ngroups
#' number of colors to create.
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

   # list with currently supported color maps
   colmaps <- list(
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
      "default" = c(
         "#2679B2", "#1C9AA8", "#379531",
         "#EED524", "#FB7F28", "#D22C2F"
      )
   )

   # define palette (if colmap has more than one value - take it as palette)
   palette <- if (length(colmap) > 1) colmap else colmaps[[colmap]]

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

   # if number of groups is larger make binning of cgroup parameter
   if (ngroups > 1) {
      cgroup <- cut(as.numeric(cgroup), ngroups, include.lowest = TRUE)
   }

   # create colors and return
   out_palette <- mdaplot.prepareColors(palette, ngroups, opacity)
   colors <- out_palette[as.numeric(cgroup)]
   attr(colors, "palette") <- out_palette

   return(colors)
}


#' Calculate limits for x-axis.
#'
#' @description
#' Calculates limits for x-axis depending on data values that have to be plotted,
#' extra plot elements that have to be shown and margins.
#'
#' @param ps
#' `plotseries` object.
#' @param xlim
#' limits provided by user
#' @param show.labels
#' logical, will data labels be shown on the plot
#' @param show.lines
#' logical or numeric with line coordinates to be shown on the plot.
#' @param show.excluded
#' logical, will excluded values be shown on the plot
#' @param bwd
#' if limits are computed for bar plot, this is a bar width (otherwise NULL)
#'
#' @return
#' Returns a vector with two limits.
#'
mdaplot.getXAxisLim <- function(ps, xlim, show.labels = FALSE, show.lines = FALSE,
   show.excluded = FALSE, bwd = 0.8) {

   # if user provided limits for x - use them
   if (!is.null(xlim)) return(xlim)

   # x axis limits in case of bar plot
   if (ps$type == "h") {
      values <- ps$x_values
      bwd <- if (length(values) == 1) 2 * bwd else bwd * min(diff(values))
      return(c(min(values) - bwd / 2, max(values) + bwd / 2))
   }

   # if excluded values must be shown - correct internal limits
   xlim <- ps$xlim

   # correct if limits are equal
   if (diff(xlim) == 0) xlim <- xlim * c(0.95, 1.05)

   if (show.excluded && !is.null(ps$x_values_excluded)) {
      xlim_excluded <- range(ps$x_values_excluded)
      xlim <- c(
         min(xlim[1], xlim_excluded[1], na.rm = TRUE),
         max(xlim[2], xlim_excluded[2], na.rm = TRUE)
      )
   }

   # if labels must be shown increase the upper limit
   if (show.labels) {
      xlim[1] <- xlim[1] - diff(xlim) * 0.05
      xlim[2] <- xlim[2] + diff(xlim) * 0.05
   }

   # find if show.lines is in use
   if (length(show.lines) == 2 && is.numeric(show.lines[1])) {
      xlim <- c(
         min(xlim[1], show.lines[1], na.rm = TRUE),
         max(xlim[2], show.lines[1], na.rm = TRUE)
      )
   }

   # add extra margins (3.5%)
   m <- diff(xlim) * 0.035
   xlim <- xlim + c(-m, m)

   return(xlim)
}

#' Calculate limits for y-axis.
#'
#' @description
#' Calculates limits for y-axis depending on data values that have to be plotted,
#' extra plot elements that have to be shown and margins.
#'
#' @param ps
#' `plotseries` object.
#' @param ylim
#' limits provided by user
#' @param show.lines
#' logical or numeric with line coordinates to be shown on the plot.
#' @param show.excluded
#' logical, will excluded values be shown on the plot
#' @param show.labels
#' logical, will data labels be shown on the plot
#' @param show.colorbar
#' logical, will colorbar be shown on the plot
#'
#' @return
#' Returns a vector with two limits.
#'
mdaplot.getYAxisLim <- function(ps, ylim, show.lines = FALSE, show.excluded = FALSE,
   show.labels = FALSE, show.colorbar = FALSE) {

   # if user provided limits for y - use them
   if (!is.null(ylim)) return(ylim)

   # get computed data limits
   ylim <- ps$ylim

   # correct if limits are equal
   if (diff(ylim) == 0) ylim <- ylim * c(0.95, 1.05)

   # if excluded values must be shown - correct internal limits
   if (show.excluded && !is.null(ps$y_values_excluded)) {
      ylim_excluded <- range(ps$y_values_excluded)
      ylim <- c(
         min(ylim[1], ylim_excluded[1], na.rm = TRUE),
         max(ylim[2], ylim_excluded[2], na.rm = TRUE)
      )
   }

   # if labels must be shown increase the upper limit
   if (show.labels) {
      ylim[2] <- ylim[2] + diff(ylim) * 0.075
   }

   # if it is a bar plot and some bars look "down" correct the bottom limit as well
   if (show.labels && ps$type == "h" && any(ps$y_values < 0)) {
      ylim[1] <- ylim[1] - diff(ylim) * 0.05
   }

   # find if show.lines is in use
   if (length(show.lines) == 2 && is.numeric(show.lines[2])) {
      ylim <- c(
         min(ylim[1], show.lines[2], na.rm = TRUE),
         max(ylim[2], show.lines[2], na.rm = TRUE)
      )
   }

   # add an extra margin to y limit if colorbar must be shown
   if (show.colorbar) {
      ylim[2] <- ylim[2] + diff(ylim) * 0.20
   }

   # add extra margins (3.5%)
   m <- diff(ylim) * 0.035
   ylim <- ylim + c(-m, m)
   return(ylim)
}

#' Prepare xticks for plot
#'
#' @param xticks
#' xticks provided by user (if any)
#' @param xlim
#' limits for x axis
#' @param x_values
#' x values for the plot data object
#' @param type
#' type of the plot
#'
#' @export
mdaplot.getXTicks <- function(xticks, xlim, x_values = NULL, type = NULL) {

   if (!is.null(xticks)) return(xticks)
   if (type != "p" && length(x_values) == 1) return(1)
   return(axisTicks(xlim, log = FALSE))
}

#' Prepare yticks for plot
#'
#' @param yticks
#' yticks provided by user (if any)
#' @param ylim
#' limits for y axis
#' @param y_values
#' y values for the plot data object
#' @param type
#' type of the plot
#'
#' @export
mdaplot.getYTicks <- function(yticks, ylim, y_values = NULL, type = NULL) {

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
#' @export
mdaplot.getXTickLabels <- function(xticklabels, xticks, excluded_cols) {

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
#' @export
mdaplot.getYTickLabels <- function(yticklabels, yticks, excluded_rows) {

   if (is.null(yticklabels)) return(TRUE)
   if (is.null(yticks)) stop("You need to specify both 'yticklabels' and 'yticks'")

   # if yticklabels were provided - remove excluded rows if any and check the length
   if (length(yticks) != length(yticklabels)) {
      stop('Number of elements in "yticks" and "yticklabels" should be the same')
   }

   return(yticklabels)
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

#' Plotting function for a single set of objects
#'
#' @description
#' \code{mdaplot} is used to make different kinds of plot for one set of data objects.
#'
#' @param data
#' a vector, matrix or a data.frame with data values.
#' @param ps
#' `plotseries` object, if NULL will be created based on the provided data values
#' @param type
#' type of the plot ("p", "d", "l", "b", "h", "e").
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
#' @param bwd
#' a width of a bar as a percent of a maximum space available for each bar.
#' @param border
#' color for border of bars (if barplot is used)
#' @param lty
#' line type
#' @param lwd
#' line width
#' @param cex
#' scale factor for the marker
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
#' a vector with text labels for data points or one of the following: "names", "indices", "values".
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
#' scheme (\code{colmap = "gray"}) and user can define its own just by specifing the needed
#' sequence of colors (e.g. \code{colmap = c("red", "yellow", "green")}, two colors is minimum).
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
#' @importFrom grDevices axisTicks dev.cur
#' @importFrom graphics abline axis grid hist image lines matlines par
#' @importFrom graphics plot plot.new points rasterImage rect segments text mtext
#'
#' @export
mdaplot <- function(data = NULL, ps = NULL, type = "p",
   pch = 16, col = NULL, bg = par("bg"), bwd = 0.8, border = NA, lty = 1, lwd = 1, cex = 1,
   cgroup = NULL, xlim = NULL, ylim = NULL, colmap = "default", labels = NULL,
   main = NULL, xlab = NULL, ylab = NULL, show.labels = FALSE,
   show.colorbar = !is.null(cgroup), show.lines = FALSE, show.grid = TRUE, grid.lwd = 0.5,
   grid.col = "lightgray", show.axes = TRUE, xticks = NULL, yticks = NULL,
   xticklabels = NULL, yticklabels = NULL,
   xlas = 0, ylas = 0, lab.col = "darkgray", lab.cex = 0.65,
   show.excluded = FALSE, col.excluded = "#C0C0C0", nbins = 60,
   force.x.values = NA, opacity = 1,
   pch.colinv = FALSE, ...) {

   if (is.null(ps)) {
      # get the data for plot
      ps <- plotseries(data, type, col = col, cgroup = cgroup, colmap = colmap, labels = labels)
   }

   # if cgroup is NULL for plot data - color grouping is not allowed
   if (is.null(ps$cgroup)) {
      show.colorbar <- FALSE
   }

   # show axes if needed
   if (show.axes) {
      # get limits for axes
      xlim <- mdaplot.getXAxisLim(ps, xlim = xlim, show.labels = show.labels,
         show.excluded = show.excluded, show.lines = show.lines,
         bwd = (if (type == "h") bwd else NULL))

      ylim <- mdaplot.getYAxisLim(ps, ylim = ylim, show.excluded = show.excluded,
         show.lines = show.lines, show.labels = show.labels, show.colorbar = show.colorbar)

      # check and prepare xticklabels
      xticklabels <- mdaplot.getXTickLabels(xticklabels, xticks, ps$excluded_cols)
      xticks <- mdaplot.getXTicks(xticks, xlim, ps$x_values, type)

      # check and prepare yticklabels
      yticklabels <- mdaplot.getYTickLabels(yticklabels, yticks, ps$excluded_rows)
      yticks <- mdaplot.getYTicks(yticks, ylim, ps$y_values, type)

      # define title and labels
      if (is.null(main)) main <- ps$name
      if (is.null(xlab)) xlab <- attr(ps$x_values, "name")
      if (is.null(ylab)) ylab <- attr(ps$y_values, "name")

      # make an empty plot with proper limits and axis labels
      mdaplot.plotAxes(xticklabels = xticklabels, yticklabels = yticklabels, xticks = xticks,
         yticks = yticks, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab,
         xlas = xlas, ylas = ylas, show.grid = show.grid, grid.lwd = grid.lwd, grid.col = grid.col
      )
   }

   # make plot for the data
   switch(type,
      "p" = plotScatter(ps, pch.colinv = pch.colinv, pch = pch, bg = bg, lwd = lwd, cex = cex,
         col.excluded = col.excluded, show.excluded = show.excluded, ...),
      "d" = plotDensity(ps, nbins = nbins, colmap = colmap),
      "l" = plotLines(ps, pch = pch, lwd = lwd, lty = lty, cex = cex, show.excluded = show.excluded,
         col.excluded = col.excluded, ...),
      "b" = plotLines(ps, pch = pch, lwd = lwd, cex = cex, show.excluded = show.excluded,
         col.excluded = col.excluded, ...),
      "h" = plotBars(ps, bwd = bwd, border = border, force.x.values = force.x.values, ...),
      "e" = plotErrorbars(ps, pch = pch, lwd = lwd, cex = cex, ...),
      stop("Wrong plot type.")
   )

   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2) {
      mdaplot.showLines(show.lines)
   }

   # show lables
   if (show.labels) {
      showLabels(ps, show.excluded = show.excluded, col = lab.col, cex = lab.cex,
         force.x.values = force.x.values, bwd = bwd)
   }

   # show colorbar if needed
   if (show.colorbar) {
      mdaplot.showColorbar(ps$cgroup, ps$colmap, lab.col = lab.col, lab.cex = lab.cex)
   }

   invisible(ps)
}

#' Create line plot with double y-axis
#'
#' @description
#' \code{mdaplotyy} create line plot for two plot series and uses separate y-axis for each.
#'
#' @param data
#' a matrix or a data.frame with two rows of values.
#' @param type
#' type of the plot ("l" or "b").
#' @param col
#' a color for markers or lines (same as \code{plot} parameter) for each series.
#' @param lty
#' line type for each series (two values)
#' @param lwd
#' line width for each series (two values)
#' @param pch
#' a character for markers (same as \code{plot} parameter) for each series (two values).
#' @param cex
#' scale factor for the markers
#' @param xlim
#' limits for the x axis (if NULL, will be calculated automatically).
#' @param ylim
#' limits for the y axis, either list with two vectors (one for each series) or NULL.
#' @param main
#' an overall title for the plot (same as \code{plot} parameter).
#' @param xlab
#' a title for the x axis (same as \code{plot} parameter).
#' @param ylab
#' a title for each of the two y axis (as a vector of two text values).
#' @param labels
#' a vector with text labels for data points or one of the following: "names", "indices", "values".
#' @param show.labels
#' logical, show or not labels for the data objects.
#' @param lab.cex
#' size for data point labels.
#' @param lab.col
#' color for data point labels.
#' @param show.grid
#' logical, show or not a grid for the plot.
#' @param grid.lwd
#' line thinckness (width) for the grid.
#' @param grid.col
#' line color for the grid.
#' @param xticks
#' values for x ticks.
#' @param xticklabels
#' labels for x ticks.
#' @param xlas
#' orientation of xticklabels.
#' @param ylas
#' orientation of yticklabels (will be applied to both y axes).
#' @param show.legend
#' logical show legend with name of each plot series or not
#' @param legend.position
#' position of legend if it must be shown
#' @param legend
#' values for the legend
#' @param ...
#' other plotting arguments.
#'
#' @details
#' This plot has properties both \code{mdaplot} and \code{mdaplotg}, so when you specify color,
#' line properties etc. you have to do it for both plot series.
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
mdaplotyy <- function(data, type = "l", col = mdaplot.getColors(2), lty = c(1, 1),
   lwd = c(1, 1), pch = (if (type == "b") c(16, 16) else c(NA, NA)), cex = 1,
   xlim = NULL, ylim = NULL, main = attr(data, "name"), xlab = attr(data, "xaxis.name"),
   ylab = rownames(data), labels = "values", show.labels = FALSE, lab.cex = 0.65,
   lab.col = "darkgray", show.grid = TRUE, grid.lwd = 0.5, grid.col = "lightgray",
   xticks = NULL, xticklabels = NULL, xlas = 0, ylas = 0, show.legend = TRUE,
   legend.position = "topright", legend = ylab, ...) {

   if (!(type %in% c("l", "b"))) {
      stop("YY line plot can only be made for type 'l' or 'b'.")
   }

   if (nrow(data) != 2) {
      stop("Matrix with two rows is required to make YY line plot.")
   }

   if (length(ylab) != 2) {
      stop("Y-axis label ('ylab') should be specified for both axes.")
   }

   if (length(col) != 2 || length(lty) != 2 || length(lwd) != 2 || length(pch) != 2) {
      stop("Color and line properties should be specified for both series.")
   }

   # create plot series
   ps1 <- plotseries(mda.subset(data, subset = 1), type = type, labels = labels, col = col[1])
   ps2 <- plotseries(mda.subset(data, subset = 2), type = type, labels = labels, col = col[2])

   # get limits for first series
   xlim <- mdaplot.getXAxisLim(ps1, xlim = xlim, show.labels = show.labels)
   ylim1 <- mdaplot.getYAxisLim(ps1, ylim = ylim[[1]], show.labels = show.labels,
      show.colorbar = TRUE)

   # check and prepare xticklabels
   xticklabels <- mdaplot.getXTickLabels(xticklabels, xticks, NULL)
   xticks <- mdaplot.getXTicks(xticks, xlim, ps1$x_values, type)

   # check and prepare yticklabels
   yticklabels <- mdaplot.getYTickLabels(NULL, NULL, NULL)
   yticks <- mdaplot.getYTicks(NULL, ylim1, ps1$y_values, type)

   # define title and labels
   if (is.null(main)) main <- ps1$name

   par(mar = c(5, 5, 5, 5))

   # make an empty plot with proper limits and axis labels
   mdaplot.plotAxes(xticklabels = xticklabels, yticklabels = yticklabels, xticks = xticks,
      yticks = yticks, xlim = xlim, ylim = ylim1, main = main, xlab = xlab, ylab = ylab[1],
      xlas = xlas, ylas = ylas, show.grid = show.grid, grid.lwd = grid.lwd, grid.col = grid.col
   )

   # show first series
   plotLines(ps1, col = col[1], lty = lty[1], lwd = lwd[1], cex = cex, pch = pch[1], ...)
   if (show.labels) {
      showLabels(ps1, show.excluded = FALSE, col = lab.col, cex = lab.cex)
   }

   # second series
   ylim2 <- mdaplot.getYAxisLim(ps2, ylim = ylim[[2]], show.labels = show.labels)
   par(new = TRUE)
   plot(ps2$x_values, ps2$y_values[1, ], type = type, col = col[2], axes = FALSE, xlab = NA,
      xlim = xlim, ylab = NA, ylim = ylim2, lty = lty[2], lwd = lwd[2], cex = cex, pch = pch[2],
      ...)
   axis(side = 4, las = ylas)
   mtext(side = 4, line = 2, ylab[2], cex = 0.9)

   if (show.labels) {
      showLabels(ps2, show.excluded = FALSE, col = lab.col, cex = lab.cex)
   }

   if (show.legend) {
      mdaplotg.showLegend(legend, col = col, pch = pch, lty = lty, lwd = lwd, cex = 0.9,
         position = legend.position)
   }

   par(mar = c(5.1, 4.1, 4.1, 2.1))
}

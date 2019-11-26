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

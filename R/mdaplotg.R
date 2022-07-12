#' Show legend for mdaplotg
#'
#' @description
#' Shows a legend for plot elements or their groups.
#'
#' @param legend
#' vector with text elements for the legend items
#' @param col
#' vector with color values for the legend items
#' @param pt.bg
#' vector with background colors for the legend items (e.g. for pch = 21:25)
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
#' @param ...
#' other parameters
#'
mdaplotg.showLegend <- function(legend, col, pt.bg = NA, pch = NULL, lty = NULL, lwd = NULL,
   cex = 1, bty = "o", position = "topright", plot = TRUE, ...) {

   # which positions need multiple columns
   onecolpos <- c("topright", "topleft", "bottomright", "bottomleft", "right", "left")
   multcolpos <- c("top", "bottom")

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
   legend(position, legend, col = col, pt.bg = pt.bg, pch = pch, lty = lty, pt.cex = cex, lwd = lwd,
          cex = 0.85, plot = plot, inset = inset, bg = "white", box.lwd = 0.75, box.col = "gray",
          ncol = ncol, ...)

}

#' Prepare data for mdaplotg
#'
#' @param data
#' datasets (in form of list, matrix or data frame)
#' @param type
#' vector with type for dataset
#' @param groupby
#' factor or  data frame with factors - used to split data matrix into groups
#'
#' @return
#' list of datasets
#'
#' The method should prepare data as a list of datasets (matrices or data frames). One list
#' element will be used to create one plot series.
#'
#' If `data` is matrix or data frame and not `groupby` parameter is provided, then every row
#' will be taken as separate set. This option is available only for line or bar plots.
#'
mdaplotg.prepareData <- function(data, type, groupby) {

   # if already a list - remove NULL elements and return
   if (is.list(data) && !is.data.frame(data)) return(data[!sapply(data, is.null)])

   if (is.null(groupby)) {
      # take every row of matrix or data frame as separate group

      if (!all(type %in% c("h", "l", "b"))) {
         stop("Group plot with matrix or data frame can be made only for types 'h', 'l' and 'b'.")
      }

      # add fake row names
      if (is.null(rownames(data))) {
         rownames(data) <- seq_len(nrow(data))
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
   if (!is.data.frame(groupby)) groupby <- as.data.frame(groupby)

   if (!all(unlist(lapply(groupby, is.factor)))) {
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

#' Check mdaplotg parameters and replicate them if necessary
#'
#' @param param
#' A parameter to check
#' @param name
#' name of the parameter (needed for error message)
#' @param is.type
#' function to use for checking parameter type
#' @param ngroups
#' number of groups (plot series)
#'
mdaplotg.processParam <- function(param, name, is.type, ngroups) {

   param <- if (length(param) == 1) rep(param, ngroups) else param

   if (!all(is.type(param))) {
      stop(paste0('Parameter "', name, '" mush be numeric!'))
   }

   if (length(param) != ngroups)
      stop(paste0('Parameter "', name, '" should be specified for each group or one for all!'))

   return(param)
}

#' Create and return vector with legend values
#'
#' @param ps
#' list with plot series
#' @param data.names
#' names of the data sets
#' @param legend
#' legend values provided by user
#'
#' @return
#' vector of text values for the legend
#'
mdaplotg.getLegend <- function(ps, data.names, legend = NULL) {

   # if legend is not provided - get legend items from data names or plotseries names
   if (is.null(legend)) {
      legend <- if (is.null(data.names)) unlist(lapply(ps, function(x) x$data_attrs$name))
         else data.names
   }

   if (is.null(legend)) {
      stop("Can not find values for the legend items.")
   }

   if (length(legend) != length(ps)) {
      stop("Number of values for 'legend' is not the same as number of plot series.")
   }

   return(legend)
}

#' Compute x-axis limits for mdaplotg
#'
#' @param ps
#' list with plotseries
#' @param xlim
#' limits provided by user
#' @param show.excluded
#' logical, will excluded values also be shown
#' @param show.legend
#' will legend be shown on the plot
#' @param show.labels
#' will labels be shown on the plot
#' @param legend.position
#' position of legend on the plot (if shown)
#' @param bwd
#' size of bar for bar plot
#'
#' @return
#' vector with two values
#'
mdaplotg.getXLim <- function(ps, xlim, show.excluded, show.legend, show.labels,
   legend.position, bwd = NULL) {

   # if user provided xlim values - use them
   if (!is.null(xlim)) {
      return(xlim)
   }

   # function which returns xlim values for given plotseries
   f <- function(p) {
      return(
         mdaplot.getXAxisLim(p, xlim = NULL, show.labels = show.labels,
            show.excluded = show.excluded, bwd = bwd)
      )
   }

   # compute limits for all plot series and combine into matrix
   xlim <- matrix(unlist(lapply(ps, f)), ncol = 2, byrow = TRUE)

   # get the smallest of min and larges of max
   xlim <- c(min(xlim[, 1]), max(xlim[, 2]))

   # add extra margins if legend must be shown
   if (show.legend) {

      # calculate margins: (10% of current limits)
      margin <- c(
         (regexpr("left", legend.position) > 0) * -0.1,
         (regexpr("right", legend.position) > 0) * 0.1
      )

      xlim <- xlim + diff(xlim) * margin
   }

   return(xlim)
}

#' Compute y-axis limits for mdaplotg
#'
#' @param ps
#' list with plotseries
#' @param ylim
#' limits provided by user
#' @param show.excluded
#' logical, will excluded values also be shown
#' @param show.legend
#' will legend be shown on the plot
#' @param legend.position
#' position of legend on the plot (if shown)
#' @param show.labels
#' logical, will data ponit labels also be shown
#'
#' @return
#' vector with two values
#'
mdaplotg.getYLim <- function(ps, ylim, show.excluded, show.legend, legend.position, show.labels) {

   # if user provided ylim values - use them
   if (!is.null(ylim)) {
      return(ylim)
   }

   # function which returns ylim values for given plotseries
   f <- function(p) {
      return(
         mdaplot.getYAxisLim(p, ylim = NULL, show.excluded = show.excluded,
            show.labels = show.labels)
      )
   }

   # compute limits for all plot series and combine into matrix
   ylim <- matrix(unlist(lapply(ps, f)), ncol = 2, byrow = TRUE)


   # get the smallest of min and larges of max
   ylim <- c(min(ylim[, 1]), max(ylim[, 2]))

   # add extra margins if legend must be shown
   if (show.legend) {

      # calculate margins: dx and dy
      margin <- c(
         (regexpr("bottom", legend.position) > 0) * -0.1,
         (regexpr("top", legend.position) > 0) * 0.1
      )

      ylim <- ylim + diff(ylim) * margin
   }

   return(ylim)
}

#' Plotting function for several plot series
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
#' from calibration and test set, scatter plot with height and weight values for women and men, and
#' so on.
#'
#' Most of the parameters are similar to \code{\link{mdaplot}}, the difference is described below.
#'
#' The data should be organized as a list, every item is a matrix (or data frame) with data for one
#' set of objects. Alternatively you can provide data as a matrix and use parameter
#' \code{groupby} to create the groups. See tutorial for more details.
#'
#' There is no color grouping option, because color is used to separate the sets. Marker symbol,
#' line style and type, etc. can be defined as a single value (one for all sets) and as a vector
#' with one value for each set.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
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
   data <- mdaplotg.prepareData(data, type, groupby)
   ngroups <- length(data)

   # check if plot.new() should be called first
   if (dev.cur() == 1) plot.new()

   type <- mdaplotg.processParam(type, "type", is.character, ngroups)
   pch <- mdaplotg.processParam(pch, "pch", is.numeric, ngroups)
   lty <- mdaplotg.processParam(lty, "lty", is.numeric, ngroups)
   lwd <- mdaplotg.processParam(lwd, "lwd", is.numeric, ngroups)
   cex <- mdaplotg.processParam(cex, "cex", is.numeric, ngroups)
   opacity <- mdaplotg.processParam(opacity, "opacity", is.numeric, ngroups)
   lab.col <- mdaplotg.processParam(lab.col, "lab.col", mdaplot.areColors, ngroups)

   # check and define colors if necessary
   if (is.null(col)) col <- mdaplot.getColors(ngroups = ngroups, colmap = colmap)
   col <- mdaplotg.processParam(col, "col", mdaplot.areColors, ngroups)

   # get plot data for each group
   ps <- vector("list", ngroups)
   for (i in seq_len(ngroups)) {
      ps[[i]] <- plotseries(data[[i]], type = type[i], col = col[i], opacity = opacity[i],
         labels = labels)
   }

   # get axis limits
   ylim <- mdaplotg.getYLim(ps, ylim, show.excluded, show.legend, legend.position, show.labels)
   xlim <- mdaplotg.getXLim(ps, xlim, show.excluded, show.legend, show.labels, legend.position, bwd)

   # check and prepare xticklabels
   xticklabels <- mdaplot.getXTickLabels(xticklabels, xticks, NULL)
   xticks <- mdaplot.getXTicks(xticks, xlim = xlim)

   # check and prepare yticklabels
   yticklabels <- mdaplot.getYTickLabels(yticklabels, yticks, NULL)
   yticks <- mdaplot.getYTicks(yticks, ylim = ylim)

   # define main title if not provided (either as "name" or as "name" attr of first dataset)
   main <- if (is.null(main)) name else main
   main <- if (is.null(main)) ps[[1]]$name else main

   # define labels for axes
   xlab <- if (is.null(xlab)) attr(ps[[1]]$x_values, "name", exact = TRUE) else xlab
   ylab <- if (is.null(ylab)) attr(ps[[1]]$y_values, "name", exact = TRUE) else ylab

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
   for (i in seq_len(ngroups)) {

      # decide if x values should be forced as group index
      force.x.values <- if (type[i] == "h") c(i, nbarplots) else NA

      # if error bars are shown and i > 1 do not show labels
      show.labels <- if (i > 1 && type[i] == "e") FALSE else show.labels

      # use mdaplot with show.axes = FALSE to create the plot
      mdaplot(ps = ps[[i]], type = type[i], col = col[i], pch = pch[i], lty = lty[i],
              lwd = lwd[i], cex = cex[i], force.x.values = force.x.values, bwd = bwd,
              show.grid = FALSE, show.labels = show.labels, opacity = opacity[i],
              lab.col = lab.col[i], lab.cex = lab.cex, show.axes = FALSE,
              show.excluded = show.excluded, ...
      )
   }

   # show legend if required
   if (show.legend == TRUE) {
      legend <- mdaplotg.getLegend(ps, names(data), legend)
      if (length(legend) != ngroups) {
         stop("Number of values for 'legend' is not the same as number of plot series.")
      }

      lty[type == "p" | type == "h"] <- 0
      pch[type == "l"] <- NA_integer_
      pch[type == "h"] <- 15

      mdaplotg.showLegend(
         legend, col = col, pch = pch, lty = lty, lwd = lwd, cex = 0.85,
         position = legend.position
      )
   }

   return(invisible(ps))
}

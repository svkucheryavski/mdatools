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
mdaplot.formatValues = function(data, round.only = F, digits = 3) {   
   if (round.only == T)
      fdata = round(data, digits)
   else
      fdata = prettyNum(data, digits = digits)
   if (!is.null(dim(data))){
      dim(fdata) = dim(data)
      dimnames(fdata) = dimnames(data)
   }
   
   fdata
}

#' Calculate axes limits
#' 
#' @description
#' Calculates axes limits depending on data values that have to be plotted, 
#' extra plot elements that have to be shown and margins. 
#' 
#' @param x.values
#' a vector with x values.
#' @param y.values
#' a vector or a matrix with y values.
#' @param lower
#' a lower margin for y limits.
#' @param upper
#' an upper margin for y limits.
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
#' 
#' @details
#' Data can be a list with several matrices or just one matrix. The matrices can have single.x
#' configuration, where first column is x values and the others are y values or normal configuration,
#' where every odd column is x values and every even is corresponding y values.
#'
#' @return
#' Returns a list with four limits for the x and y axes.
#'  
mdaplot.getAxesLim = function(x.values, y.values, lower = NULL, upper = NULL,
                              show.colorbar = F, show.lines = F, 
                              legend = NULL, show.legend = F, legend.position = 'topright', 
                              show.labels = F) {
   if (is.null(x.values) || is.null(y.values))
      return(NULL)
   
   if (is.null(show.labels))
      show.labels = FALSE
   
   scale = 0.075 # scale for margins - 7.5% of plot width or height
   xmin = min(x.values)
   xmax = max(x.values)
   if (is.null(lower) || is.null(upper)) {
      ymin = min(y.values, y.values, na.rm = TRUE)
      ymax = max(y.values, y.values, na.rm = TRUE)
   } else {
      ymin = min(y.values, y.values - lower, na.rm = TRUE)
      ymax = max(y.values, y.values + upper, na.rm = TRUE)
   }
   
   if (xmin == xmax) {
      xmin = xmin - 0.1
      xmax = xmax + 0.1
   }
   
   if (ymin == ymax) {
      ymin = ymin - 0.1
      ymax = ymax + 0.1
   }
   
   # correct limits if some lines that have to be shown are outside the data points cloud
   if (is.numeric(show.lines)) {
      if (!is.na(show.lines[1])) {   
         xmax = max(xmax, show.lines[1], na.rm = TRUE)
         xmin = min(xmin, show.lines[1], na.rm = TRUE)
      }
      
      if (!is.na(show.lines[2])) {   
         ymax = max(ymax, show.lines[2], na.rm = TRUE)
         ymax = max(ymax, show.lines[2], na.rm = TRUE)
      }   
   }
   
   # calculate margins: dx and dy
   dx = (xmax - xmin) * scale
   dy = (ymax - ymin) * scale
   
   # define limits with margins
   xlim = c(xmin - dx, xmax + dx)
   ylim = c(ymin - dy, ymax + dy)
   
   # add an extra margin to y limit if colorbar must be shown
   if (show.colorbar == T)
      ylim[2] = ylim[2] + dy * 3
   
   # add extra margins if legend must be shown
   if (show.legend == T && !is.null(legend)) {
      legend.size = mdaplot.showLegend(legend, position = legend.position, plot = F)
      w = legend.size$rect$w
      h = legend.size$rect$h
      
      if (legend.position == 'topleft') {
         xlim[1] = xlim[1] - dx * 0.2
         ylim[2] = ylim[2] + dy * 0.2     
      } else if (legend.position == 'top') {
         ylim[2] = ylim[2] + h + 0.2 * dy  
      } else if (legend.position == 'topright') {
         xlim[2] = xlim[2] + dx * 0.2
         ylim[2] = ylim[2] + dy * 0.2     
      } else if (legend.position == 'bottomleft') {
         xlim[1] = xlim[1] - dx * 0.2
         ylim[1] = ylim[1] - dy * 0.2   
      } else if (legend.position == 'bottom') {
         ylim[1] = ylim[1] - h - 0.2 * dy
      } else if (legend.position == 'bottomright') {
         xlim[2] = xlim[2] + dx * 0.2
         ylim[1] = ylim[1] - dy * 0.2
      }   
   }  
   
   # add extra margins if labels must be shown
   if (show.labels == T) {
      ylim[1] = ylim[1] - dy 
      ylim[2] = ylim[2] + dy
   }   
   
   # return the limits
   lim = list(
      xlim = xlim,
      ylim = ylim
   )   
   
   lim
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
#' which colormap to use ('default', 'gray', or user defined in form c('color1', 'color2', ...)).
#' @param opacity
#' opacity for colors (between 0 and 1)
#' 
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' 
#' @return
#' Returns vector with generated color values
#' 
#' @export
mdaplot.getColors = function(ngroups = 1, cgroup = NULL, colmap = 'default', opacity = 1) {   
   if (length(colmap) == 1) {   
      # colormap is a name      
      if (colmap == 'gray') {   
         # use grayscale colormap
         palette = c("#E8E8E8", "#D6D6D6", "#C4C4C4", "#B2B2B2", "#9A9A9A", "#808080", "#484848", 
                     "#101010")
         
         # if only one color is needed reorder pallete so the black is first
         if (is.null(cgroup) && ngroups == 1)
            palette = palette[length(palette):1]      
      } else if (colmap == 'jet') {
         palette = c("#00007F", "blue", "#007FFF", "cyan",
           "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      } else {   
         # use default colormap for colorbrew
         palette = c("#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", 
                     "#D53E4F")
      }   
   } else {
      # assuming that user colormap has been provided
      if (sum(mdaplot.areColors(colmap) == F) == 0)
         palette = colmap
      else
         stop("Parameter 'colmap' must contains valid color values or name of palette!")
   }   
   
   colfunc = colorRampPalette(palette, space = 'Lab')
   
   if (!is.null(cgroup)) {   
      if (is.factor(cgroup)) {
         nlevels = length(attr(cgroup, 'levels'))
         ngroups = nlevels
         col = colfunc(ngroups)      
         colvals = col[as.numeric(cgroup)]
      } else { 
         cfactor = factor(cgroup)
         nlevels = length(attr(cfactor, 'levels'))
         if (nlevels < 8)
            ngroups = nlevels
         else
            ngroups = 8
         
         col = colfunc(ngroups)      
         cgroup = cut(as.vector(as.numeric(cgroup)), ngroups, labels = 1:ngroups)
         colvals = col[cgroup]
      }
   } else {
      colvals = colfunc(ngroups)      
   }   
   
   if (opacity < 1) {
      opacity = format(as.hexmode(round(255 * opacity)), width = 2);
      colvals = paste(colvals, opacity, sep = '')
   }
   colvals
}

#' Plot grid
#' 
#' @description
#' Shows grid for a plot
#' 
#' @param lwd
#' line width for the grid 
#' 
mdaplot.showGrid = function(lwd = 0.5) {
   # Shows a grid for plots   
   grid(lwd = lwd)   
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
mdaplot.showColorbar = function(cgroup, colmap = 'default', lab.col = 'darkgray', lab.cex = 0.65) {
   # get number of levels for the cgroup
   if (!is.factor(cgroup) && length(unique(cgroup)) > 12) {
         # get colors for 8 groups based on colormap
         col = mdaplot.getColors(ngroups = 12, colmap = colmap)
         ncol = length(unique(col))
      
         # split values to intervals
         cgroupl = levels(cut(as.vector(cgroup), ncol))
      
         # get left and right values for the intervals
         lvals = as.numeric( sub("\\((.+),.*", "\\1", cgroupl) ) 
         rvals = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", cgroupl) ) 
      
         # correct issue with first element
         if (min(cgroup) != lvals[1])
            lvals[1] = min(cgroup)
      
         # combine values and define matrix for labels
         vals = c(lvals, rvals[ncol])
         labels = matrix(0, ncol = 2, nrow = ncol + 1)
   } else {
      if (!is.factor(cgroup))
         cgroup = factor(cgroup)
   
      nlevels = length(attr(cgroup, 'levels'))
      
      # no splitting is needed, just use factors as labels
      col = mdaplot.getColors(ngroups = nlevels, colmap = colmap)
      ncol = length(unique(col))
      vals = levels(cgroup)
      labels = matrix(0, ncol = 2, nrow = ncol)
   }  
   
   # use formatted values as rownames for labels matrix
   rownames(labels) = mdaplot.formatValues(vals)
   
   # get size of the plotting area and calculate size for color bar elements
   lim = par('usr') 
   
   dx = lim[2] - lim[1]
   dy = lim[4] - lim[3]
   
   w = (dx * 0.8)/ncol
   h = dy * 0.015
   
   x = lim[1] + dx * 0.1
   y = lim[4] - (h + 0.1 * h);
   
   # show colorbar and define coordinates for labels
   for (i in 1:ncol)
   {
      rect(x + w * (i - 1), y, x + w * i, y - h, col = col[i], border = NA)      
      labels[i, ] = c(x + w * (i - 1), y - h)
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
#' @param bty
#' border type for the legend
#' @param position
#' legend position ('topright', 'topleft', 'bottomright', 'bottomleft', 'top', 'bottom')
#' @param plot
#' logical, show legend or just calculate and return its size
#'
mdaplot.showLegend = function(legend, col, pch = NULL, lty = NULL, lwd = NULL, bty = 'o', 
                              position = 'topright', plot = T)
{
   # which positions need multiple columns
   onecolpos = c('topright', 'topleft', 'bottomright', 'bottomleft')
   multcolpos = c('top', 'bottom', 'right', 'left')
   
   if (position %in% onecolpos)
      ncol = 1
   else if (position %in% multcolpos)
      ncol = length(legend)
   else   
      stop("Wrong values for 'legend.position' argument!")
   
   # calculate inset values depending on a ration between width and height of a plot
   lim = par('plt') 
   dx = lim[2] - lim[1]
   dy = lim[4] - lim[3]   
   inset = c(0.02, 0.02 * (dx/dy))
   
   # show legend
   legend(position, legend, col = col,  pch = pch, lty = lty, lwd = lwd, plot = plot, 
          cex = 0.8, inset = inset, bg = 'white', box.lwd = 0.75, box.col = 'gray',
          ncol = ncol)   
}


#' Plot labels
#' Shows labels for data elements (points, bars) on a plot. 
#'
#' @param x.values
#' a vector with x-values
#' @param y.values
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
mdaplot.showLabels = function(x.values, y.values, labels, pos = 3, 
                              cex = 0.65, col = 'darkgray', type = NULL) {
   
   if (!is.null(dim(y.values)) && nrow(y.values) > 1 && type != 'p')
      y.values = apply(y.values, 2, max)
   
   # show labels
   x.values = as.vector(x.values)
   y.values = as.vector(y.values)
   
   if (is.numeric(labels))
      labels = mdaplot.formatValues(labels)
   
   if (!any(is.nan(x.values) || is.nan(y.values))) {
      
      if (!is.null(type) && type == 'h') {   
         # show labels properly for bars with positive and negative values
         neg = y.values < 0
         if (sum(neg) > 0)
            text(x.values[neg], y.values[neg], labels[neg], cex = cex, pos = 1, col = col)   
         if (sum(!neg) > 0)
            text(x.values[!neg], y.values[!neg], labels[!neg], cex = cex, pos = 3, col = col)   
      } else {
         text(x.values, y.values, labels, cex = cex, pos = pos, col = col)   
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
#'  If it is needed to show only one line, the other coordinate shall be set to NA.
#'
mdaplot.showLines = function(point, lty = 2, lwd = 0.75, col = rgb(0.2, 0.2, 0.2)) {
   
   if (!is.na(point[2]))
      abline(h = point[2], lty = lty, lwd = lwd, col = col)
   if (!is.na(point[1]))
      abline(v = point[1], lty = lty, lwd = lwd, col = col)
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
mdaplot.showRegressionLine = function(data, lty = 1, lwd = 1, colmap = 'default', col = NULL) {   
   if (is.list(data))
   {   
      ngroups = length(data)
      
      if (is.null(col))         
         col = mdaplot.getColors(ngroups, colmap = colmap)
      else if (length(col) == 1)
         col = rep(col, ngroups)
      
      for (i in 1:ngroups)   
      {   
         x = data[[i]][, 1]
         y = data[[i]][, 2]
         abline(lm(y ~ x), lty = lty, lwd = lwd, col = col[i])
      }          
   }  
   else
   {
      if (is.null(col))         
         col = mdaplot.getColors(1, colmap = colmap)
      
      x = data[, 1]
      y = data[, 2]
      
      abline(lm(y ~ x), lty = lty, lwd = lwd, col = col)                
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
bars = function(x, y, col = NULL, bwd = 0.8, border = NA) {
   # Show bars on a plot area
   #
   # Arguments:
   #   x: vector with x values (where center of bar will be located)
   #   y: vector with y values (height of the bar)
   #   col: color of the bars
   #   bwd: width of the bars
   
   if (length(bwd) == 1)
      bwd = matrix(bwd, ncol = length(x))
  
   if (length(col) != length(y))
      col = rep(col, length.out = length(y))
   
   for (i in 1:length(x))
   {
      rect(x[i] - bwd[i]/2, 0, x[i] + bwd[i]/2, y[i], col = col[i], border = border)      
   }
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
#' 
errorbars = function(x, lower, upper, y = NULL, col = NULL, pch = 16)
{
   e2 = (max(x) - min(x)) / 50
   e1 = (max(x) - min(x)) / (length(x) * 5)
   e = min(e1, e2)
   
   segments(x, y - lower, x, y + upper, col = col)
   segments(x - e, y - lower, x + e, y - lower, col = col)
   segments(x - e, y + upper, x + e, y + upper, col = col)
   
   if (!is.null(y))
      points(x, y, col = col, pch = pch)
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
#' 
mdaplot.plotAxes = function(xticklabels = NULL, yticklabels = NULL, xticks = NULL, yticks = NULL, 
                            lim = NULL, main = NULL, xlab = NULL, ylab = NULL, xlas = 0, ylas = 0) {
   
   # make plot without ticks
   plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = lim$xlim, ylim = lim$ylim, 
        xaxt = 'n', yaxt = 'n')

   # generate x ticks
   if (is.null(xticks)) {
      xticks = axisTicks(lim$xlim, log = FALSE)
      # if xticklabels were provided it is expected that xticks should be integers
      if (!is.null(xticklabels)) 
         xticks = unique(round(xticks[xticks > 0 & xticks <= length(xticklabels)]))
   }
   
   # check if xtikclabels are provided and show x-axis
   if (is.null(xticklabels)) {
      axis(1, at = xticks, las = xlas)
   } else {
      lxtl = length(xticklabels)
      if (lxtl == length(xticks))
         axis(1, at = xticks, labels = xticklabels, las = xlas)
      else if (lxtl > length(xticks) && lxtl >= max(xticks) && min(xticks) > 0)
         axis(1, at = xticks, labels = xticklabels[xticks], las = xlas)
      else
         axis(1, at = xticks, las = xlas)
   }
   
   # generate y ticks
   if (is.null(yticks)) {
      yticks = axisTicks(lim$ylim, log = FALSE)
   }
   
   # check if yticklabels are provided and show y-axis
   if (is.null(yticklabels)) {
      axis(2, at = yticks, las = ylas)
   } else {
      lytl = length(yticklabels)
      if (lytl == length(yticks))
         axis(2, at = yticks, labels = yticklabels, las = ylas)
      else if (lytl > length(yticks) && lytl >= max(xticks) && min(yticks) > 0)
         axis(2, at = yticks, labels = yticklabels[xticks], las = ylas)
      else
         axis(2, at = yticks, las = ylas)
   }
}

prepare.plot.data = function(data, type, xlim, ylim, bwd, show.excluded, show.colorbar, show.labels, 
                             show.lines, show.axes) {
   
   # convert data frame to a matrix
   if (is.data.frame(data)) {
      data = mda.df2mat(data)
   }
   
   # check plot type 
   valid.types = c('p', 'l', 'b', 'h', 'e', 'd')
   if (!(type %in% valid.types))
      stop('Wrong plot type!')
   
   density = FALSE
   if (type == 'd') {
      if (ncol(data) < 2)
         stop('Density plot can be created for dataset with at least two columns!')
      type = 'p'
      density = TRUE
   }
   
   if (is.null(data) || length(data) < 1)
      stop('The provided dataset is empty!')
   
   # correct dimension if one row or one column is provided
   if (!is.null(dim(data))) {
      if (nrow(data) == 1 && ncol(data) > 2 && type == 'p') {
         data = mda.t(data)
         attr(data, 'xaxis.name') = 
            ifelse(is.null(attr(data, 'xaxis.name')), 'Variables', attr(data, 'xaxis.name'))      
      } else if (ncol(data) == 1 && !(type == 'p' || type == 'e')) {
         data = mda.t(data)
         attr(data, 'xaxis.name') = 
            ifelse(is.null(attr(data, 'xaxis.name')), 'Objects', attr(data, 'xaxis.name'))      
      }
   }
   
   # get data attributes
   data.attr = attributes(data)
  
   # if data has no dimension   
   if (is.null(dim(data))) {
      names = names(data)
      if (type == 'p') {
         data = matrix(data, ncol = 1)
         rownames(data) = names         
      } else {
         data = matrix(data, nrow = 1)
         colnames(data) = names
      }
   } 

   data = as.matrix(data)
   
   # prepare vector with x.values for non-scatter plots 
   # we do it here in order to process excluded columns correctly
   excluded.cols = data.attr$exclcols
   if (length(excluded.cols) > 0) {
      excluded.cols = mda.getexclind(excluded.cols, colnames(data), ncol(data))      
      data = data[, -excluded.cols, drop = F]    
   }
   
   x.values = NULL
   if (type != 'p' && ncol(data) > 0) {
      if (is.null(data.attr$xaxis.values)) {
         x.values = 1:ncol(data)
         names(x.values) = colnames(data)
      } else {
         x.values = data.attr$xaxis.values
         if (length(excluded.cols) > 0)
            x.values = x.values[-excluded.cols]
      }
   }


   # process excluded columns, broken.lines is needed to show line plots with 
   # excluded columns correctly
   #broken.lines = FALSE
   
   # process excluded rows
   excluded.rows = data.attr$exclrows
   if (!is.null(excluded.rows) && length(excluded.rows) > 0) {
      excluded.rows = mda.getexclind(excluded.rows, rownames(data), nrow(data))      
      
      # if no excluded rows found do not exclude anything
      if (length(excluded.rows) == 0) {
         excluded.rows = NULL
         show.excluded = FALSE
      }
   } else {
      show.excluded = FALSE
   }
   
   # split data to x and y 
   lower = NULL
   upper = NULL
   xaxis.name = data.attr$xaxis.name
   yaxis.name = data.attr$yaxis.name
   if (type == 'p') {
      if (ncol(data) > 1) {
         x.values = data[, 1]
         attr(x.values, 'name') = colnames(data)[1]
         y.values = data[, 2, drop = F]
         attr(y.values, 'name') = colnames(data)[2]
      } else {
         y.values = data[, 1, drop = F]
         attr(y.values, 'name') = colnames(data)[1]
         x.values = 1:length(y.values)
         attr(x.values, 'name') = 
            ifelse(is.null(data.attr$yaxis.name), 'Objects', data.attr$yaxis.name)
      }
   } else {
      # the x.values have been defined earlier     
      attr(x.values, 'name') = ifelse(is.null(xaxis.name), 'Variables', xaxis.name)
      
      # for bar plot we show only the first row
      if (type == 'h') {
         y.values = data[1, , drop = F]
      }
      
      # for other types we plot all rows as lines/polygons
      if (type == 'l' || type == 'b') {
         y.values = data
      }
      
      # prepare data for errorbar plot
      if (type == 'e') {
         if (nrow(data) < 2)
            stop('For errorbar plot data should contain at least two rows!')
         y.values = data[1, , drop = F]
         attr(y.values, 'name') = rownames(data)[1]        
         lower = data[2, ]
         if (nrow(data) > 2)
            upper = data[3, ]
         else
            upper = lower
      }
      
      attr(y.values, 'name') = ifelse(is.null(yaxis.name), '', yaxis.name)
   }
   
   # if y.values have no dimension correct this
   if ( is.null(dim(y.values)) ){
      y.values = matrix(y.values, nrow = length(y.values))   
   }
   
   # process excluded rows if they have to be shown
   y.values.excludedrows = NULL
   x.values.excludedrows = NULL
   if (length(excluded.rows) > 0 && !(type == 'h' || type == 'e')) {
      y.values.name = attr(y.values, 'name')
      y.values.excludedrows = y.values[excluded.rows, , drop = F]
      y.values = y.values[-excluded.rows, , drop = F]         
      attr(y.values, 'name') = y.values.name
      
      if (type == 'p') {
         x.values.name = attr(x.values, 'name')
         x.values.excludedrows = x.values[excluded.rows]
         x.values = x.values[-excluded.rows]
         attr(x.values, 'name') = x.values.name
      }
   }
   
   if (show.axes) {
      # calculate limits       
      if (show.excluded) {
         lim = mdaplot.getAxesLim(c(x.values, x.values.excludedrows), 
                                  rbind(y.values, y.values.excludedrows),  
                                  lower = lower,
                                  upper = upper,
                                  show.colorbar = show.colorbar, 
                                  show.labels = show.labels,
                                  show.lines = show.lines)
         
      } else {
         lim = mdaplot.getAxesLim(x.values, 
                                  y.values,  
                                  lower = lower,
                                  upper = upper,
                                  show.colorbar = show.colorbar, 
                                  show.labels = show.labels,
                                  show.lines = show.lines)
      }
   } else {
      lim = NULL
   }
   
   res = list()
   res$density = density
   res$lower = lower
   res$upper = upper
   res$x.values = x.values
   res$y.values = y.values
   res$x.values.excludedrows = x.values.excludedrows
   res$y.values.excludedrows = y.values.excludedrows
   res$lim = lim
   res$show.excluded = show.excluded
   res$excluded.rows = excluded.rows
   res$excluded.cols = excluded.cols
   res$data.attr = data.attr
   
   res
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
#' by a user (if NULL it will be created automatically) 
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
#' @param lty  
#' the line type (same as \code{plot} parameter).
#' @param lwd  
#' the line width (thickness) (same as \code{plot} parameter).
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
#' @param show.axes  
#' logical, make a normal plot or show only elements (markers, lines, bars) without axes.
#' @param xticks
#' values for x ticks
#' @param yticks 
#' values for y ticks
#' @param xticklabels  
#' labels for x ticks.
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
#' @param col.excluded
#' color for the excluded objects (rows)
#' @param nbins
#' if scatter density plot is shown, number of segments to split the plot area into 
#' (see also ?smoothScatter)
#' @param colramp
#' Colramp function for density scatter plot
#' @param force.x.values
#' vector with corrected x-values for a bar plot (do not specify this manually)
#' @param opacity
#' opacity for plot colors (value between 0 and 1)
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
mdaplot = function(data = NULL, plot.data = NULL, type = 'p', pch = 16, col = NULL, lty = 1, 
                   lwd = 1, bwd = 0.8,
                   cgroup = NULL, xlim = NULL, ylim = NULL, colmap = 'default', labels = NULL, 
                   main = NULL, xlab = NULL, ylab = NULL, show.labels = F, 
                   show.colorbar = T, show.lines = F, show.grid = T, show.axes = T, 
                   xticks = NULL, yticks = NULL, xticklabels = NULL, yticklabels = NULL, 
                   xlas = 0, ylas = 0, lab.col = 'darkgray', lab.cex = 0.65, 
                   show.excluded = FALSE, col.excluded = '#E0E0E0', nbins = 256,
                   colramp = mdaplot.getColors, force.x.values = NA, opacity = 1, ...)
{   
   if (is.null(plot.data)) {
      plot.data = prepare.plot.data(data, type, xlim, ylim, bwd, show.excluded, show.colorbar, 
                                    show.labels, show.lines, show.axes)
   }
   
   plot.data$lower -> lower
   plot.data$upper -> upper
   plot.data$x.values -> x.values
   plot.data$y.values -> y.values
   plot.data$x.values.excludedrows -> x.values.excludedrows
   plot.data$y.values.excludedrows -> y.values.excludedrows
   plot.data$lim -> lim
   plot.data$excluded.rows -> excluded.rows
   plot.data$excluded.cols -> excluded.cols
   plot.data$data.attr -> data.attr
   plot.data$show.excluded -> show.excluded
   plot.data$density -> density
   
   # if some rows are excluded remove part of values from cgroup   
   if (length(excluded.rows) > 0 && !is.null(cgroup) && length(cgroup) > 1)
      cgroup = cgroup[-excluded.rows]
   
   # processs labels and ticklabels
   if (show.axes == T) {  
      # if some columns are excluded and xticklabels is provided for all columns - exclude some of 
      # the values
      if (!is.null(excluded.cols) && !is.null(xticklabels) && length(xticklabels) == ncol(data))
         xticklabels = xticklabels[-excluded.cols]
      
      nxvals = length(x.values)
      
      # if number of x-values is up to 12 show all of them via xticks
      #if (is.null(xticks) && length(x.values) < 13)
      #   xticks = x.values
      
      # define main label
      if (is.null(main)) {
         if (type == 'h' ) {
            row.names = data.attr$dimnames[[1]]
            if (!is.null(row.names)) {
               if (length(excluded.rows) > 0) {
                  row.names = row.names[-excluded.rows]
               } 
               main = row.names[1]
            }
         } else {
            main = data.attr[["name"]]
         }
      }
      
      # define label for x-axis
      if (is.null(xlab))
         xlab = attr(x.values, 'name')         
      
      # define label for y-axis
      if (is.null(ylab))
         ylab = attr(y.values, 'name')
      
      # correct xticklabels and xticks if only one variable is provided for line or bar plot      
      if (type != 'p' && nxvals == 1) {
         xticks = 1
         if (is.null(xticklabels))
            xticklabels = x.values
      }
      
      # correct x axis limits and bar width if it is a bar plot
      if (type == 'h') {
         if (nxvals == 1) {
            lim$xlim = c(x.values - bwd, x.values + bwd)
         } else {
            bwd = bwd * min(diff(x.values))
            lim$xlim = c(min(x.values) - bwd/2, max(x.values) + bwd/2)
         }
      }
      
      # redefine x axis limits if user specified values
      if (!is.null(xlim))
         lim$xlim = xlim
      
      # redefine y axis limits if user specified values
      if (!is.null(ylim))
         lim$ylim = ylim
      
      # if no x-values were provided to line or bar plot we generate them as integer numbers
      if ( (type %in% c('h', 'l', 'b')) && is.null(xticks) && is.null(data.attr$xaxis.values)) {
         xticks = axisTicks(lim$xlim, log = FALSE)
         xticks = unique(round(xticks[xticks > 0 & xticks <= max(x.values)]))
      }
         
      # make an empty plot with proper limits and axis labels
      mdaplot.plotAxes(xticklabels, yticklabels, xticks, yticks, lim, main, xlab, ylab, xlas, ylas) 
   }
   
   #  get proper colors     
   if (density == FALSE && !is.null(cgroup) && !(type == 'e')) {   
      # show color groups according to cdata values
      col = mdaplot.getColors(cgroup = cgroup, colmap = colmap, opacity = opacity)
   } else {
      # show all points with the same color
      if (is.null(col))
         col = mdaplot.getColors(1, colmap = colmap, opacity = opacity)
      # set cgroup to NULL so method will not try to show color groups or colorbar legend
      cgroup = NULL
   }
   
   # show grid if needed
   if (show.grid == T)
      mdaplot.showGrid()
   
   # correct x.values if they were forced by bwd
   if (is.numeric(force.x.values)) {
      x.values = x.values - bwd/2 + (force.x.values[1] - 0.5) * bwd/force.x.values[2]
      bwd = bwd/force.x.values[2]
   }
   
   # make plot for the data 
   if (type == 'p' && density == FALSE)
      points(x.values, y.values, col = col, pch = pch, lwd = lwd, ...)
   if (type == 'd' && density == TRUE)
      {par(new = T);smoothScatter(x.values, y.values, nbin = nbins, colramp = colramp, axes = F, ann = F,...)}
   else if (type == 'l' || type == 'b')
      matlines(x.values, t(y.values), type = type, col = col, pch = pch, lty = lty, lwd = lwd, ...)
   else if (type == 'h')
      bars(x.values, y.values, col = col, bwd = bwd)
   else if (type == 'e')
      errorbars(x.values, lower, upper, y = y.values, col = col, pch = pch)
   
   # show excluded rows
   if (show.excluded) {
      if (type == 'p')
         points(x.values.excludedrows, y.values.excludedrows, type = type, col = col.excluded, 
                pch = pch, lwd = lwd, ...)
      else if (type == 'l' || type == 'b')
         matlines(x.values, t(y.values.excludedrows), type = type, col = col.excluded, pch = pch, 
                  lty = lty, lwd = lwd, ...)
   }  
   
   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplot.showLines(show.lines)
   
   # show labels if needed
   if (show.labels == T) {
      # compute vector with y-values for labels (line and linescatter plot)
      if (type == 'l' || type == 'b') {
         if (show.excluded == TRUE && length(excluded.rows) > 0){
            y.values.labels = apply(rbind(y.values, y.values.excludedrows), 2, max)
         } else {
            y.values.labels = apply(y.values, 2, max)
         }
      }
      
      labels.excl = NULL
      if (is.null(labels) || (length(labels) == 1 && labels == 'names')) {
         # if labels were not provided, by default use names
         
         # generate names if they are empty
         if (is.null(names(x.values)))
            names(x.values) = 1:length(x.values)
         
         labels.incl = names(x.values)
         if (length(excluded.rows) > 0) 
            labels.excl = names(x.values.excludedrows)
      } else if (length(labels) == 1 && labels == 'values') {
         if (type == 'p' || type == 'h') {
            labels.incl = y.values
            if (length(excluded.rows) > 0)
               labels.excl = y.values.excludedrows
         } else {
            labels.incl = y.values.labels               
         }
      } else if (length(labels) == 1 && labels == 'indices') {
         if (type == 'p') {
            # if indices are provided via attribute - use them
            if (!is.null(data.attr$labels))
               ind = data.attr$labels
            else
               ind = 1:(nrow(y.values) + 
                           ifelse(
                              is.null(y.values.excludedrows), 
                              0, 
                              nrow(y.values.excludedrows)
                           )
               )
            if (length(excluded.rows) > 0) {
               labels.incl = ind[-excluded.rows]
               labels.excl = ind[excluded.rows]
            } else {
               labels.incl = ind
            }
         } else {
            labels.incl = 1:length(x.values)
         }
      } else {
         # labels were provided
         if (length(excluded.rows) > 0) {
            labels.incl = labels[-excluded.rows]
            labels.excl = labels[excluded.rows]
         } else {
            labels.incl = labels
         }
      }
      
      if (type == 'l' || type == 'b') {
         y.values = y.values.labels
      }
      
      if (is.null(labels.incl) || length(labels.incl) != length(x.values))
         stop('No correct labels or label type was provided!')
      
      mdaplot.showLabels(x.values, y.values, labels.incl, type = type, col = lab.col, cex = lab.cex)   
      
      if (show.excluded && !is.null(labels.excl) && type == 'p') {
         mdaplot.showLabels(x.values.excludedrows, y.values.excludedrows, labels.excl, 
                            type = type, col = lab.col, cex = lab.cex)   
      }        
   }      
   
   # show colorbar if needed
   if (!is.null(cgroup) && show.colorbar == T)
      mdaplot.showColorbar(cgroup, colmap, lab.col = lab.col, lab.cex = lab.cex)   
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
#' @param colmap  
#' a colormap to use for coloring the plot items.
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
mdaplotg = function(data, groupby = NULL, type = 'p', pch = 16,  lty = 1, lwd = 1, bwd = 0.8,
                    legend = NULL, xlab = NULL, ylab = NULL, main = NULL, labels = NULL, 
                    ylim = NULL, xlim = NULL, colmap = 'default', legend.position = 'topright', 
                    show.legend = T, show.labels = F, show.lines = F, show.grid = T, 
                    xticks = NULL, xticklabels = NULL, yticks = NULL, yticklabels = NULL, 
                    show.excluded = FALSE, lab.col = 'darkgray', lab.cex = 0.65, xlas = 1, ylas = 1, ...) {   

   # name for the plot (mdain)
   name = NULL
   
   # prepare list with groups of objects
   if (is.matrix(data) || is.data.frame(data)) {
      if (is.null(groupby)) {
         # take every line as a group
         name = attr(data, 'name', exact = TRUE)
         
         if (!all(type %in% c('h', 'l', 'b')))
            stop('Group plot with just one matrix or dataframe is available only for types "h", "l" and "b"!')
         
         # split data into a list of subsets for each group
         data.list = list()
         for (i in 1:nrow(data)) {
            data.list[[i]] = mda.subset(data, subset = i)
         }
         
         # get rownames as legend values
         if (is.null(legend))
            legend = rownames(data)
         
         # redefine the data with list
         data = data.list
      } else {
         # split rows into groups using data frame with factors
         
         if (is.data.frame(groupby))
            is.groupby.factor = all(unlist(lapply(groupby, is.factor)))
         else
            is.groupby.factor = is.factor(groupby)
         
         if (is.groupby.factor == FALSE)
            stop('Parameter "groupby" should be a factor or data frame with several factors')
         
         attrs = mda.getattr(data)
         data = as.data.frame(data)
         
         # in this case if labels = indices generate labels for each case
         data$row.ind = 1:nrow(data)
         data = split(data, groupby)
         for (i in 1:length(data)) {
            row.ind = data[[i]]$row.ind
            data[[i]] = subset(data[[i]], select = -row.ind)
            data[[i]] = mda.setattr(data[[i]], attrs)
            attr(data[[i]], 'exclrows') = which(row.ind %in% attrs$exclrows)
            attr(data[[i]], 'labels') = row.ind
         }
      }
   } 
   
   # check if plot.new() should be called first
   if (dev.cur() == 1)
      plot.new()
   
   ngroups = length(data)
   
   # if plot type is not specified for each group multply default value
   if (length(type) == 1)
      type = rep(type, ngroups)
   else if (length(type) != ngroups)
      stop('Parameter "type" should be specified for each group or be common for all!')
   #else if (!(sum(type == 'h') == 0 | sum(type == 'h') == length(pch)))
   #   stop('Barplot (type = "h") for groups can not be combined with other plots!');
   
   # if marker symbol is not specified for each group multply default value
   if (!is.numeric(pch))
      stop('Parameter "pch" mush be numeric!')   
   else if (length(pch) == 1)
      pch = rep(pch, ngroups)
   else if (length(pch) != ngroups)
      stop('Parameter "pch" should be specified for each group or be common for all!')
   
   # if line type is not specified for each group multply default value
   if (!is.numeric(lty))
      stop('Parameter "lty" mush be numeric!')   
   else if (length(lty) == 1)
      lty = rep(lty, ngroups)
   else if (length(lty) != ngroups)
      stop('Parameter "lty" should be specified for each group or be common for all!')
   
   # if line width is not specified for each group multply default value
   if (!is.numeric(lwd))
      stop('Parameter "lwd" mush be numeric!')   
   else if (length(lwd) == 1)
      lwd = rep(lwd, ngroups)
   else if (length(lwd) != ngroups)
      stop('Parameter "lwd" hould be specified for each group or be common for all!')
   
   # if label color is not specified for each group multply default value
   if (length(lab.col) == 1)
      lab.col = rep(lab.col, ngroups)
   else if (length(lab.col) != ngroups)
      stop('Parameter "lab.col" should be specified for each group or be common for all!')
   
   # get plot data for each group 
   pd = list()
   for (i in 1:ngroups) {
      pd[[i]] = prepare.plot.data(data[[i]], type[i], xlim, ylim, bwd, show.excluded, 
                                  show.colorbar = FALSE, show.labels, 
                                  show.lines, show.axes = TRUE)
   }
   
   # process legend
   if (is.null(legend)) {
      legend = names(data)
      if (is.null(legend)) {
         if (all(type == 'h'))
            legend = unlist(lapply(pd, function(x) {rownames(x$y.values)[1]}))
         else
            legend = unlist(lapply(pd, function(x) {x$data.attr$name}))
      }
   }
   
   # compute limits   
   loc.xlim = matrix(unlist(lapply(pd, function(x) {x$lim$xlim})), ncol = 2, byrow = TRUE)
   loc.ylim = matrix(unlist(lapply(pd, function(y) {y$lim$ylim})), ncol = 2, byrow = TRUE)
   
   loc.xlim = c(min(loc.xlim[, 1]), max(loc.xlim[, 2]))
   loc.ylim = c(min(loc.ylim[, 1]), max(loc.ylim[, 2]))
   lim = list(xlim = loc.xlim, ylim = loc.ylim)

   # change limits to user defined if any
   if (!is.null(ylim))
      lim$ylim = ylim
   
   # correct x limits if bar plot should be shown
   if (any(type == 'h')) {
      dx = 0.35
      lim$xlim = lim$xlim + c(-dx, dx)
   }
   
   if (!is.null(xlim))
      lim$xlim = xlim
   
   col = mdaplot.getColors(ngroups, colmap = colmap)
   
   if (is.null(main) && !is.null(name))
      main = name

   # define main label
   if (is.null(main))
      main = pd[[1]]$data.attr[['name']]
   
   # define label for x-axis
   if (is.null(xlab))
      xlab = attr(pd[[1]]$x.values, 'name', exact = TRUE)         
   
   # define label for y-axis
   if (is.null(ylab))
      ylab = attr(pd[[1]]$y.values, 'name', exact = TRUE)
   
   # make an empty plot with proper limits and axis labels
   mdaplot.plotAxes(xticklabels, yticklabels, xticks, yticks, lim, main, xlab, ylab, xlas, ylas)
   
   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplot.showLines(show.lines)
   
   # show initial exes
   if (show.grid == T)
      mdaplot.showGrid()
   
   nbarplots =  sum(type == 'h')
   
   # make a plot for each group   
   for (i in 1:ngroups) {
      if (type[i] == 'h')
         force.x.values = c(i, nbarplots)
      else
         force.x.values = NA
      
      if (i > 1 && type[i] == 'e')
         show.labels = F
      
      mdaplot(plot.data = pd[[i]], type = type[i], col = col[i], pch = pch[i], lty = lty[i],
              lwd = lwd[i], force.x.values = force.x.values, bwd = bwd,
              labels = labels, show.grid = F, show.labels = show.labels,
              lab.col = lab.col[i], lab.cex = lab.cex, show.axes = FALSE, ...
      )
   }  
   
   
   # show legend if needed
   if (!is.null(legend) && length(legend) > 0 && show.legend == T)
   {
      lty[type == 'p' | type == 'h'] = 0
      pch[type == 'l'] = NA_integer_
      pch[type == 'h'] = 15
      mdaplot.showLegend(legend, col, pch = pch, lty = lty, lwd = lwd, position = legend.position)      
   }   
}

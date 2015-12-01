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
mdaplot.formatValues = function(data, round.only = F, digits = 3)
{   
   if (round.only == T)
      fdata = round(data, digits)
   else
      fdata = prettyNum(data, digits = digits)
   
   dimnames(fdata) = dimnames(data)
   fdata
}

#' Calculate axes limits
#' 
#' @description
#' Calculates axes limits depending on data values that have to be plotted, 
#' extra plot elements that have to be shown and margins. 
#' 
#' @param data
#' a matrix or list with data values (values to be plotted).
#' @param single.x
#' logical, has data matrix (matrices) one column for X and many for Y or not.
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
#' @param xticks
#' values for x ticks to show
#' @param yticks
#' values for yticks to show 
#' 
#' @details
#' Data can be a list with several matrices or just one matrix. The matrices can have single.x
#' configuration, where first column is x values and the others are y values or normal configuration,
#' where every odd column is x values and every even is corresponding y values.
#'
#' @return
#' Returns a list with four limits for the x and y axes.
#'  
mdaplot.getAxesLim = function(data, single.x = T, show.colorbar = F, show.lines = F, 
                              legend = NULL, show.legend = F, legend.position = 'topright', 
                              show.labels = F, xticks = NULL, yticks = NULL)
{   
   if (is.null(data))
      return(NULL)
   
   scale = 0.075 # scale for margins - 7.5% of plot width or height
   
   if (single.x == F )
   {   
      # if in every data matrix there are x and y columns for each series
      
      if (is.list(data))
      {
         xmax = max(data[[1]][, seq(1, ncol(data[[1]]), 2)])
         xmin = min(data[[1]][, seq(1, ncol(data[[1]]), 2)])
         ymax = max(data[[1]][, seq(2, ncol(data[[1]]), 2)])
         ymin = min(data[[1]][, seq(2, ncol(data[[1]]), 2)])
         
         for (i in 1:length(data))
         {
            xmax = max(xmax, data[[i]][, seq(1, ncol(data[[i]]), 2)])
            xmin = min(xmin, data[[i]][, seq(1, ncol(data[[i]]), 2)])
            ymax = max(ymax, data[[i]][, seq(2, ncol(data[[i]]), 2)])
            ymin = min(ymin, data[[i]][, seq(2, ncol(data[[i]]), 2)])         
         }   
      }  
      else
      {   
         xmax = max(data[, seq(1, ncol(data), 2)])
         xmin = min(data[, seq(1, ncol(data), 2)])
         ymax = max(data[, seq(2, ncol(data), 2)])
         ymin = min(data[, seq(2, ncol(data), 2)])
      }
   }
   else
   {
      # if in every data matrix first column is x and the other columns are y values
      
      if (is.list(data))
      {
         xmax = max(data[[1]][, 1])
         xmin = min(data[[1]][, 1])
         ymax = max(data[[1]][, 2:ncol(data[[1]])])
         ymin = min(data[[1]][, 2:ncol(data[[1]])])
         
         for (i in 1:length(data))
         {
            xmax = max(xmax, data[[i]][, 1])
            xmin = min(xmin, data[[i]][, 1])
            ymax = max(ymax, data[[i]][, 2:ncol(data[[i]])])
            ymin = min(ymin, data[[i]][, 2:ncol(data[[i]])])         
         }   
      }  
      else
      {   
         xmax = max(data[, 1])
         xmin = min(data[, 1])
         ymax = max(data[, 2:ncol(data)])
         ymin = min(data[, 2:ncol(data)])
      }      
   }   
   
   # correct limits if some lines that have to be shown are outside the data points cloud
   if (is.numeric(show.lines))
   {
      if (!is.na(show.lines[1])) 
      {   
         xmax = max(xmax, show.lines[1])
         xmin = min(xmin, show.lines[1])
      }
      
      if (!is.na(show.lines[2])) 
      {   
         ymax = max(ymax, show.lines[2])
         ymax = max(ymax, show.lines[2])
      }   
   }
  
   # correct limits for manual tick values
   if (!is.null(xticks))
   {
      xmax = max(xmax, xticks)
      xmin = min(xmin, xticks)
   }   

   if (!is.null(yticks))
   {
      ymax = max(ymax, yticks)
      ymin = min(ymin, yticks)
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
   if (show.legend == T && !is.null(legend))
   {
      legend.size = mdaplot.showLegend(legend, position = legend.position, plot = F)
      w = legend.size$rect$w
      h = legend.size$rect$h
      
      if (legend.position == 'topleft')
      {
         xlim[1] = xlim[1] - dx * 0.2
         ylim[2] = ylim[2] + dy * 0.2     
      }   
      else if (legend.position == 'top')
      {
         ylim[2] = ylim[2] + h + 0.2 * dy  
      }   
      else if (legend.position == 'topright')
      {
         xlim[2] = xlim[2] + dx * 0.2
         ylim[2] = ylim[2] + dy * 0.2     
      }   
      else if (legend.position == 'bottomleft')
      {
         xlim[1] = xlim[1] - dx * 0.2
         ylim[1] = ylim[1] - dy * 0.2   
      }   
      else if (legend.position == 'bottom')
      {
         ylim[1] = ylim[1] - h - 0.2 * dy
      }   
      else if (legend.position == 'bottomright')
      {
         xlim[2] = xlim[2] + dx * 0.2
         ylim[1] = ylim[1] - dy * 0.2
      }   
   }  
   
   # add extra margins if labels must be shown
   if (show.labels == T)
   {
      ylim[1] = ylim[1] - dy 
      ylim[2] = ylim[2] + dy
   }   
      
   # return the limits
   lim = list(
      xlim = xlim,
      ylim = ylim
   )   
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
#' 
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' 
#' @return
#' Returns vector with generated color values
#' 
mdaplot.getColors = function(ngroups = 1, cgroup = NULL, colmap = 'default')
{   
   if (length(colmap) == 1)
   {   
      # colormap is a name      
      if (colmap == 'gray')
      {   
         # use grayscale colormap
         palette = c("#E8E8E8", "#D6D6D6", "#C4C4C4", "#B2B2B2", "#9A9A9A", "#808080", "#484848", 
                     "#101010")
         
         # if only one color is needed reorder pallete so the black is first
         if (is.null(cgroup) && ngroups == 1)
            palette = palette[length(palette):1]      
      }   
      else
      {   
         # use default colormap for colorbrew
         palette = c("#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", 
                     "#D53E4F")
      }   
   }
   else
   {
      # assuming that user colormap has been provided
      if (sum(mdaplot.areColors(colmap) == F) == 0)
         palette = colmap
      else
         stop("Parameter 'colmap' must contains valid color values or name of palette!")
   }   
   
   colfunc = colorRampPalette(palette, space = 'Lab')
   
   if (!is.null(cgroup))
   {   
      if (is.factor(cgroup))
      {
         nlevels = length(attr(cgroup, 'levels'))
         ngroups = nlevels
         col = colfunc(ngroups)      
         colvals = col[as.numeric(cgroup)]
      }   
      else
      { 
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
   }
   else
   {
      colvals = colfunc(ngroups)      
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
mdaplot.showGrid = function(lwd = 0.5)
{
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
#' 
mdaplot.showColorbar = function(cgroup, colmap = 'default')
{
   # get number of levels for the cgroup
   cfactor = factor(cgroup)
   nlevels = length(attr(cfactor, 'levels'))
   
   if (nlevels > 8)
   {
      # more than 8 levels, treat cgroup as non factor variable
      
      # get colors for 8 groups based on colormap
      col = mdaplot.getColors(ngroups = 8, colmap = colmap)
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
   }
   else
   {      
      # no splitting is needed, just use factors as labels
      col = mdaplot.getColors(ngroups = nlevels, colmap = colmap)
      ncol = length(unique(col))
      vals = levels(cfactor)
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
   mdaplot.showLabels(labels, pos = 1)
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
#' @param data
#' data matrix with coordinates of the points (x, y) 
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
mdaplot.showLabels = function(data, pos = 3, cex = 0.65, col = 'darkgray', type = NULL)
{

   if (is.null(rownames(data)))
      rownames(data) = 1:nrow(data)
   
   # show labels
   if (!any(is.nan(data)))
   {
      if (!is.null(type) && type == 'h')
      {   
         # show labels properly for bars with positive and negative values
         neg = data[, 2] < 0
         if (sum(neg) > 0)
            text(data[neg, 1], data[neg, 2], rownames(data)[neg], cex = cex, pos = 1, col = col)   
         if (sum(!neg) > 0)
            text(data[!neg, 1], data[!neg, 2], rownames(data)[!neg], cex = cex, pos = 3, col = col)   
      }
      else 
         text(data[, 1], data[, 2], rownames(data), cex = cex, pos = pos, col = col)   
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
mdaplot.showLines = function(point, lty = 2, lwd = 0.75, col = rgb(0.2, 0.2, 0.2))
{
   
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
mdaplot.showRegressionLine = function(data, lty = 1, lwd = 1, colmap = 'default', col = NULL)
{   
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
bars = function(x, y, col = NULL, bwd = 0.8, border = NA)
{
   # Show bars on a plot area
   #
   # Arguments:
   #   x: vector with x values (where center of bar will be located)
   #   y: vector with y values (height of the bar)
   #   col: color of the bars
   #   bwd: width of the bars
   
   if (length(bwd) == 1)
      bwd = matrix(bwd, ncol = length(x))
   
   for (i in 1:length(x))
   {
      rect(x[i] - bwd[i]/2, 0, x[i] + bwd[i]/2, y[i], col = col, border = border)      
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
   
   segments(x, lower, x, upper, col = col)
   segments(x - e, lower, x + e, lower, col = col)
   segments(x - e, upper, x + e, upper, col = col)
   
   if (!is.null(y))
      points(x, y, col = col, pch = pch)
}  

#' Create axes plane
#' 
#' @description
#' Creates an empty axes plane for given parameters
#' 
#' @param xticks
#' tick values for x axis
#' @param xticklabels
#' labels for x ticks
#' @param yticks
#' tick values for y axis
#' @param yticklabels
#' labels for y ticks
#' @param lim
#' vector with limits for x and y axis
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' 
mdaplot.plotAxes = function(xticks, xticklabels, yticks, yticklabels, lim, main, xlab, ylab)
{
   if ((!is.null(xticklabels) || !is.null(xticks)) && (!is.null(yticklabels) || !is.null(yticks)))
   {
      # tick labels for both x and y axes are provided
      
      if (is.null(yticklabels))
         yticklabels = yticks
      
      if (is.null(xticklabels))
         xticklabels = xticks
      
      plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = lim$xlim, ylim = lim$ylim, 
           xaxt = 'n', yaxt = 'n')
      axis(1, at = xticks, labels = xticklabels)
      axis(2, at = yticks, labels = yticklabels)
   }
   else if (!is.null(xticklabels))
   {
      # tick labels for x only provided
      
      if (is.null(xticklabels))
         xticklabels = xticks
      
      plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = lim$xlim, ylim = lim$ylim, 
           xaxt = 'n')
      axis(1, at = xticks, labels = xticklabels)
   }  
   else if (!is.null(yticklabels) || !is.null(yticks))
   {
      # tick labels for y only provided
     
      if (is.null(yticklabels))
         yticklabels = yticks

      plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = lim$xlim, ylim = lim$ylim, 
           yaxt = 'n')
      axis(2, at = yticks, labels = yticklabels)
   }   
   else
   {   
      # no tick labels provided
      plot(0, 0, type = 'n', main = main, xlim = lim$xlim, ylim = lim$ylim, xlab = xlab, ylab = ylab)
   }   
}

#' Plotting function for a single set of objects
#'
#' @description
#' \code{mdaplot} is used to make a scatter, line or a bar plot for one set of data objects.
#' 
#' @param data  
#' a matrix with data values.
#' @param type  
#' type of the plot ('p', 'l', 'b', 'h').
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
#' a vector with text labels for data points (if NULL, row names will be used).
#' @param single.x  
#' logical, is first column of data matrix used for x values (or every odd)?
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
#' tick values for x axis.
#' @param xticklabels  
#' labels for x ticks.
#' @param yticks  
#' tick values for y axis.
#' @param yticklabels  
#' labels for y ticks.
#' @param lab.col  
#' color for data point labels.
#' @param lab.cex  
#' size for data point labels.
#' @param ...  
#' other plotting arguments.
#' 
#' @details 
#' Most of the parameters are similar to what are used with standard \code{plot} function. The 
#' differences are described below.
#'   
#' The function makes a plot of one set of objects. It can be a set of points (scatter plot), 
#' bars, lines or scatter-lines. The data is organized as a matrix. For scatter and bar plot only 
#' two columns are needed, but to plot set of lines there are two ways to organise the data matrix.
#' 
#' If parameter \code{single.x = T} it means first column of the matrix contains x values, as the 
#' other columns contain y values (one column - one line, e.g. [x y1 y2 y3]). If \code{single.x = F}
#' it means every odd column of matrix contains x values and every even column - y values of a
#' line (e.g. [x1 y1 x2 y2 x3 y3]). This option is needed if lines have different x values.
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
#' y axis, see examples for more details.
#'
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#' 
#' @seealso 
#' \code{\link{mdaplotg}} - to make for several sets of data objects (groups of objects).
#'
#' @examples
#' library(mdatools)
#' ### Examples of using mdaplot function
#'  
#' ## 1. Make a line plot for y = x^2 with different x values for each line
#' 
#' x1 = seq(-10, 10, length = 100) 
#' y1 = x1^2
#' x2 = seq(-9, 9, length = 100) + 1
#' y2 = x2^2 * 1.2
#' x3 = seq(-8, 8, length = 100) + 2
#' y3 = x3^2 * 1.4
#' 
#' mdaplot(cbind(x1, y1, x2, y2, x3, y3), type = 'l', single.x = FALSE)
#' 
#' ## 2. Change tick values and labels for x axis
#' mdaplot(cbind(x1, y1, x2, y2, x3, y3), type = 'l', single.x = FALSE,
#'         xticks = c(-8, 0, 8), xticklabels = c('Negative', 'Zero', 'Positive'))
#' 
#' 
#' ## 3. Make a line plot of the spectra with coloring by concentration of first component
#' ## using different color maps
#' 
#' data(simdata)
#' 
#' # first column of matrix val will contain the wavelength
#' # and the last columns - spectra from calibration set
#' val = cbind(simdata$wavelength, t(simdata$spectra.c))
#' 
#' # concentration will be used for color groups
#' c1 = simdata$conc.c[, 1] 
#' 
#' par(mfrow = c(2, 2))
#' mdaplot(val, type = 'l', cgroup = c1)
#' mdaplot(val, type = 'l', cgroup = c1, colmap = 'gray', show.colorbar = FALSE)
#' mdaplot(val, type = 'l', cgroup = c1, colmap = c('red', 'green'))
#' mdaplot(val, type = 'l', cgroup = c1, colmap = c('#ffff00', '#00ffff'))
#' par(mfrow = c(1, 1))
#' 
#' ## 3. Show scatter plots from spectral data with color groups and other parameters
#' ## see how limits are adjusted if show.lines option is used
#' 
#' nobj = 30
#' # concentration is used for color groups
#' c1 = simdata$conc.c[1:nobj, 1]
#' # x values are absorbance of waveband 100, y values - absorbance for waveband 110
#' pdata = cbind(
#'    simdata$spectra.c[1:nobj, 100, drop = FALSE], 
#'    simdata$spectra.c[1:nobj, 110, drop = FALSE])
#' 
#' par(mfrow = c(2, 2))
#' mdaplot(pdata, cgroup = c1, main = 'Spectra', show.lines = c(0.1, 0.06))
#' mdaplot(pdata, cgroup = c1, main = 'Spectra', show.lines = c(0.06, 0.01))
#' mdaplot(pdata, cgroup = c1, xlab = '309 nm', ylab = '319 nm', main = 'Spectra', show.labels = TRUE)
#' mdaplot(pdata, col = 'red', pch = 17, show.labels = TRUE, labels = paste('Obj', 1:nobj))
#' par(mfrow = c(1, 1))
#' 
#' @export
mdaplot = function(data, type = 'p', pch = 16, col = NULL, lty = 1, lwd = 1, bwd = 0.8,
                   cgroup = NULL, xlim = NULL, ylim = NULL, colmap = 'default', labels = NULL, 
                   main = NULL, xlab = NULL, ylab = NULL, single.x = T, show.labels = F, 
                   show.colorbar = T, show.lines = F, show.grid = T, show.axes = T, 
                   xticks = NULL, xticklabels = NULL, yticks = NULL, yticklabels = NULL, 
                   lab.col = 'darkgray', lab.cex = 0.65, ...)
{   
   data = as.matrix(data)
   
   if (is.null(dim(data)) || nrow(data) < 1)
   {   
      warning('Data matrix is empty!')
      return
   }
   
   if (show.axes == T)
   {  

      if (is.null(xticks) && !is.null(xticklabels))
         xticks = 1:length(xticklabels)
      
      if (is.null(yticks) && !is.null(yticklabels))
         yticks = 1:length(yticklabels)
      
      # calculate limits and get proper colors      
      if (!is.null(cgroup))
      {   
         # show color groups according to cdata values
         lim = mdaplot.getAxesLim(data, show.colorbar = show.colorbar, show.labels = show.labels,
                                  show.lines = show.lines, single.x = single.x)
      }
      else   
      {
         # show all points with the same color
         lim = mdaplot.getAxesLim(data, show.lines = show.lines, show.labels = show.labels,
                                  single.x = single.x, xticks = xticks, yticks = yticks)
      }
      
      # Correct x axis limits if it is a bar plot
      if (type == 'h')
         lim$xlim = c(1 - bwd/2 - 0.1, nrow(data) + bwd/2 + 0.1)
      
      # Redefine x axis limits if user specified values
      if (!is.null(xlim))
         lim$xlim = xlim
      
      # Redefine y axis limits if user specified values
      if (!is.null(ylim))
         lim$ylim = ylim
      
      # if user did not specified x axis label use 1st column name for it
      if (is.null(xlab))
         xlab = colnames(data)[1]
      
      # if user did not specified y axis label use 2nd column name for it
      if (is.null(ylab))
         ylab = colnames(data)[2]
      # make an empty plot with proper limits and axis labels
      if (!is.list(data) && nrow(data) < 10 && type != 'p' && single.x == T && is.null(xticklabels))
      {   
         xticks = data[, 1]
         xticklabels = data[, 1]
      }
      
      # if number of variables is small set up proper tick labels for x      
      if (nrow(data) < 10 && single.x == T && type != 'p' && is.null(xticklabels))
      {
         xticks =  data[, 1];
         xticklabels = data[, 1]
      }
      
      # make an empty plot with proper limits and axis labels
      mdaplot.plotAxes(xticks, xticklabels, yticks, yticklabels, lim, main, xlab, ylab)
   }
   
   #  get proper colors      
   if (!is.null(cgroup))
   {   
      # show color groups according to cdata values
      col = mdaplot.getColors(cgroup = cgroup, colmap = colmap)
   }
   else   
   {
      # show all points with the same color
      if (is.null(col))
         col = mdaplot.getColors(1, colmap = colmap)
   }
   
   # get x and y values from data
   if (type == 'e')
   {
      if (ncol(data) != 4)
         stop('Number of columns in data for errorbar plot should be equal to four!')   
      
      x = data[, 1, drop = F]
      y = data[, 2, drop = F]      
      lower = data[, 3, drop = F]      
      upper = data[, 4, drop = F]      
   }   
   else if (single.x == T)
   {   
      x = data[, 1, drop = F]
      y = data[, 2:ncol(data), drop = F]
   }
   else   
   {
      x = data[, seq(1, ncol(data), 2), drop = F]
      y = data[, seq(2, ncol(data), 2), drop = F]
   }
   
   # set up labels
   if (!is.null(labels) && !is.logical(labels))
      rownames(data) = labels
   
   # show grid if needed
   if (show.grid == T)
      mdaplot.showGrid()
   
   # make plot for the data 
   if (type == 'p')
      points(x, y, type = type, col = col, pch = pch, lwd = lwd, ...)
   else if (type == 'l' || type == 'b')
      matlines(x, y, type = type, col = col, pch = pch, lty = lty, lwd = lwd, ...)
   else if (type == 'h')
      bars(x, y, col = col, bwd = bwd)
   else if (type == 'e')
      errorbars(x, lower, upper, y, col = col)
   
   # show labels if needed
   if ((show.labels == T || !is.null(labels)) && type != 'e')
      mdaplot.showLabels(data, type = type, col = lab.col, cex = lab.cex)   
   
   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplot.showLines(show.lines)
   
   # show colorbar if needed
   if (!is.null(cgroup) && show.colorbar == T)
      mdaplot.showColorbar(cgroup, colmap)   
}

#' Plotting function for several sets of objects
#'
#' @description 
#' \code{mdaplotg} is used to make scatter, line or bar plots or their combination for several sets 
#' of objects.
#'
#' @param data  
#' a list with data values (see details below).
#' @param type  
#' type of the plot ('p', 'l', 'b', 'h').
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
#' a matrix with text labels for data points (if NULL, row names for each group will be used).
#' @param ylim  
#' limits for the y axis (if NULL, will be calculated automatically).
#' @param xlim  
#' limits for the x axis (if NULL, will be calculated automatically).
#' @param colmap  
#' a colormap to use for coloring the plot items.
#' @param legend.position  
#' position of the legend ('topleft', 'topright', 'top', 'bottomleft', 'bottomright', 'bottom').
#' @param single.x  
#' logical, is first column of data matrix used for x values (or every odd)?
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
#' @param lab.col  
#' color for data point labels.
#' @param lab.cex  
#' size for data point labels.
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
#' See examples for details.
#' 
#' There is no color grouping option, because color is used to separate the sets. Marker symbol, 
#' line style and type, etc. can be defined as a single value (one for all sets) and as a vector 
#' with one value for each set.
#' 
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#' 
#' @examples 
#' ### Examples of using mdaplotg
#' 
#' ## 1. Show scatter plot for with height and weight values 
#' ## for males and females
#' 
#' data(people)
#' 
#' sex = people[, 'Sex']
#' fdata = cbind(people[sex == 1, 'Height'], people[sex == 1, 'Weight'])
#' mdata = cbind(people[sex == -1, 'Height'], people[sex == -1, 'Weight'])
#' colnames(fdata) = colnames(mdata) = c('Height', 'Weight')
#' pdata = list(fdata, mdata)
#' 
#' par(mfrow = c(2, 2))
#' mdaplotg(pdata, legend = c('female', 'male'))
#' mdaplotg(pdata, legend = c('female', 'male'), pch = c(17, 19), colmap = c('green', 'orange'))
#' mdaplotg(pdata, legend = c('female', 'male'), legend.position = 'topleft', colmap = 'gray')
#' mdaplotg(pdata, legend = c('female', 'male'), legend.position = 'topleft', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 2. Change tick values for x axis
#' mdaplotg(pdata, legend = c('female', 'male'), xticks = c(160, 175, 190), 
#'          xticklabels = c('Short', 'Medium', 'Tall'))
#' 
#' ## 3. Show line plots with spectra from calibration and test set of the Simdata
#' 
#' data(simdata)
#' 
#' cdata = cbind(simdata$wavelength, t(simdata$spectra.c))
#' tdata = cbind(simdata$wavelength, t(simdata$spectra.t))
#' pdata = list(cdata, tdata)
#' 
#' par(mfrow = c(2, 2))
#' mdaplotg(pdata, type = 'l', legend = c('cal', 'test'))
#' mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), xlab = 'Wavelength, nm', ylab = 'Log(1/R)')
#' mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), lty = c(3, 2))
#' mdaplotg(pdata, type = 'l', legend = c('cal', 'test'), lwd = 2, colmap = c('green', 'orange'))
#' par(mfrow = c(1, 1))
#' 
#' @importFrom graphics abline axis grid hist lines matlines par plot points rect segments text
#' 
#' @export
mdaplotg = function(data, type = 'p', pch = 16,  lty = 1, lwd = 1, bwd = 0.8,
                    legend = NULL, xlab = NULL, ylab = NULL, main = NULL, labels = NULL, 
                    ylim = NULL, xlim = NULL, colmap = 'default', legend.position = 'topright', 
                    single.x = T, show.legend = T, show.labels = F, show.lines = F, show.grid = T, 
                    xticks = NULL, xticklabels = NULL, yticks = NULL, yticklabels = NULL, 
                    lab.col = 'darkgray', lab.cex = 0.65, ...)
{   
   #   If data for bar plot is organized as a list, the X values should be contioneous,
   #   e.g. 1:20 for first item, 21:35, for second, 36:42 for third, etc.
      
   # get number of groups
   if (!is.list(data))
   {   
      if (single.x == T)
         ngroups = ncol(data) - 1
      else
         ngroups = round(ncol(data)/2)
   }
   else
   {   
      ngroups = length(data)
   }

   if (is.null(xticks) && !is.null(xticklabels))
      xticks = 1:length(xticklabels)
   
   if (is.null(yticks) && !is.null(yticklabels))
      yticks = 1:length(yticklabels)
   
   # calculate limits and get colors
   lim = mdaplot.getAxesLim(data, show.lines = show.lines, single.x = single.x, 
                            show.legend = show.legend, legend = legend, show.labels = show.labels,
                            legend.position = legend.position, xticks = xticks, yticks = yticks)

   if (!is.null(ylim))
      lim$ylim = ylim
   
   if (!is.null(xlim))
      lim$xlim = xlim
   
   col = mdaplot.getColors(ngroups, colmap = colmap)
   
   # if plot type is not specified for each group multply default value
   if (length(type) == 1)
      type = rep(type, ngroups)
   else if (length(type) != ngroups)
      stop('Parameter "type" should be specified for each group or be common for all!')
   else if (!(sum(type == 'h') == 0 | sum(type == 'h') == length(pch)))
      stop('Barplot (type = "h") for groups can not be combined with other plots!');
      
   # Correct x axis limits if it is a bar plot
   if (type[1] == 'h')
   {
      if (is.list(data))      
         lim$xlim = c(1 - bwd/2 - 0.1, lim$xlim[2] + bwd/2 + 0.1)
      else   
         lim$xlim = c(1 - bwd/2 - 0.1, nrow(data) + bwd/2 + 0.1)
   }   
   
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
   
   # define axis labels as data column names if they are not specified
   if (is.null(xlab) && !is.null(colnames(data[[1]])))
       xlab = colnames(data[[1]])[1]
   
   if (is.null(ylab) && !is.null(colnames(data[[1]])))
      ylab = colnames(data[[1]])[2]

   # if data is small set up proper x tick labels
   if (!is.list(data) && nrow(data) < 10 && type != 'p' && single.x == T && is.null(xticklabels))
   {   
      xticks = data[, 1]
      xticklabels = data[, 1]
   }
   
   # make an empty plot with proper limits and axis labels
   mdaplot.plotAxes(xticks, xticklabels, yticks, yticklabels, lim, main, xlab, ylab)
   
   if (show.grid == T)
      mdaplot.showGrid()
   
   # make a plot for each group   
   if (is.list(data))
   {   
      for (i in 1:ngroups)
      {
         if (!is.null(labels) && type[i] != 'e')
         {
            show(labels)
            show(type)
            if (is.list(labels))
               slabels = labels[[i]]
            else
               slabels = labels[, i]
         }
         else
         {
            slabels = NULL
         }   
         
         mdaplot(data[[i]], type = type[i], col = col[i], pch = pch[i], lty = lty[i],
                 lwd = lwd[i],
                 labels = slabels, show.grid = F, show.axes = F, show.labels = show.labels,
                 lab.col = lab.col, lab.cex = lab.cex)
      }
   }
   else
   {
      if (single.x == T)
      {   
         gbwd = bwd/ngroups            
         for (i in 1:ngroups)
         {
            if (type[1] == 'h')
            {
               x = (1:nrow(data)) + gbwd * (i - 1) - gbwd/2 * (ngroups - 1)   
            }  
            else
            {   
               x = data[, 1, drop = F]
            }
            
            y = data[, i + 1, drop = F]
            
            mdaplot(cbind(x, y), type = type[i], col = col[i], pch = pch[i], lty = lty[i],
                    lwd = lwd[i],
                    bwd = 0.9 * gbwd, labels = labels[, i], show.labels = show.labels,
                    show.grid = F, show.axes = F, lab.col = lab.col, lab.cex = lab.cex)
         }
      }         
      else
      {
         for (i in 1:ngroups)
         {
            x = data[, 2 * i - 1, drop = F]
            y = data[, 2 * i, drop = F]
            mdaplot(cbind(x, y), type = type[i], col = col[i], pch = pch[i], lty = lty[i],
                    lwd = lwd[i],
                    labels = labels[, i], show.grid = F, show.axes = F, show.labels = show.labels,
                    lab.col = lab.col, lab.cex = lab.cex)
         }
      }
   }  
   
   # show lines if needed
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplot.showLines(show.lines)
   
   # show legend if needed
   if (!is.null(legend) && length(legend) > 0 && show.legend == T)
   {
      lty[type == 'p' | type == 'h'] = 0
      pch[type == 'l'] = NA_integer_
      pch[type == 'h'] = 15
      mdaplot.showLegend(legend, col, pch = pch, lty = lty, lwd = lwd, position = legend.position)      
   }   
}

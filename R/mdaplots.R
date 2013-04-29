# class and methods for plotting #
mdaplots = function(scores, loadings, data, ...) UseMethod("mdaplots")

mdaplots.getAxesLim = function(data, show.colorbar = F, multi.y = F, show.legend = F, show.lines = F)
{
   scale = 0.05 # 5%
   if (multi.y == F )
   {   
      if (is.list(data))
      {
         xmax = max(data[[1]][, 1])
         xmin = min(data[[1]][, 1])
         ymax = max(data[[1]][, 2])
         ymin = min(data[[1]][, 2])
         
         for (i in 1:length(data))
         {
            xmax = max(xmax, data[[i]][, 1])
            xmin = min(xmin, data[[i]][, 1])
            ymax = max(ymax, data[[i]][, 2])
            ymin = min(ymin, data[[i]][, 2])         
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
      if (is.list(data))
      {
      }  
      else
      {   
         xmax = max(data[, 1])
         xmin = min(data[, 1])
         ymax = max(data[, 2:ncol(data)])
         ymin = min(data[, 2:ncol(data)])
      }      
   }   
   
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
   
   dx = (xmax - xmin) * scale
   dy = (ymax - ymin) * scale
   
   xlim = c(xmin - dx, xmax + dx)
   ylim = c(ymin - dy, ymax + dy)
   
   if (show.colorbar == T)
      ylim[2] = ylim[2] + dy * 2
   
   if (show.legend == T)
   {
      xlim[2] = xlim[2] + dx
      ylim[2] = ylim[2] + dy
   }   
   
   lim = list(
      xlim = xlim,
      ylim = ylim
   )   
}

mdaplots.showColorbar = function(col, cgroup)
{
   # TODO(svk): format data for colorbar legend
   
   col = mdaplots.getColors(length(unique(col)))
   ncol = length(col)
   cgroup = levels(cut(as.vector(cgroup), ncol))
   lvals = as.numeric( sub("\\((.+),.*", "\\1", cgroup) ) 
   rvals = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", cgroup) )  
   
   lim = par('usr') 
   dx = lim[2] - lim[1]
   dy = lim[4] - lim[3]
   
   w = (dx * 0.8)/ncol
   h = dy * 0.015
   
   x = lim[1] + dx * 0.1
   y = lim[4] - (h + 0.1 * h);
   
   for (i in 1:length(col))
   {
      rect(x + w * (i - 1), y, x + w * i, y - h, col = col[i], border = NA)
      text(x + w * (i - 1), y - h, format(lvals[i], digits = 3), 
           cex = 0.6, pos = 1, col = 'gray')
   }	
   text(x + w * i, y - h, rvals[i], cex = 0.6, pos = 1, col = 'gray')
}

mdaplots.getColors = function(n = 1, cgroup = NULL)
{
   rgbcol = c(213, 62, 79, 
              244, 109, 67, 
              253, 174, 97, 
              254, 224, 139, 
              230, 245, 152, 
              171, 221, 164,    
              102, 194, 165,
              50, 136, 189)
   rgbcol = matrix(rgbcol, ncol = 3, byrow = T)
   ngroups = nrow(rgbcol)
   col8 = rgb(rgbcol[, 1], rgbcol[, 2], rgbcol[, 3], maxColorValue = 255)
   col8 = col8[seq(length(col8), 1, -1)]
   
   if (!is.null(cgroup))
   {   
      cgroup = cut(as.vector(cgroup), ngroups, labels = 1:ngroups)
      col = col8[cgroup]
   }
   else
   {   
      if (length(n) == 1)
      {   
         if (n == 1)
            col = col8[1]
         else if (n == 2)
            col = col8[c(1, 8)]
         else if (n == 3)
            col = col8[c(1, 5, 8)]
         else
            col = col8
      }
   }
   
   return (col)
}

mdaplots.showLegend = function(legend, col, pch = NULL, lty = NULL, pos = 'topright')
{
   legend(pos, legend, col = col,  pch = pch, lty = lty,
          cex = 0.8, inset = 0.01, bg = 'white')
}

mdaplots.showLabels = function(data)
{
   if (is.null(rownames(data)))
      rownames(data) = 1:nrow(data)
   
   text(data[, 1], data[, 2], rownames(data), cex = 0.6, pos = 3, col = 'gray')   
}  

mdaplots.showLines = function(point)
{
   if (!is.na(point[2]))
      abline(h = point[2], lty = 2, col = 'darkgray')
   if (!is.na(point[1]))
      abline(v = point[1], lty = 2, col = 'darkgray')
}  

mdaplots.scatterg = function(data, pch = 16,  legend = NULL, xlab = NULL,
                             ylab = NULL, show.labels = F, show.lines = F, ...)
{   
   lim = mdaplots.getAxesLim(data, show.lines = show.lines)
   col = mdaplots.getColors(length(data))
   
   if (is.null(xlab))
      xlab = colnames(data[[1]])[1]
   if (is.null(ylab))
      ylab = colnames(data[[1]])[2] 
   
   for (i in 1:length(data))
   {
      if (i == 1)
      {   
         plot(data[[i]][, 1], data[[i]][, 2], type = 'p', 
              col = col[i], pch = pch, 
              xlim = lim$xlim,
              ylim = lim$ylim,
              xlab = xlab,
              ylab = ylab,
              ...)
      }
      else
      {
         points(data[[i]][, 1], data[[i]][, 2], type = 'p', 
                col = col[i], pch = pch)       
      }   
      if (show.labels == T)
         mdaplots.showLabels(data[[i]])
      
      if (!is.null(legend) && length(legend) > 1)
         mdaplots.showLegend(legend, col, pch = pch)
   }
   
   grid()
   
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplots.showLines(show.lines)
   
}

mdaplots.scatter = function(data, pch = 16, cgroup = NULL, show.labels = F, 
                            show.colorbar = T, show.lines = F, ...)
{   
   data = as.matrix(data)
   
   if (!is.null(cgroup))
   {   
      lim = mdaplots.getAxesLim(data, show.colorbar = show.colorbar, show.lines = show.lines)
      col = mdaplots.getColors(cgroup = cgroup)
   }
   else   
   {
      lim = mdaplots.getAxesLim(data, show.lines = show.lines)
      col = mdaplots.getColors(1)
   }
   
   plot(data[, 1], data[, 2], type = 'p', 
        col = col, pch = pch, 
        xlab = colnames(data)[1],
        ylab = colnames(data)[2], 
        xlim = lim$xlim,
        ylim = lim$ylim,
        ...)
   
   if (show.labels == T)
      mdaplots.showLabels(data)
   
   grid()
   
   if (is.numeric(show.lines) && length(show.lines) == 2 )
      mdaplots.showLines(show.lines)
   
   if (!is.null(cgroup) && show.colorbar == T)
      mdaplots.showColorbar(col, cgroup)   
}

mdaplots.line = function(data, type = 'l', show.legend = F, pch = 16, 
                         xlab = NULL, ylab = NULL, show.labels = F, ...)
{
   data = as.matrix(data)
   ny = ncol(data) - 1
   col = mdaplots.getColors(ny)
   lim = mdaplots.getAxesLim(data, multi.y = T)

   if (!(type == 'l' | type == 'b'))
      type = 'l'
   
   if (is.null(xlab))
      xlab = colnames(data)[1]

   if (is.null(ylab))
      ylab = colnames(data)[2] 
   
   for (i in 1:ny)
   { 
      if (i == 1)
      {   
         if (nrow(data) <= 10)
            plot(data[, 1], data[, i + 1], col = col[i], type = type, 
                 xlab = xlab,
                 ylab = ylab, 
                 ylim = lim$ylim,
                 xaxt = 'n',
                 pch = pch,
                 ...)
         else
            plot(data[, 1], data[, i + 1], col = col[i], type = type, 
                 xlab = xlab,
                 ylab = ylab, 
                 ylim = lim$ylim,
                 pch = pch,
                 ...)
         
      }
      else
      {
         lines(data[, 1], data[, i + 1], col = col[i], type = type, pch = pch)         
      }         
   }   
   
   if (nrow(data) < 10)
      axis(side = 1, at = seq(1, nrow(data)), rownames(data))
   
   if (show.legend == T && ny > 1)
      mdaplots.showLegend(colnames(data[, 2:ncol(data)]), col, lty = 1)      
   grid()      
}   

mdaplots.lineg = function(data, pch = 16, type = 'l', show.labels = F,
                          legend = NULL, xlab = NULL, ylab = NULL, ...)
{
   lim = mdaplots.getAxesLim(data)
   col = mdaplots.getColors(length(data))
   
   if (is.null(xlab))
      xlab = colnames(data[[1]])[1]
   if (is.null(ylab))
      ylab = colnames(data[[1]])[2] 
   
   for (i in 1:length(data))
   {
      if (i == 1)
      {   
         plot(data[[i]][, 1], data[[i]][, i + 1], col = col[i], type = type, 
              xlab = xlab,
              ylab = ylab, 
              ylim = lim$ylim,
              pch = pch,
              ...)
      }
      else
      {
         lines(data[[i]][, 1], data[[i]][, 2], col = col[i], type = type, pch = pch)       
      }   
      
      if (show.labels == T)
         mdaplots.showLabels(data[[i]])
      
   }
   
   if (!is.null(legend) && length(legend) > 1)
      mdaplots.showLegend(legend, col, pch = pch)
   grid()   
}   

mdaplots.linescatter = function(data, pch = 16, show.labels = F, ...)
{
   data = as.matrix(data)
   col = mdaplots.getColors(1)
   lim = mdaplots.getAxesLim(data)
   
   plot(data[, 1], data[, 2], type = 'b', 
        col = col, pch = pch, 
        xlab = colnames(data)[1],
        ylab = colnames(data)[2], 
        xlim = lim$xlim,
        ylim = lim$ylim,
        ...)
   if (show.labels == T)
      mdaplots.showLabels(data)
   grid()
}   

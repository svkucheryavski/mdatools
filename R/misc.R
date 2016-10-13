#' Convert image to data matrix
#' 
#' @param img
#' an image (3-way array)
#' 
#' @export
mda.im2data = function(img) {
   width = dim(img)[2]
   height = dim(img)[1]
   nchannels = dim(img)[3]
   
   npixels = width * height
   dim(img) = c(npixels, nchannels)
   attr(img, 'width') = width
   attr(img, 'height') = height
   
   img
}

#' Convert data matrix to an image
#' 
#' @param data
#' data matrix 
#' 
#' @export
mda.data2im = function(data) {
   width = attr(data, 'width', exact = TRUE)
   height = attr(data, 'height', exact = TRUE)
   bgpixels = attr(data, 'bgpixels', exact = TRUE)
  
   if (length(bgpixels) > 0) {
      img = matrix(NA, nrow = nrow(data) + length(bgpixels), ncol = ncol(data))
      img[-bgpixels, ] = data
   } else {
      img = data
   }
      
   dim(img) = c(height, width, ncol(data))
   img
}

#' Remove background pixels from image data
#' 
#' @param data
#' a matrix with image data
#' @param bgpixels
#' vector with indices or logical values corresponding to background pixels
#' 
#' @export
mda.setimbg = function(data, bgpixels) {
   attrs = mda.getattr(data)

   if (length(attrs$exclrows) > 0)
      stop('You can not set background pixels if some of them were excluded!')
   
   # unfold bgpixels to a vector 
   dim(bgpixels) = NULL
   
   # get indices instead of logical values
   if (is.logical(bgpixels))
      bgpixels = which(bgpixels)
   
   # correct indices of bgpixels if some of the pixels were already removed   
   if (length(attrs$bgpixels) > 0) {
      npixels = attrs$width * attrs$height
      row.ind = 1:npixels
      row.ind = row.ind[-attrs$bgpixels]
      bgpixels = row.ind[bgpixels]
   }

   # remove corresponding rows and correct attributes   
   data = data[-bgpixels, ]      
   attrs$bgpixels = unique(c(attrs$bgpixels, bgpixels))

   data = mda.setattr(data, attrs)
   data
}

#' show image data as an image
#' 
#' @param data
#' data with image
#' @param channels
#' indices for one or three columns to show as image channels
#' @param show.excluded
#' logical, if TRUE the method also shows the excluded (hidden) pixels
#' @param main
#' main title for the image
#' @param colmap
#' colormap using to show the intensity levels
#' 
#' @export
imshow = function(data, channels = 1, show.excluded = FALSE, main = NULL, colmap = 'jet') {
   names = colnames(data)
   attrs = mda.getattr(data)
   
   data = mda.subset(data, select = channels)
   data = (data - min(data)) / (max(data) - min(data))
   data = mda.data2im(data)
   
   bg = is.na(data)
   
   nrows = dim(data)[1]
   ncols = dim(data)[2]
   
   if (is.character(colmap) && length(colmap) == 1) {
      if (colmap == 'gray')
         colmap = colorRampPalette(c('#000000', '#ffffff'), space = 'Lab')(256)
      else
         colmap = mdaplot.getColors(256, NULL, colmap)
   }

   if (is.null(main) && !is.null(names))
      main = paste(' ', names[channels], sep = '', collapse = '')
 
   if (length(channels) == 1) {
      nrows = nrow(data)
      image(t(data[seq(nrows, 1, -1), , 1]), xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), 
            main = main, useRaster = T, col = colmap, axes = FALSE)
      if (any(bg)) {
         bgimg = matrix(NA, nrows, ncols)
         bgimg[bg[, , 1]] = 0
         rasterImage(bgimg, 0, 0, 1, 1)
      }
         
   } else {
      if (any(bg))
         data[bg] = 0
      plot(0, main = main, type = 'n', xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', 
           xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
      rasterImage(data, 0, 0, 1, 1)
   } 
  
   # hide excluded pixels with dark gray color
   if (show.excluded == FALSE && length(attrs$exclrows) > 0) {
      ind = 1:(nrows * ncols)
      if (length(attrs$bgpixels) > 0)
         ind = ind[-attrs$bgpixels]
      eximage = rep(NA, nrows * ncols)
      eximage[ind[attrs$exclrows]] = 0.25
      dim(eximage) = c(nrows, ncols)
      rasterImage(eximage, 0, 0, 1, 1)
   }
   
}

#' Wrapper for show() method
#' 
#' @param x
#' data set
#' @param n
#' number of rows to show
#' 
#' @export
mda.show = function(x, n = 50) {
   excl.rows = attr(x, 'exclrows', exact = TRUE)
   excl.cols = attr(x, 'exclcols', exact = TRUE)
  
   name = attr(x, 'name', exact = TRUE)
   xaxis.name = attr(x, 'xaxis.name', exact = TRUE)
   yaxis.name = attr(x, 'yaxis.name', exact = TRUE)
  
   if (!is.null(name) && nchar(name) > 0) {
      cat(sprintf('%s\n%s\n', name, paste(rep('-', nchar(name)), collapse = '')))
   }
   
   if (!is.null(excl.rows))
      x = x[-excl.rows, ]
   
   if (!is.null(excl.cols))
      x = x[, -excl.cols, drop = FALSE]
  
   if (n > nrow(x))
      n = nrow(x)
   
   show(x[1:n, ])   
}

#' A wrapper for subset() method with proper set of attributed
#' 
#' @param x
#' dataset (data frame or matrix)
#' @param subset
#' which rows to keep (indices, names or logical values)
#' @param select
#' which columns to select (indices, names or logical values)
#'  
#' @return 
#' a data with the subset
#'
#' @details 
#' The method works similar to the standard \code{subset()} method, with minor differences. First of all 
#' it keeps (and correct, if necessary) all important attributes. If only columns are selected, it keeps 
#' all excluded rows as excluded. If only rows are selected, it keeps all excluded columns. If both rows
#' and columns are selected it removed all excluded elements first and then makes the subset.
#' 
#' The parameters \code{subset} and \code{select} may each be a vector with numbers or nanes without excluded 
#' elements, or a logical expression.
#' 
#' @export  
mda.subset = function(x, subset = NULL, select = NULL) {
   if (is.null(x))
      return(NULL)
   
   attrs = mda.getattr(x)
   
   if (!is.null(subset)) {
      if (is.logical(subset) & !is.null(attrs$exclrows))
         subset = subset[-attrs$exclrows]
      
      # remove excluded rows first
      if (!is.null(attrs$exclrows))
         x = x[-attrs$exclrows, , drop = F]
      
      # get numeric indices for the rows and subset them
      subset = mda.getexclind(subset, rownames(x), nrow(x))
      x = x[subset, ]
      
      # correct attributes
      if (!is.null(attrs$xaxis.values)) {
         if (!is.null(attrs$exclrows))
            attrs$yaxis.values = attrs$yaxis.values[-attrs$exclrows]
         attrs$yaxis.values = attrs$yaxis.values[subset]
      }
      attrs$exclrows = NULL
   }

   if (!is.null(select)) {
      if (is.logical(select) && !is.null(attrs$exclcols))
         select = select[-attrs$exclcols]
      
      # remove excluded rows first
      if (!is.null(attrs$exclcols))
         x = x[, -attrs$exclcols, drop = F]
      
      # get numeric indices for the rows and subset them
      select = mda.getexclind(select, colnames(x), ncol(x))
      x = x[, select, drop = F]
      
      # correct attributes
      if (!is.null(attrs$xaxis.values)) {
         if (!is.null(attrs$exclcols))
            attrs$xaxis.values = attrs$xaxis.values[-attrs$exclcols]
         attrs$xaxis.values = attrs$xaxis.values[select]
      }
      attrs$exclcols = NULL
   }
  
   x = mda.setattr(x, attrs)
   x
}

# mda.sortrows = function(x, cols = 1) {
#    attrs = mda.getattr(x)
#    
# }
# 
# mda.sortcols = function(x, rows = 1) {
#    attrs = mda.getattr(x)
#    
# }
# 
# mda.nrows = function(x) {
#    
# }
# 
# mda.ncols = function(x) {
#    
# }

#' A wrapper for rbind() method with proper set of attributes
#' 
#' @param ...
#' datasets (data frames or matrices) to bind
#' 
#' @return 
#' the merged datasets
#' 
#' @export
mda.rbind = function(...) {
   objects = list(...)
   nobj = length(objects)
   
   attrs = mda.getattr(objects[[1]])
   out.exclrows = attrs$exclrows
   out.yaxis.values = attrs$yaxis.values
   
   out.x = objects[[1]]
   for (i in 2:nobj) {
      x = objects[[i]]
      exclrows = attr(x, 'exclrows', exact = TRUE)
      yaxis.values = attr(x, 'yaxis.values')
      if (!is.null(exclrows))
         out.exclrows = c(out.exclrows, exclrows + nrow(out.x))
      if (is.null(out.yaxis.values) || is.null(yaxis.values))
         out.yaxis.values = NULL
      else
         out.yaxis.values = c(out.yaxis.values, yaxis.values)
      out.x = rbind(out.x, x)
   }   
   
   out.x = mda.setattr(out.x, attrs)
   attr(out.x, 'exclrows') = out.exclrows
   attr(out.x, 'yaxis.values') = out.yaxis.values
   
   out.x
}

#' A wrapper for cbind() method with proper set of attributes
#' 
#' @param ...
#' datasets (data frames or matrices) to bind
#' 
#' @return 
#' the merged datasets
#' 
#' @export
mda.cbind = function(...) {
   objects = list(...)
   nobj = length(objects)
   
   attrs = mda.getattr(objects[[1]])
   out.exclcols = attrs$exclcols
   out.xaxis.values = attrs$xaxis.values
   
   out.x = objects[[1]]
   for (i in 2:nobj) {
      x = objects[[i]]
      exclcols = attr(x, 'exclcols')
      xaxis.values = attr(x, 'xaxis.values')
      if (!is.null(exclcols))
         out.exclcols = c(out.exclcols, exclcols + ncol(out.x))
      if (is.null(out.xaxis.values) || is.null(xaxis.values))
         out.xaxis.values = NULL
      else
         out.xaxis.values = c(out.xaxis.values, xaxis.values)
      out.x = cbind(out.x, x)
   }   
   
   out.x = mda.setattr(out.x, attrs)
   attr(out.x, 'exclcols') = out.exclcols
   attr(out.x, 'xaxis.values') = out.xaxis.values
   
   out.x
   
}

#' A wrapper for t() method with proper set of attributes
#' 
#' @param x
#' dataset (data frames or matrices) to transpose
#' 
#' @return 
#' the transposed dataset
#' 
#' @export
mda.t = function(x) {
   attrs = mda.getattr(x)
   out.attrs = attrs
   out.attrs$exclrows = attrs$exclcols
   out.attrs$exclcols = attrs$exclrows
   out.attrs$xaxis.name = attrs$yaxis.name
   out.attrs$yaxis.name = attrs$xaxis.name
   out.attrs$xaxis.values = attrs$yaxis.values 
   out.attrs$yaxis.values = attrs$xaxis.values 

   x = t(x)
   x = mda.setattr(x, out.attrs) 
}

#' Exclude/hide rows in a dataset
#' 
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' indices of rows to exclude (numbers, names or logical values)
#' 
#' @return 
#' dataset with excluded rows
#' 
#' @details 
#' The method assign attribute 'exclrows', which contains number of rows, which should be excluded/hidden
#' from calculations and plots (without removing them physically). The argument \code{ind} should contain 
#' rows numbers (excluding already hidden), names or logical values.
#'  
#' @export 
mda.exclrows = function(x, ind) {

   if(is.null(ind))
      return(x)
   
   excl.rows = attr(x, 'exclrows', exact = TRUE)
   nrows.tot = nrow(x)
   nrows.excl = length(excl.rows)
   
   if (nrows.excl == 0) {
      # no objects are excluded yet
      attr(x, 'exclrows') = mda.getexclind(ind, rownames(x), nrows.tot)
   } else {
      # some objects were excluded before
      if (is.logical(ind)) 
         ind = ind[-excl.rows]
      ind = mda.getexclind(ind, rownames(x)[-excl.rows], nrows.tot - nrows.excl)
      ind.tot = 1:nrows.tot
      ind.tot = ind.tot[-excl.rows]
      attr(x, 'exclrows') = sort(unique(c(ind.tot[ind], excl.rows)))
   }
   
   x
}

#' include/unhide the excluded rows
#' 
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' number of excluded rows to include
#' 
#' @return 
#' dataset with included rows
#' 
#' @description 
#' include rows specified by user (earlier excluded using mda.exclrows)
#' 
#' @export 
mda.inclrows = function(x, ind) {
   excl.rows = attr(x, 'exclrows', exact = TRUE)
   ind.log = excl.rows %in% ind
   attr(x, 'exclrows') = excl.rows[!ind.log]
   
   x
}

#' Exclude/hide columns in a dataset
#' 
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' indices of columns to exclude (numbers, names or logical values)
#' 
#' @return 
#' dataset with excluded columns
#' 
#' @details 
#' The method assign attribute 'exclcols', which contains number of columns, which should be excluded/hidden
#' from calculations and plots (without removing them physically). The argument \code{ind} should contain 
#' column numbers (excluding already hidden), names or logical values.
#'  
#' @export 
mda.exclcols = function(x, ind) {
   if(is.null(ind))
      return(x)

   excl.cols = attr(x, 'exclcols', exact = TRUE)
   ncols.tot = ncol(x)
   ncols.excl = length(excl.cols)
   
   if (ncols.excl == 0) {
      # no objects are excluded yet
      attr(x, 'exclcols') = mda.getexclind(ind, colnames(x), ncols.tot)
   } else {
      # some objects were excluded before
      if (is.logical(ind)) 
         ind = ind[-excl.cols]
      ind = mda.getexclind(ind, colnames(x)[-excl.cols], ncols.tot - ncols.excl)
      ind.tot = 1:ncols.tot
      ind.tot = ind.tot[-excl.cols]
      attr(x, 'exclcols') = sort(unique(c(ind.tot[ind], excl.cols)))
   }
   
   x
}

#' Include/unhide the excluded columns
#' 
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' number of excluded columns to include
#' 
#' @return 
#' dataset with included columns.
#' 
#' @description 
#' include colmns specified by user (earlier excluded using mda.exclcols) 
#' 
#' @export 
mda.inclcols = function(x, ind) {
   excl.cols = attr(x, 'exclcols', exact = TRUE)
   ind.log = excl.cols %in% ind
   attr(x, 'exclcols') = excl.cols[!ind.log]
   
   x
}

#' Set data attributes
#' 
#' @description 
#' Set most important data attributes (name, xvalues, excluded rows and columns, etc.) to a dataset
#' 
#' @param x
#' a dataset
#' @param attrs
#' list with attributes
#' @param type
#' a text variable telling which attributes to set ('all', 'row', 'col')
#' 
#' @export
mda.setattr = function(x, attrs, type = 'all') {
   attr(x, 'name') = attrs$name
   attr(x, 'width') = attrs$width
   attr(x, 'height') = attrs$height
   attr(x, 'bgpixels') = attrs$bgpixels
   
   if (type == 'row' || type == 'all') {
      attr(x, 'yaxis.name') = attrs$yaxis.name
      attr(x, 'yaxis.values') = attrs$yaxis.values
      attr(x, 'exclrows') = attrs$exclrows
   }
   
   if (type == 'col' || type == 'all') {
      attr(x, 'xaxis.name') = attrs$xaxis.name
      attr(x, 'xaxis.values') = attrs$xaxis.values
      attr(x, 'exclcols') = attrs$exclcols
   }
   
   x
}

#' 
#' Get data attributes
#' 
#' @description 
#' Returns a list with important data attributes (name, xvalues, excluded rows and columns, etc.)
#' 
#' @param x
#' a dataset
#' 
#' @export
mda.getattr = function(x) {
   attrs = list()
   
   attrs$name = attr(x, 'name', exact = TRUE) 
   attrs$exclrows = attr(x, 'exclrows', exact = TRUE) 
   attrs$exclcols = attr(x, 'exclcols', exact = TRUE) 
   attrs$xaxis.values = attr(x, 'xaxis.values', exact = TRUE) 
   attrs$yaxis.values = attr(x, 'yaxis.values', exact = TRUE) 
   attrs$xaxis.name = attr(x, 'xaxis.name', exact = TRUE) 
   attrs$yaxis.name = attr(x, 'yaxis.name', exact = TRUE) 
   attrs$width = attr(x, 'width', exact = TRUE)
   attrs$height = attr(x, 'height', exact = TRUE)
   attrs$bgpixels = attr(x, 'bgpixels', exact = TRUE)
   
   attrs
}

#' Get indices of excluded rows or columns
#' 
#' @param excl
#' vector with excluded values (logical, text or numbers)
#' @param names
#' vector with names for rows or columns
#' @param n
#' number of rows or columns
#' 
#' @export
mda.getexclind = function(excl, names, n) {
   nitems = ifelse( is.logical(excl), sum(excl), length(excl))
   
   if (is.character(excl))
      excl = which(names %in% excl)
   if (is.logical(excl))
      excl = which(excl)

   if (length(excl) < nitems)
      stop('At least one index or name is incorrect!')
   
   if (!(is.null(excl) || length(excl) == 0) && (!is.numeric(excl) || min(excl) < 1 || max(excl) > n))
      stop('At least one index or name is incorrect!')
   
#   if (length(excl) >= n)
#      stop('All elements of the dataset were excluded, nothing to use!')
   
   excl 
}

#' Convert data frame to a matrix
#' 
#' @description 
#' The function converts data frame to a numeric matrix. 
#' 
#' @param x
#' a data frame
#' @param full
#' logical, if TRUE number of dummy variables for a factor will be the same as number of levels, 
#' otherwise by one smaller 
#'
#' @details 
#' If one or several columns of the data frame are factors they will be converted to a set of dummy 
#' variables. If any columns/rows were hidden in the data frame they will remain hidden in the matrix. If
#' there are factors among the hidden columns, the corresponding dummy variables will be hidden as well.
#' 
#' All other attributes (names, axis names, etc.) will be inherited.
#'   
#' @return 
#' a numeric matrix 
#' 
#' @export
mda.df2mat = function(x, full = FALSE) {
   attrs = mda.getattr(x)

   if (is.null(x) || is.matrix(x) || is.vector(x))
      return(x)
   
   if (is.factor(x))
      x = data.frame(x)
   
   col.fac = unlist(lapply(x, is.factor))
   col.num = which(!col.fac)
   col.fac = which(col.fac)
   
   dummy = function(i, x, full = FALSE) {
      name = colnames(x)[i]
      x = x[, i]
      names = levels(x)
      
      if (full == TRUE)
         n = nlevels(x)
      else
         n = nlevels(x) - 1
      
      d = matrix(0, nrow = length(x), ncol = n)
      colnames = rep('', n)
      for (k in 1:n){
         d[, k] = as.numeric(x) == k
         colnames[k] = names[k]
      }
      colnames(d) = colnames 
      attr(d, 'cols.info') = c(i, n)    
      
      d
   }
   
   if (is.null(col.fac) || length(col.fac) == 0) {
      # no factors among columns - easy job
      x = as.matrix(x)
      x = mda.setattr(x, attrs)
   } else {
      if (!is.null(attrs$exclcols)) {
         if (is.character(attrs$exclcols))
            attrs$exclcols = which(colnames(x) %in% attrs$exclcols)
         if (is.logical(attrs$exclcols))
            attrs$exclcols = which(attrs$exclcols)
         
         exclcols.fac.ind = which(col.fac %in% attrs$exclcols) # hidden factors
         exclcols.num.ind = which(col.num %in% attrs$exclcols) # hidden numeric columns
      } else {
         exclcols.fac.ind = NULL
         exclcols.num.ind = NULL
      }
      
      # split data to numeric columns and factors
      if (length(col.fac) < ncol(x))
         num.data = as.matrix(x[, -col.fac, drop = FALSE])
      else
         num.data = NULL
      
      fac.data = x[, col.fac, drop = FALSE]
      if (!is.null(exclcols.fac.ind) && length(exclcols.fac.ind) > 0) {
         fac.data.hidden = fac.data[, exclcols.fac.ind, drop = FALSE]
         fac.data = fac.data[, -exclcols.fac.ind, drop = FALSE]
      } else {
         fac.data.hidden = NULL
      }
      
      # convert all non-excluded factors to dummy variables
      fac.data = lapply(1:ncol(fac.data), dummy, x = fac.data, full = full)
      fac.data = do.call(cbind, fac.data)
      
      # convert all excluded factors to numeric values
      if (!is.null(fac.data.hidden)) {
         fac.data.hidden = as.matrix(as.data.frame(lapply(fac.data.hidden, as.numeric)))
         n.incl.col = ncol(num.data) + ncol(fac.data)
         exclcols.fac.ind = (n.incl.col + 1):(n.incl.col + ncol(fac.data.hidden))
      } else {
         exclcols.fac.ind = NULL
      }
      
      # combine the data values and set attributes
      x = cbind(num.data, fac.data, fac.data.hidden)
      
      # correct and set arguments
      attrs$exclcols = c(exclcols.num.ind, exclcols.fac.ind)
      x = mda.setattr(x, attrs)
   }
   
   x
}

#' Get selected components
#' 
#' @description
#' returns number of components depending on a user choice
#' 
#' @param obj
#' an MDA model or result object (e.g. \code{pca}, \code{pls}, \code{simca}, etc)
#' @param ncomp
#' number of components to select, provided by user
#' 
#' @details
#' Depedning on a user choice it returns optimal number of component for the model (if 
#' use did not provide any value) or check the user choice for correctness and returns
#' it back
#'  
getSelectedComponents = function(obj, ncomp = NULL)
{
   if (is.null(ncomp))
   {   
      if (is.null(obj$ncomp.selected))
         ncomp = 1
      else
         ncomp = obj$ncomp.selected
   }   
   
   ncomp
}

#' Get main title
#' 
#' @description
#' returns main title for a plot depending on a user choice
#' 
#' @param main
#' main title of a plot, provided by user
#' @param ncomp
#' number of components to select, provided by user
#' @param default
#' default title for the plot
#' 
#' @details
#' Depedning on a user choice it returns main title for a plot
#'  
getMainTitle = function(main, ncomp, default)
{
   if (is.null(main))
   {  
      if (is.null(ncomp))
         main = default
      else
         main = sprintf('%s (ncomp = %d)', default, ncomp)         
   }   
   
   main
}

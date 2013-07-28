prep = 
   setRefClass('prep',
               fields = list(methods = 'list'),
               methods = list(
                  add = function(name, ...)
                  {
                     p = as.list(match.call(expand.dots = TRUE)[-1])
                     methods <<- c(methods, c(name, p))
                  }
               )
   )


prep.autoscale = function(data, center = T, scale = F)
{   
   # Autoscale (mean center and standardize) data matrix.
   #
   # Arguments:
   #   data: a matrix with data values    
   #   center: a logical value or vector with numbers for centering
   #   scale: a logical value or vector with numbers for weighting
   #
   # Returns:
   #   data: preprocessed data values
   
   # define values for centering
   if (is.logical(center) && center == T )
      center = apply(data, 2, mean)
   else if (is.numeric(center))
      center = center   
   
   # define values for weigting
   if (is.logical(scale) && scale == T)
      scale = apply(data, 2, sd)
   else if(is.numeric(scale))
      scale = scale         
   
   # make autoscaling and attach preprocessing attributes
   data = scale(data, center = center, scale = scale)
   attr(data, 'scaled:center') = NULL
   attr(data, 'scaled:scale') = NULL
   attr(data, 'prep:center') = center
   attr(data, 'prep:scale') = scale
   
   data
}

prep.snv = function(data)
{
   # Makes standard normal variate (SNV) preprocessing.
   #
   # Arguments:
   #   data: a matrix with data values    
   #
   # Returns:
   #   data: preprocessed data values
   
   data = t(scale(t(data), center = T, scale = T))
} 

prep.savgol = function(data, width = 3, porder = 1, dorder = 0)
{
   # Apply Savytzky-Golay filter to the data values.
   #
   # Arguments:
   #   data: a matrix with data values    
   #   width: a width of the filter
   #   porder: a polinomial order
   #   dorder: a derivative order
   #
   # Returns:
   #   data: preprocessed data values
   
   nobj = nrow(data)
   nvar = ncol(data)
   
   pdata = matrix(0, ncol = nvar, nrow = nobj)
   
   for (i in 1:nobj)
   {
      d = data[i, ]
      
      w = (width - 1)/2                        
      f  = pinv(outer(-w:w, 0:porder, FUN = "^"))  
      
      d = convolve(d, rev(f[dorder + 1, ]), type = "o")      
      pdata[i, ] = d[(w + 1) : (length(d) - w)] 
   }  
   
   pdata
}

pinv = function(data)
{
   # Calculates pseudo-inverso of data matrix
   s = svd(data)
   s$v %*% diag(1/s$d) %*% t(s$u)
}
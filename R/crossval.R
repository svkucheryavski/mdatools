#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented 
#' cross-validation
#' 
#' @param nobj
#' number of objects in a dataset
#' @param nseg
#' number of segments to split the data to
#' 
#' @return
#' matrix with object indices for each segment
#' 
crossval = function(nobj, nseg = NULL)
{
   if (is.null(nseg))
   {
      if (nobj < 24) { nseg = nobj}
      else if (nobj >= 24 && nobj < 40) { nseg = 8}
      else if (nobj > 40) { nseg = 4 }
   }   
   else if (nseg == 1)
   {
      nseg = nobj
   }   
   
   seglen = ceiling(nobj / nseg)
   fulllen = seglen * nseg
   
   idx = c(sample(1:nobj), rep(NA, fulllen - nobj))
   idx = matrix(idx, nrow = nseg, byrow = T)   
   
   return (idx)        
}   
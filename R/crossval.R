# function returns indices for cross-validation loop
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
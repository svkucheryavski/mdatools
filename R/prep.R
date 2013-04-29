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
   if (is.logical(center) && center == T )
      center = apply(data, 2, mean)
   else if (is.numeric(center))
      center = center   
   
   if (is.logical(scale) && scale == T)
      scale = apply(data, 2, sd)
   else if(is.numeric(scale))
      scale = scale         
   
   data = scale(data, center = center, scale = scale)
   attr(data, 'scaled:center') = NULL
   attr(data, 'scaled:scale') = NULL
   attr(data, 'prep:center') = center
   attr(data, 'prep:scale') = scale
   
   return (data)
}

prep.snv = function(X)
{
   X = t(X)
   X = scale(X, center = T, scale = T)
   X = t(X)
} 

prep.savgol <- function(TT, fl, forder = 1, idorder = 0)
{
   nobj = nrow(TT)
   TT2 = matrix(0, ncol = ncol(TT), nrow = nrow(TT))
   
   for (i in 1:nobj)
   {
      T = TT[i,]
      m <- length(T)
      dorder = idorder + 1
      
      fc <- (fl - 1)/2                        
      X  <- outer(-fc:fc, 0:forder, FUN="^")  
      
      Y  <- pinv(X);                                
      T2 <- convolve(T, rev(Y[dorder, ]), type="o")
      T2 <- T2[(fc+1):(length(T2)-fc)]
      TT2[i,] = T2
   }  
   TT2
}

pinv <- function(A)
{
   s <- svd(A)
   s$v %*% diag(1/s$d) %*% t(s$u)
}
#' Autoscale values
#' 
#' @description
#' Autoscale (mean center and standardize) values in columns of data matrix.
#' 
#' @param data
#' a matrix with data values
#' @param center
#' a logical value or vector with numbers for centering
#' @param scale
#' a logical value or vector with numbers for weighting
#' 
#' @return
#' data matrix with processed values
#' 
prep.autoscale = function(data, center = T, scale = F)
{   
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

#' Standard Normal Variate transformation
#' 
#' @description
#' Applies Standard Normal Variate (SNV) transformation to the rows of data matrix
#' 
#' @param data
#' a matrix with data values
#' 
#' @return
#' data matrix with processed values
#' 
#' @details
#' SNV is a simple preprocessing to remove scatter effects (baseline offset and slope) from 
#' spectral data, e.g. NIR spectra.
#'  
#'  @examples
#'  
#'  ### Apply SNV to spectra from simdata
#'  
#'  library(mdatools)
#'  data(simdata)
#'  
#'  spectra = simdata$spectra.c
#'  wavelength = simdata$wavelength
#'  
#'  cspectra = prep.snv(spectra)
#'  
#'  par(mfrow = c(2, 1))
#'  mdaplot(cbind(wavelength, t(spectra)), type = 'l', main = 'Before SNV')
#'  mdaplot(cbind(wavelength, t(cspectra)), type = 'l', main = 'After SNV')
#'
prep.snv = function(data)
{
   data = t(scale(t(data), center = T, scale = T))
} 

#' Savytzky-Golay filter
#' 
#' @description
#' Applies Savytzky-Golay filter to the rows of data matrix
#' 
#' @param data
#' a matrix with data values
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
#' 
prep.savgol = function(data, width = 3, porder = 1, dorder = 0)
{
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

#' Multiplicative Scatter Correction transformation
#' 
#' @description
#' Applies Multiplicative Scatter Correction (MSC) transformation to data matrix (spectra)
#' 
#' @param spectra
#' a matrix with spectra values
#' @param mspectrum
#' mean spectrum (if NULL will be calculated from \code{spectra})
#' 
#' @return
#' list with two fields - preprocessed spectra and calculated mean spectrum
#' 
#' @details
#' MSC is used to remove scatter effects (baseline offset and slope) from 
#' spectral data, e.g. NIR spectra.
#'  
#'  @examples
#'  
#'  ### Apply MSC to spectra from simdata
#'  
#'  library(mdatools)
#'  data(simdata)
#'  
#'  spectra = simdata$spectra.c
#'  wavelength = simdata$wavelength
#'  
#'  res = prep.msc(spectra)
#'  cspectra = res$cspectra
#'  
#'  par(mfrow = c(2, 1))
#'  mdaplot(cbind(wavelength, t(spectra)), type = 'l', main = 'Before MSC')
#'  mdaplot(cbind(wavelength, t(cspectra)), type = 'l', main = 'After MSC')
#'
prep.msc = function(spectra, mspectrum = NULL)
{
   if (is.null(mspectrum))
      mspectrum = apply(spectra, 2, mean)   
   
   cspectra = matrix(0, nrow = nrow(spectra), ncol = ncol(spectra))
   for (i in 1:nrow(spectra))
   {
      coef = coef(lm(spectra[i, ] ~ mspectrum))
      cspectra[i, ] = (spectra[i, ] - coef[1]) / coef[2]   
   }
   
   list(cspectra = cspectra, mspectrum = mspectrum)
}  

# prep.alsbasecor = function(spectra, penalty = 0.1, smoothness = 10^5)
# {
#    m = ncol(spectra)
#    
#    if (m == 1)
#    {
#       
#    }  
#       
#    D = diff(speye(m), 2);
#    w = matrix(1, nrow = m, ncol = 1)
#    
#    for (n in 1:nrow(spectra))
#    {  
#       spectrum = spectra[n, ]
#       
#       for (it in 1:20)
#       {   
#          W = spdiags(w, 0, m, m);
#          C = chol(W + smoothness * t(D) * D)
#          baseline = C \ (t(C) \ (w * spectrum));
#          w = penalty * (spectrum > baseline) + (1 - penalty) * (spectrum < baseline);
#       }
#    
#       spectra[n, ] = spectrum - t(baseline) 
#    }                
# }   

#' Pseudo-inverse matrix
#' 
#' @description
#' Computes pseudo-inverse matrix using SVD
#' 
#' @param data
#' a matrix with data values to compute inverse for
#' 
pinv = function(data)
{
   # Calculates pseudo-inverse of data matrix
   s = svd(data)
   s$v %*% diag(1/s$d) %*% t(s$u)
}
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
#' @param max.cov
#' columns that have coefficient of variation (in percent) below `max.cv` will not be scaled
#' 
#' @return
#' data matrix with processed values
#' 
#' @export
prep.autoscale = function(data, center = T, scale = F, max.cov = 0.1) {   
   
   attrs = mda.getattr(data)
   dimnames = dimnames(data)
   
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
   
   if(is.numeric(scale)) {
      if (!is.numeric(center))
         m = apply(data, 2, mean)
      else
         m = center
      cv = scale/m * 100
      scale[is.nan(cv) | cv < 0.1] = 1
   } 
   data = scale(data, center = center, scale = scale)
   
   data = mda.setattr(data, attrs)
   attr(data, 'scaled:center') = NULL
   attr(data, 'scaled:scale') = NULL
   attr(data, 'prep:center') = center
   attr(data, 'prep:scale') = scale
   dimnames(data) = dimnames
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
#' @export
prep.snv = function(data){
   attrs = mda.getattr(data)
   dimnames = dimnames(data)
   
   data = t(scale(t(data), center = T, scale = T))
   data = mda.setattr(data, attrs)
   dimnames(data) = dimnames
   data
} 

#' Normalization
#' 
#' @description
#' Normalizes signals (rows of data matrix) to unit area or unit length
#' 
#' @param data
#' a matrix with data values
#' @param type
#' type of normalization \code{'area'} or \code{'length'}
#' 
#' @return 
#' data matrix with normalized values
#' 
#' @export
prep.norm = function(data, type = 'area') {
   attrs = mda.getattr(data)
   dimnames = dimnames(data)
   
   if (type == 'area')
   {   
      w = apply(abs(data), 1, sum)
   }   
   else if (type == 'length')
   {
      w = apply(data^2, 1, sum)
      w = sqrt(w)
   }   
   else
   {   
      stop('Wrong value for argument "type"!')
   }   
   
   data = sweep(data, 1, w, '/')
   data = mda.setattr(data, attrs)
   dimnames(data) = dimnames
   data
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
#' @export
prep.savgol = function(data, width = 3, porder = 1, dorder = 0) {
   attrs = mda.getattr(data)
   dimnames = dimnames(data)
   
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
   
   pdata = mda.setattr(pdata, attrs)
   dimnames(pdata) = dimnames
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
#' @export
prep.msc = function(spectra, mspectrum = NULL) {
   attrs = mda.getattr(spectra)
   dimnames = dimnames(spectra)
   
   if (is.null(mspectrum))
      mspectrum = apply(spectra, 2, mean)   
   
   cspectra = matrix(0, nrow = nrow(spectra), ncol = ncol(spectra))
   for (i in 1:nrow(spectra))
   {
      coef = coef(lm(spectra[i, ] ~ mspectrum))
      cspectra[i, ] = (spectra[i, ] - coef[1]) / coef[2]   
   }
   
   cspectra = mda.setattr(cspectra, attrs)
   attr(cspectra, 'mspectrum') = mspectrum
   dimnames(cspectra) = dimnames
   cspectra
}  

#' Pseudo-inverse matrix
#' 
#' @description
#' Computes pseudo-inverse matrix using SVD
#' 
#' @param data
#' a matrix with data values to compute inverse for
#' 
#' @export
pinv = function(data)
{
   # Calculates pseudo-inverse of data matrix
   s = svd(data)
   s$v %*% diag(1/s$d) %*% t(s$u)
}


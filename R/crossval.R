#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented 
#' cross-validation
#' 
#' @param nobj
#' number of objects in a dataset
#' @param cv
#' cross-validation settings, can be a number or a list. If cv is a number, it will be 
#' used as a number of segments for random cross-validation (if cv = 1, full cross-validation 
#' will be preformed), if it is a list, the following syntax can be used: cv = list('rand', nseg, nrep) 
#' for random repeated cross-validation with nseg segments and nrep repetitions or cv = list('ven', nseg) 
#' for systematic splits to nseg segments ('venetian blinds').  
#' 
#' @return
#' matrix with object indices for each segment
#' 
#' @export
crossval = function(nobj, cv = NULL) {
   methods = c('rand', 'ven', 'loo')
   
   if (is.null(cv))
   {
      type = 'rand'
      rep = 1
      if (nobj < 24) { nseg = nobj}
      else if (nobj >= 24 && nobj < 40) { nseg = 8}
      else if (nobj >= 40) { nseg = 4 }
   }   
   else if (is.numeric(cv))
   {
      type = 'rand'
      nrep = 1
      if (cv == 1)
         nseg = nobj
      else
         nseg = cv
   }
   else
   {
      type = cv[[1]]
      
      if ( type == 'loo' )
      {
         type = 'rand'
         nseg = nobj
         nrep = 1
      }   
      else
      {   
         nseg = cv[[2]]
      
         if (length(cv) > 2)
            nrep = cv[[3]]
         else
            nrep = 1
      }   
   }   
   
   if ( !(type %in% methods) ) 
      stop('Wrong name for cross-valudation method!')      
   
   if (type != 'rand')
      nrep = 1
   
   seglen = ceiling(nobj / nseg)
   fulllen = seglen * nseg
   
   idx = array(0, dim = c(nseg, ceiling(nobj / nseg), nrep))
      
   if (type == 'rand')
   {   
      for (irep in 1:nrep)
      {   
         v = c(sample(1:nobj), rep(NA, fulllen - nobj))
         idx[, , irep] = matrix(v, nrow = nseg, byrow = T)   
      }
   }
   else if (type == 'ven')
   {
      v = c(1:nobj, rep(NA, fulllen - nobj))
      idx[, , 1] = matrix(v, nrow = nseg, byrow = F)               
   }   
   return (idx)        
}   

#' String with description of cross-validation method
#'
#' @param cv 
#' a list with cross-validation settings
#' 
#' @return
#' a string with the description text
#'
crossval.str = function(cv)
{
   str = '';
   if (is.numeric(cv))
   {
      if (cv == 1)
         str = 'full cross-validation'
      else
         str = sprintf('random cross-validation with %d segments', cv)
   }
   else
   {
      type = cv[[1]]
      if ( type == 'loo' )
         str = 'full cross-validation'
      else if ( type == 'ven')
      {   
         str = sprintf('venetian blinds with %d segments', cv[[2]])
      }   
      else
      {
         str = sprintf('random cross-validation with %d segments', cv[[2]])
      }   
   }   
   
   str
}   
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

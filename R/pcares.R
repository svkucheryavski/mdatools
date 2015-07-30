#' Results of PCA decomposition
#' @description 
#' \code{pcares} is used to store results for PCA decomposition of data.
#'
#'
#' @param scores
#' matrix with score values (nobj x ncomp).
#' @param loadings
#' matrix with loading values (nvar x ncomp).
#' @param residuals
#' matrix with data residuals.
#' @param totvar
#' total variance for data matrix.
#' @param tnorm
#' vector with singular values for scores normalization.
#' @param ncomp.selected
#' selected number of components for decomposition.
#' @param ...
#' other arguments supported by \code{ldecomp}.
#' 
#' @details 
#' In fact \code{pcares} is a wrapper for \code{\link{ldecomp}} - general class for storing
#' results for linear decomposition X = TP' + E. So, most of the methods, arguments and 
#' returned values are inherited from \code{ldecomp}.
#'   
#' There is no need to create a \code{pcares} object manually, it is created automatically when 
#' build a PCA model (see \code{\link{pca}}) or apply the model to a new data (see 
#' \code{\link{predict.pca}}). The object can be used to show summary and plots for the results.
#' 
#' @return 
#' Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q }{matrix with Q residuals (nobj x ncomp).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#' 
#' @seealso 
#' Methods for \code{pcares} objects:
#' \tabular{ll}{
#'  \code{print.pcares} \tab shows information about the object.\cr
#'  \code{summary.pcares} \tab shows statistics for the PCA results.\cr
#' }
#' 
#' Methods, inherited from \code{\link{ldecomp}} class:   
#' \tabular{ll}{
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#'  \code{\link{plotResiduals.ldecomp}} \tab makes Q vs. T2 residuals plot.\cr
#' }
#' 
#' Check also \code{\link{pca}} and \code{\link{ldecomp}}.    
#'
#' @examples
#' ### Examples for PCA results class
#' 
#' library(mdatools)
#' 
#' ## 1. Make a model for every odd row of People data
#' ## and apply it to the objects from every even row
#' 
#' data(people)
#' x = people[seq(1, 32, 2), ]
#' x.new = people[seq(1, 32, 2), ]
#' 
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#' 
#' res = predict(model, x.new)
#' summary(res)
#' plot(res)
#' 
#' ## 1. Make PCA model for People data with autoscaling
#' ## and full cross-validation and get calibration results
#' 
#' 
#' data(people)
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#' 
#' res = model$calres
#' summary(res)
#' plot(res)
#' 
#' ## 2. Show scores plots for the results
#' par(mfrow = c(2, 2))
#' plotScores(res)
#' plotScores(res, cgroup = people[, 'Beer'], show.labels = TRUE)
#' plotScores(res, comp = c(1, 3), show.labels = TRUE)
#' plotScores(res, comp = 2, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 3. Show residuals and variance plots for the results
#' par(mfrow = c(2, 2))
#' plotVariance(res, type = 'h')
#' plotCumVariance(res, show.labels = TRUE, legend.position = 'bottomright')
#' plotResiduals(res, show.labels = TRUE, cgroup = people[, 'Sex'])
#' plotResiduals(res, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export 
pcares = function(scores, loadings, residuals, totvar, tnorm = NULL, ncomp.selected = NULL, ...)
{
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.
   
   res = ldecomp(scores, loadings, residuals, totvar, tnorm = tnorm, ncomp.selected = ncomp.selected
                 , ...)
   class(res) = c('pcares', 'ldecomp')   
   
   res
}   


#' Plot method for PCA results object
#' 
#' @description
#' Show several plots to give an overview about the PCA results
#' 
#' @param x
#' PCA results (object of class \code{pcares})
#' @param comp
#' which components to show the scores plot for (can be one value or vector with two values).
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @export 
plot.pcares = function(x, comp = c(1, 2), show.labels = T, ...)
{   
   par(mfrow = c(2, 2))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, show.labels = show.labels, ...)
   plotVariance(x, show.labels = show.labels, ...)
   plotCumVariance(x, show.labels = show.labels, ...)
   par(mfrow = c(1, 1))
}


#' Summary method for PCA results object
#' 
#' @description
#' Shows some statistics (explained variance, eigenvalues) about the results.
#' 
#' @param object
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#' 
#' @export
summary.pcares = function(object, ...)
{
   summary.ldecomp(object, 'Summary for PCA results', ...)
}


#' Print method for PCA results object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#'
#' @export 
print.pcares = function(x, ...)
{   
   print.ldecomp(x, 'Results for PCA decomposition (class pcares)', ...)
   cat('\n')
}
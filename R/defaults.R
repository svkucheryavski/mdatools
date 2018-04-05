#' Confusion matrix for classification results
#' 
#' @details 
#' Returns confusion matrix for classification results represented by the object.
#' 
#' @param obj
#' classification results (object of class \code{simcares}, \code{simcamres}, etc)
#' @param ...
#' other parameters.
#' 
#' @export
getConfusionMatrix = function(obj, ...) {
   UseMethod("getConfusionMatrix")
}  

#' Plot for class belonging probability
#' 
#' @description
#' Makes a plot with class belonging probabilities for each object of the classification results. 
#' Works only with classification methods, which compute this probability (e.g. SIMCA). 
#' 
#' @param obj
#' an object with classification results (e.g. SIMCA)
#' @param ...
#' other parameters
#' 
#' @export
plotProbabilities = function(obj, ...) {
   UseMethod("plotProbabilities")
}  

#' Get class belonging probability
#' 
#' @description
#' Compute class belonging probabilities for classification results.
#' 
#' @param obj
#' an object with classification results (e.g. SIMCA)
#' @param ...
#' other parameters
#' 
#' @export
getProbabilities = function(obj, ...) {
   UseMethod("getProbabilities")
}  

#' Set residual limits for PCA model
#' 
#' @description
#' Calculates and set critical limits for residuals of PCA model
#' 
#' @param obj
#' a SIMCA model
#' @param ...
#' other parameters
#' 
#' @export
setResLimits = function(obj, ...) {
   UseMethod("setResLimits")
}  


#' Shows extreme plot for SIMCA model
#' 
#' @description
#' Generic function for creating extreme plot for SIMCA model
#' 
#' @param obj
#' a SIMCA model
#' @param ...
#' other parameters
#' 
#' @export
plotExtreme = function(obj, ...) {
   UseMethod("plotExtreme")
}  

#' Get regression coefficients
#' 
#' @description
#' Generic function for getting regression coefficients from PLS model
#' 
#' @param obj
#' a PLS model
#' @param ...
#' other parameters
#' 
#' @export
getRegcoeffs = function(obj, ...) {
   UseMethod("getRegcoeffs")
}  

#' VIP scores plot
#' 
#' @description
#' Generic function for plotting VIP scores values for regression model (PCR, PLS, etc)
#' 
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#' 
#' @export
plotVIPScores = function(obj, ...) {
   UseMethod("plotVIPScores")
}  

#' VIP scores 
#' 
#' @description
#' Generic function for returning VIP scores values for regression model (PCR, PLS, etc)
#' 
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#' 
#' @export
getVIPScores = function(obj, ...) {
   UseMethod("getVIPScores")
}  

#' Selectivity ratio plot
#' 
#' @description
#' Generic function for plotting selectivity ratio values for regression model (PCR, PLS, etc)
#' 
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#' 
#' @export
plotSelectivityRatio = function(obj, ...) {
   UseMethod("plotSelectivityRatio")
}  

#' Selectivity ratio 
#' 
#' @description
#' Generic function for returning selectivity ratio values for regression model (PCR, PLS, etc)
#' 
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#' 
#' @export
getSelectivityRatio = function(obj, ...) {
   UseMethod("getSelectivityRatio")
}  

#' Select optimal number of components for a model
#' 
#' @description
#' Generic function for selecting number of components for multivariate models (e.g. PCA, PLS, ...)
#' 
#' @param model
#' a model object
#' @param ncomp
#' number of components to select
#' 
#' @export
selectCompNum = function(model, ncomp) {   
   UseMethod("selectCompNum")
}   

#' Cooman's plot
#' 
#' @details
#' Generic function for Cooman's plot
#' 
#' @param obj
#' classification model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotCooman = function(obj, ...)
{
   UseMethod("plotCooman")
}

#' Model distance plot
#' 
#' @details
#' Generic function for plotting distance from object to a multivariate model
#' 
#' @param obj
#' a model object
#' @param ...
#' other arguments
#' 
#' @export
plotModelDistance = function(obj, ...)
{
   UseMethod("plotModelDistance")
}

#' Discrimination power plot
#' 
#' @details
#' Generic function for plotting discrimination power values for classification model
#' 
#' @param obj
#' a model object
#' @param ...
#' other arguments
#' 
#' @export
plotDiscriminationPower = function(obj, ...)
{
   UseMethod("plotDiscriminationPower")      
}  

#' Calibration data
#' 
#' @details
#' Generic function getting calibration data from a linear decomposition model (e.g. PCA)
#' 
#' @param obj
#' a model object
#' @param ...
#' other arguments
#' 
#' @export
getCalibrationData = function(obj, ...)
{
   UseMethod("getCalibrationData")   
}  

#' Modelling power plot
#' 
#' @details
#' Generic function for plotting modelling power values for classification model
#' 
#' @param obj
#' a model object
#' @param ...
#' other arguments
#' 
#' @export
plotModellingPower = function(obj, ...)
{
   UseMethod("plotModellingPower")   
}  

#' Misclassification ratio plot
#' 
#' @details
#' Generic function for plotting missclassification values for classification model or results
#' 
#' @param obj
#' a model or a result object
#' @param ...
#' other arguments
#' 
#' @export
plotMisclassified = function(obj, ...)
{
   UseMethod('plotMisclassified')
}

#' Specificity plot
#' 
#' @details
#' Generic function for plotting specificity values for classification model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotSpecificity = function(obj, ...)
{
   UseMethod("plotSpecificity")   
}  

#' Sensitivity plot
#' 
#' @details
#' Generic function for plotting sensitivity values for classification model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotSensitivity = function(obj, ...)
{
   UseMethod("plotSensitivity")   
}

#' Classification performance plot
#' 
#' @details
#' Generic function for plotting classification performance for model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotPerformance = function(obj, ...)
{
   UseMethod("plotPerformance")   
}  

#' Predictions
#' 
#' @details
#' Generic function for showing predicted values for classification or regression model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
showPredictions = function(obj, ...)
{
   UseMethod("showPredictions")   
}  

#' X residuals plot
#' 
#' @details
#' Generic function for plotting x residuals for classification or regression model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXResiduals = function(obj, ...)
{   
   UseMethod("plotXResiduals")
}   

#' Y residuals plot
#' 
#' @details
#' Generic function for plotting y residuals for classification or regression model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotYResiduals = function(obj, ...)
{   
   UseMethod("plotYResiduals")
}   

#' X variance plot
#' 
#' @details
#' Generic function for plotting explained variance for decomposition of x data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXVariance = function(obj, ...)
{   
   UseMethod("plotXVariance")
}   

#' Y variance plot
#' 
#' @details
#' Generic function for plotting explained variance for decomposition of y data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotYVariance = function(obj, ...)
{   
   UseMethod("plotYVariance")
}   

#' Biplot
#' 
#' @details
#' Generic function for biplot
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotBiplot = function(obj, ...)
{   
   UseMethod("plotBiplot")
}   

#' Scores plot
#' 
#' @details
#' Generic function for scores values for data decomposition
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotScores = function(obj, ...)
{   
   UseMethod("plotScores")
}   

#' X scores plot
#' 
#' @details
#' Generic function for plotting scores values for decomposition of x data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXScores = function(obj, ...)
{   
   UseMethod("plotXScores")
}   

#' XY scores plot
#' 
#' @details
#' Generic function for plotting scores values for decomposition of x and y data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXYScores = function(obj, ...)
{   
   UseMethod("plotXYScores")
}   

#' Selected intervals plot
#' 
#' @details
#' Generic function for plotting selected intervals or variables
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotSelection = function(obj, ...)
{   
   UseMethod("plotSelection")
}   

#' RMSE plot
#' 
#' @details
#' Generic function for plotting RMSE values vs. complexity of a regression model
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotRMSE = function(obj, ...)
{   
   UseMethod("plotRMSE")
}   

#' Variance plot
#' 
#' @details
#' Generic function for plotting explained variance for data decomposition 
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotCumVariance = function(obj, ...)
{   
   UseMethod("plotCumVariance")
}   

#' X cumulative variance plot
#' 
#' @details
#' Generic function for plotting cumulative explained variance for decomposition of x data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXCumVariance = function(obj, ...)
{   
   UseMethod("plotXCumVariance")
}   

#' Y cumulative variance plot
#' 
#' @details
#' Generic function for plotting cumulative explained variance for decomposition of y data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotYCumVariance = function(obj, ...)
{   
   UseMethod("plotYCumVariance")
}   

#' Loadings plot
#' 
#' @details
#' Generic function for plotting loadings values for data decomposition
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotLoadings = function(obj, ...)
{   
   UseMethod("plotLoadings")
}   

#' Predictions plot
#' 
#' @details
#' Generic function for plotting predicted values for classification or regression model or results
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotPredictions = function(obj, ...)
{   
   UseMethod("plotPredictions")
}   

#' Regression coefficients plot
#' 
#' @details
#' Generic function for plotting regression coefficients values for a regression model
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotRegcoeffs = function(obj, ...)
{   
   UseMethod("plotRegcoeffs")
}   

#' Residuals plot
#' 
#' @details
#' Generic function for plotting residual values for data decomposition
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotResiduals = function(obj, ...)
{   
   UseMethod("plotResiduals")
}   

#' Variance plot
#' 
#' @details
#' Generic function for plotting explained variance for data decomposition
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotVariance = function(obj, ...)
{   
   UseMethod("plotVariance")
}   

#' X loadings plot
#' 
#' @details
#' Generic function for plotting loadings values for decomposition of x data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXLoadings = function(obj, ...)
{   
   UseMethod("plotXLoadings")
}   

#' X loadings plot
#' 
#' @details
#' Generic function for plotting loadings values for decomposition of x and y data
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotXYLoadings = function(obj, ...)
{   
   UseMethod("plotXYLoadings")
}   

#' Statistic histogram 
#' 
#' @details
#' Generic function for plotting statistic histogram plot
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotHist = function(obj, ...)
{   
   UseMethod("plotHist")
}   

#' Correlation plot 
#' 
#' @details
#' Generic function for correlation plot
#' 
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#' 
#' @export
plotCorr = function(obj, ...)
{   
   UseMethod("plotCorr")
}   

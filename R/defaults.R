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
plotVIPScores = function(obj, ...)
{
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
getVIPScores = function(obj, ...)
{
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
plotSelectivityRatio = function(obj, ...)
{
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
getSelectivityRatio = function(obj, ...)
{
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
selectCompNum = function(model, ncomp)
{   
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
plotYVariance = function(obj, ...)
{   
   UseMethod("plotYVariance")
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
plotXYScores = function(obj, ...)
{   
   UseMethod("plotXYScores")
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
plotCorr = function(obj, ...)
{   
   UseMethod("plotCorr")
}   

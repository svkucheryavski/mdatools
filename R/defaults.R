# default function for interclass methods #

getSelectedComponents = function(obj, ncomp = NULL)
{
   UseMethod("getSelectedComponents")
}  

plotCooman = function(object, ...)
{
   UseMethod("plotCooman")
}

plotModelDistance = function(object, ...)
{
   UseMethod("plotModelDistance")
}

plotDiscriminationPower = function(object, ...)
{
   UseMethod("plotDiscriminationPower")      
}  

getCalibrationData = function(object, ...)
{
   UseMethod("getCalibrationData")   
}  

plotModellingPower = function(object, ...)
{
   UseMethod("plotModellingPower")   
}  

plotMisclassified = function(object, ...)
{
   UseMethod('plotMisclassified')
}

plotSpecificity = function(object, ...)
{
   UseMethod("plotSpecificity")   
}  

plotSensitivity = function(object, ...)
{
   UseMethod("plotSensitivity")   
}

plotPerformance = function(object, ...)
{
   UseMethod("plotPerformance")   
}  

showPredictions = function(object, ...)
{
   UseMethod("showPredictions")   
}  

selectCompNum = function(object, ...)
{   
   UseMethod("selectCompNum")
}   

plotXResiduals = function(object, ...)
{   
   UseMethod("plotXResiduals")
}   

plotYResiduals = function(object, ...)
{   
   UseMethod("plotYResiduals")
}   

plotXVariance = function(object, ...)
{   
   UseMethod("plotXVariance")
}   

plotYVariance = function(object, ...)
{   
   UseMethod("plotYVariance")
}   

plotScores = function(object, ...)
{   
   UseMethod("plotScores")
}   

plotXScores = function(object, ...)
{   
   UseMethod("plotXScores")
}   

plotXYScores = function(object, ...)
{   
   UseMethod("plotXYScores")
}   

plotRMSE = function(object, ...)
{   
   UseMethod("plotRMSE")
}   

plotCumVariance = function(object, ...)
{   
   UseMethod("plotCumVariance")
}   

plotXCumVariance = function(object, ...)
{   
   UseMethod("plotXCumVariance")
}   

plotYCumVariance = function(object, ...)
{   
   UseMethod("plotYCumVariance")
}   

plotLoadings = function(object, ...)
{   
   UseMethod("plotLoadings")
}   

plotPredictions = function(object, ...)
{   
   UseMethod("plotPredictions")
}   

plotRegcoeffs = function(object, ...)
{   
   UseMethod("plotRegcoeffs")
}   

plotResiduals = function(object, ...)
{   
   UseMethod("plotResiduals")
}   

plotVariance = function(object, ...)
{   
   UseMethod("plotVariance")
}   

plotXLoadings = function(object, ...)
{   
   UseMethod("plotXLoadings")
}   

plotXYLoadings = function(object, ...)
{   
   UseMethod("plotXYLoadings")
}   


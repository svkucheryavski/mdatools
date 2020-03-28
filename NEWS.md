v.0.10.3
========

* Fixed bug `#85` when using y-values as data frame gave an error in PLS regression

* Fixed bug `#86` and changed the way PLS limits maximum number of components to avoid problems with singular matrices. Now if PLS algorithm finds during calculations that provided number of components is too large, it gives a warning and reduces this number.

* Code refactoring and tests for preprocessing methods

v.0.10.2
========

* Fixed a bug in `categorize.pls()` method, which could give wrong results for test set measurements (see issue #82).

v.0.10.1
========

* Small improvements to  `plotExtreme.pca()` so user can specify additional parameters, such as, for example `cex`. If plot is made for several components, you can now specify just one value for all points (e.g. color of points or marker symbol).

* Parameter `show.limits` in methods `plotResiduals.pca()`, `plotXResiduals.pls()`, `plotXYResiduals.pls()` can now take two logical values — first for extreme limit and second for outlier limit. So, you can show only one of the two limits on the plot. If one value is specified it will be taken for both limits.

* New function `plotHotellingEllipse()` adds Hotelling T^2^ ellipse to any scatter plot (of course it is made first of all for PCA and PLS scores plots). The function works similar to `plotConvexHull()` and `plotConfidenceEllipse()`, see help for examples.

* Fixed a bug in `summary()` method for PLS, which worked incorrectly in case of several response variables (PLS2).

v.0.10.0
========

Many changes have been made in this version, but most of them are under the hood. Code has been refactored significantly in order to improve its efficiency and make future support easier. Some functionality has been re-written from the scratch. **Most** of the code is backward compatible, which means your old scripts should have no problem to run with this version. However, some changes are incompatible and this can lead to occasional errors and warning messages. All details are shown below, pay a special attention to **breaking changes** part.

Another important thing is the way cross-validation works starting from this version. It was decided to use cross-validation only for computing performance statistics, e.g. error of predictions in PLS or classification error in SIMCA or PLS-DA. Decomposition results, such as explained variance or residual distances are not computed for cross-validation anymore. It was a bad idea from the beginning, as the way it has been implemented is not fully correct — distances and variances measured for different local models should not be compared directly. After a long consideration it was decided to implement this part in a more correct and conservative way.

Finally, all model results (calibration, cross-validation and test set validation), are now combined
into a single list, `model$res`. This makes a lot of things easier. However, the old way of
accessing the result objects (e.g. `model$calres` or `model$cvres`) still works, you can access e.g. calibration results both using `model$res$cal` and `model$calres`, so this change will not break the compatibility.

Below is more detailed list of changes. The [tutorial](http://mdatools.com/docs/) has been updated accordingly.

## Breaking changes

Here are changes which can potentially lead to error messages in previously written code.

* Cross-validation results are no more available for PCA (as mentioned above), so any use of `model$cvres` object for PCA model will lead to an error. For the same reason `pca()` does not take the `cv` parameter anymore.

* Method `plotModellingPower()` is no longer available (was used for SIMCA models).

* Method `plotResiduals()` is no longer available for SIMCAM models (multiclass SIMCA), use
corresponding method for individual models instead.

* Selectivity ratio and VIP scores are not a part of PLS model anymore. This is done to make the calibration of models faster. Use `selratio()` and `vipscores()` to compute them. Functions `plotSelectivityRatio()` and `plotVIPScores()` are still available but they both compute the values first, which may take a bit of time on large datasets. This change makes parameter `light` superfluous and it is no more supported in `pls()`.

* Other two parameters, which are no more needed when you use `pls()`, are `coeffs.ci` and `coeffs.alpha`. Jack-Knifing based confidence intervals for regression coefficients now automatically computed every time you use cross-validation. You can specify the significance level for the intervals when you either visualize them using `plot.regcoeffs()` or `plotRegcoeffs()` for PLS model or when you get the values by using `getRegcoeffs()`.

* When you make prediction plot for any classification model, you should specify name of result
object to show the predictions for. In old versions the name of results were `"calres"`, `"cvres"`,
`"testres"`. From this version they have been changed to `"cal"`, `"cv"` and `"test"`
correspondingly.

* In PLS-DA there was a possibility to show predictions not for classification results but for regression model the PLS-DA is built upon using the following code: `plotPredictions(structure(model, class = "pls"))`. From this version you should use `plotPredictions(structure(model, class = "regmodel"))` instead, as the `plotPredictions()` function for regression has been moved from `pls` class to its parent, more general class, `regmodel`.

* In methods `plotCorr()` and `plotHist()` for randomization test, parameter `comp` has been
renamed to `ncomp`. Parameter `comp` assumes a possibility to specify several values as a vector,
while `ncomp` assumes only one value, which is the case for these two plots.

* In regression coefficients plot logical parameter `show.line` has been replaced with more general `show.lines` from `mdaplot()`.

* `plotPredictions()` method for models and results is now based on `mdaplot` (not `mdaplotg()` as  before) and does not support arguments for e.g. legend position, etc.

## General
* Code coverage with tests has been extended significantly.
* Added Travis CI integration so you can see how safe it is to install the latest GitHub version. Every time I push a new version to GitHub repository Travis will test and check the code similarly how it is done on CRAN and if check passed, you will see `build:passed` on bage in GitHub

## Plotting functions
* `mdaplot()` now returns object with plot data (`plotseries` class), which can be used for extra options (e.g. adding convex hull).
* New, more contrast default color palette for plots (use `colmap="old"` if you don't like it).
* New method `plotConvexHull()` adds convex hull for groups of points on any scatter plot.
* New method `plotConfidenceEllipse()` adds confidence ellipse for groups of points on any scatter plot.
* Parameter `opacity` can now be used with `mdaplotg()` plots and be different for each group.
* Both `mdaplot()` and `mdaplotg()` based plots now can take parameters `grid.col` and `grid.lwd` for tuning the grid look.
* Better handling of scatter plots with `pch=21...25` using `col` and `bg` parameters.
* Density plot (`type="d"`) is now based on hexagonal binning - very fast for large data (>100 000 rows).
* New method `mdaplotyy()` to create a line plot for two line series with separate y-axis for each.
* Bar plot is now much faster in case of many variables.

## PCA
As mentioned above, the biggest change which can potentially lead to some issues with your old code is that cross-validation is no more available for PCA models.

Other changes:
* Default value for `lim.type` parameter is `"ddmoments"` (before it was `"jm"`). This changes default method for computing critical limits for orthogonal and score distances.
* Added new tools for assessing complexity of model (e.g. DoF plots, see tutorial for details).
* More options available for analysis of residual distances (e.g marking objects as extremes, etc.).
* Method `setResLimits()` is renamed to `setDistanceLimits()` and has an extra parameter, `lim.type` which allows to change the method for critical limits calculation without rebuilding the PCA model itself.
* Extended output for `summary()` of PCA model including DoF for distances (*Nh* and *Nq*).
* `plotExtreme()` is now also available for PCA model (was used only for SIMCA models before).
* For most of PCA model plots you can now provide list with result objects to show the plot for. This makes possible to combine, for example, results from calibration set and new predictions on the same plot.
* You can now add convex hull or confidence ellipse to groups of points on scores or residuals plot made for a result object.
* New method `categorize()` allowing to categorize data rows as "regular", "extreme" or "outliers" based on residual distances and corresponding critical limits.

## SIMCA/SIMCAM

* Calculation of distance between models has been corrected.
* Model distance plot now shows model/class names as x-tick labels by default.
* `plotResiduals.simcam()` and `plotResiduals.simcamres ()` are not available anymore (both were a shortcut for `plotResiduals.simca()` which was superfluous.
* Summary information now is shown as a single matrix with extra column containing number of selected components in each model.

## Regression coefficients
* Added a new method `confint()` which returns confidence interval (if corresponding statistics are available).
* Minor improvements to regression coefficients plot (e.g. logical parameter `show.line` is replaced with `show.lines` from `mdaplot()`).

## PLS regression
As mentioned above, the PLS calibration has been simplified, thus selectivity ratio and VIP scores are not computed automatically when PLS model is created. This makes the calibration faster and makes parameter `light` unnecessary (removed). Also Jack-Knifing is used every time you apply cross-validation, there is no need to specify parameters `coeffs.alpha` and `coeffs.ci` anymore (both parameters have been removed). It does not lead to any additional computational time and therefore it was decided to do it automatically.

Other changes are listed below:

* `summary()` output has been slightly improved.
* New method `plotWeights()` for creating plot with PLS weights.
* Selectivity ratio values can be computed on demand by using new function `selratio()`.
* Function `getSelectivityRatio()` is deprecated and shows warning (use `selratio()` instead).
* Function `plotSelectivityRatio()` computes the ratio values first, which makes it a bit slower.
* VIP scores values can be computed on demand by using new function `vipscores()`.
* Function `getVIPScores()` is deprecated and shows warning (use `vipscores()` instead).
* Function `plotVIPScores()` computes the score values first, which makes it a bit slower.
* Systematic cross-validation (`"ven"`) now takes into account the order of response values, so there is no need to order data rows in advance.
* From this version critical limits are computed for orthogonal and score distances in decomposition of x-data (predictors). It is done using similar to PCA way and same methods, which can be specified by `lim.type` parameter (default value `"ddsimca"`). X-residuals plot show the limits.
* New method `plotXYResiduals()` showing distance/residuals plot for both X (full distance) and Y.
* New method `categorize()` allowing to categorize data rows based on PLS results and critical limits computed for X- and Y-distance.


v.0.9.6
=======
* fixed a bug related to wrong calculation of R2 in case of multiple response variables (PLS2)
* refactoring of `regres` methods
* added tests for some of the `regres` methods

v.0.9.5
=======
* better description of cross-validation settings in help text (parameter `cv`)
* added column R2 (coefficient of determination) to PLS summary as it is not always identical to `Y cumexpvar`
* better use of `cex` parameter for group plots (can be specified differently for each group)
* if `cex` is specified it will be also applied for legend items

v.0.9.4
=======
* fixed a bug leading to wrong results when using parameter `max.cov` in `prep.autoscale()` (#59)

v.0.9.3
=======
* fixed a bug leading to wrong results in multiclass PLS-DA if class labels in reference class variable (factor) were not in alphabetical order

v.0.9.2
=======
* improvements to `ipls()` method plus fixed a bug preventing breaking the selection loop (#56)
* fixed a bug in `selectCompNum()` related to use of Wold criterion (#57)
* fixed a bug with using of `max.cov` parameter in `prep.autoscale()` (#58)
* default `max.cov` value in `prep.autoscale()` is set to 0 (to avoid scaling only of constant variables)
* code refactoring and small improvements
* added tests for `prep.autoscale()`

v.0.9.1
=======
* all plot functions have new `opacity` parameter for semi-transparent colors
* several improvements to PLS-DA method for one-class discrimination
* fixed a bug with wrong estimation of maximum number of components for PCA/SIMCA with cross-validation
* added chapter on PLS-DA to the tutorial (including last improvements)

v.0.9.0
=======
* added randomized PCA algorithm (efficient for datasets with large number of rows)
* added option to inherit and show critical limits on residuals plot for PCA/SIMCA results
* added support for data driven approach to PCA/SIMCA (DD-SIMCA)
* added calculation of class belongings probability for SIMCA results
* added `plotExtreme()` method for SIMCA models
* added `setResLimits()` method for PCA/SIMCA models
* added `plotProbabilities()` method for SIMCA results
* added `getConfusionMatrix()` method for classification results
* added option to show prediction statistics using `plotPrediction()` for PLS results
* added option to use equal axes limits in `plotPrediction()` for PLS results
* the tutorial has been amended and extended correspondingly

v.0.8.4
=======
* small improvements to calculation of statistics for regression coefficients
* `pls.getRegCoeffs()` now also returns standard error and confidence intervals calculated for unstandardized variables
* new method `summary()` for object with regression coefficients (`regcoeffs`)
* fixed a bug with double labels on regression coefficients plot with confidence intervals
* fixed a bug in some PLS plots where labels for cross-validated results forced to be numbers
* when using `mdaplot` for data frame with one or more factor columns, the factors are now transofrmed to dummy variables (before it led to an error)

v.0.8.3
=======
* fixed a bug in `mdaplots` when using factor with more than 8 levels for color grouping led to an error
* fixed a bug in `pca` with wrong calculation of eigenvalues in NIPALS algorithm
* bars on a bar plot now can be color grouped

v.0.8.2
=======
* parameters `lab.cex` and `lab.col` now are also applied to colorbar labels

v.0.8.1
=======
* fixed a bug in PCA when explained variance was calculated incorrectly for data with excluded rows
* fixed several issues with SIMCA (cross-validation) and SIMCAM (Cooman's plot)
* added a chapter about SIMCA to the tutorial

v.0.8.0
=======
* tutorial has been moved from GitBook to Bookdown and fully rewritten
* GitHub repo for the package has the tutorial as a static html site in `docs` folder
* the `mdaplot()` and `mdaplotg()` were rewritten completely and now are more easy to use (check tutorial)
* new color scheme 'jet' with jet colors
* new plot type (`'d'`) for density scatter plot
* support for `xlas` and `ylas` in plots to rotate axis ticks
* support for several data attributes to give extra functionality for plots (including manual x-values for line plots)
* rows and columns can be now hidden/excluded via attributes
* factor columns of data frames are now converted to dummy variables automatically when model is created/applied
* scores and loadings plots show % of explained variance in axis labels
* biplot is now available for PCA models (`plotBiplot()`)
* scores plot for PCA model can be now also shown with color grouping (`cgroup`) if no there is no test set
* cross-validation in PCA and PLS has been improved to make it faster
* added a posibility to exclude selected rows and columns from calculations
* added support for images (check tutorial)

v.0.7.2
=======
* corrected a typo in title of selectivity ratio plot
* `prep.autoscale()` now do not scale columns with coefficient of variation below given threshold

v.0.7.1
=======
* fixed an issue lead to plot.new() error in some cases
* documentation was regenerated with new version of Roxygen
* file People.RData was renamed to people.RData
* NIPALS method for PCA has been added
* code optimization to speed calculations up

v.0.7.0
=======
* interval PLS variable selection (iPLS) is implemented
* normalization was added to preprocessing methods (`prep.norm`)
* method `getRegcoeffs` was added to PLS model
* automatic selection of optimal components in PLS (Wold's criterion and first local min)
* parameter `cgroup` for plots now can work with factors correctly (including ones with text levels)
* all documentation was converted to roxygen2 format
* NAMESPACE file is generated by roxygen2
* fixed several small bugs and typos

v.0.6.2
========
* Q2 residuals renamed to Q (Squared residual distance)
* All plots have parameters `lab.col` and `lab.cex` for changing color and font size for data point labels

v.0.6.1
========
* fixed a bug led to incorrect calculation of specificity
* jack-knife confidence intervals now also calculated for PLS-DA models
* confidence intervals for regression coefficients are shown by default if calculated

v.0.6.0
========
* randomization test for PLS has been added, see `?randtest`
* systematic and repeated random cross-validation are available, see `?crossval`
* fixed bug with labels on bar plot with confidence intervals
* fixed bug in PLS when using maximum number of components lead to NA values in weights

v. 0.5.3
========
* fixed several small bugs
* improvemed documentation for basic methods

v. 0.5.2
========
* fixed bug for computing classification performance for numeric class names
* improvements to SIMCA implementation

v. 0.5.1
========
* added more details to documentation
* bug fixes for variable selection methods

v. 0.5.0
========
* all documentation has been rewritten using `roxygen2` package
* added extra preprocessing methods
* added VIP scores calculation and plot for PLS and PLS-DA models
* added Selectivity ratio calculation and plot for PLS and PLS-DA models
* added calculation of confidence intervals for PLS regression coefficient using jack-knife
* bug fixes and small improvements
* the first release available in CRAN

v. 0.4.0
========
* New `classres` class for representation and visualisation of classification results
* in PCA model, limits for T2 and Q2 now are calculated for all available components
* in PCA results, limits for T2 and Q2 calculated for a model are kept and shown on residuals plot
* added parameters `xticklabels` and `yticklabels` to `mdaplot` and `mdaplotg` functions
* New `simca` and `simcares` classes for one-class SIMCA model and results
* New `simcam` and `simcamres` classes for multiclass SIMCA model and results
* New `plsda`and `plsdares`classes for PLS-DA model and results
* bug fixes and improvements

v. 0.3.2
========
* Enhancements in group bar plot
* Fixed bugs with wrong labels of bar plot with negative values

v. 0.3.1
========
* Corrected errors and typos in README.md and small bg fixes

v. 0.3.0
========
* PLS and all related methods were rewritten from the scratch to make them faster, more efficient
and also to follow the same code conventions as previously rewritten PCA. Here are main changes
you need to do in your code if you used mdatools PLS before: `selectNumComp(model, ncomp)` instead
of `pls.selectncomp(model, ncomp)`, `test.x` ad `test.y` instead of `Xt` and `yt`, finally separate logical
arguments `center` and `scale` are used instead of previously used `autoscale`. By default `scale = F` and `center = T`.
* PLS and all related methods are now well documented (see `?pls`)
* plotting tools for all classes and methods were rewritten completely. Now all plotting methods
use either `mdaplot` or `mdaplotg` functions, which extend basic functionality of R plots. For example,
they allow to make color groups and colorbar legend, calculate limits automatically depending on
elements on a plot, make automatic legend and many other things.


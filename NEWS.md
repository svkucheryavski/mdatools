v. 0.14.2
=========
* fixed bug [#118](https://github.com/svkucheryavski/mdatools/issues/118) (thanks to @gongyh).
* added additional sanity checks to preprocessing methods (most of them work correctly only with matrices).
* added automatic data frame to matrix conversion to methods for model training.

v. 0.14.1
=========

* Added `cv.scope` parameter for PLS, PLS-DA and iPLS methods. The parameter sets the scope for center/scale operations inside cross-validation loop: `"global"` — centering and scaling will be done using globally computed means and standard deviations, `"local"` — centering and scaling will be done using locally computed means and standard deviations (for each local calibration set). In other words, in case of the global scope, all cross-validation local models will have the same center as the global one, in case of the local scope, each local model will have its own center in the variable space. The default value is `"local"`, as it was before, so this change will not break your previous code.

* Fixed several minor bugs (#111, #112, #114) and added small updates and improvements to [documentation](https://mda.tools/docs/index.html).


v. 0.14.0
==========

The changes are relatively small, but some of them can be potentially breaking, hence the version is bumped up to 0.14.0.

* Procrustes cross-validation method, `pcv()`, has been recently improved and extended. It was decided to move it to a separate dedicated R package, `pcv`. Check [GitHub repo](https://github.com/svkucheryavski/pcv) for details. The documentation chapter has been updated accordingly.

* Fixed a bug related to generating segment indices for Venetian blinds cross-validation for regression. In case of regression, the indices are generating by taking into account the order of the response values. There was a small bug in this implementation, now it is fixed. Remember, that you can always provide manually generated vector og segment indices as value of `cv` argument.

* Made small changes in `prep.alsbasecorr()` to meet new requirements of the `Matrix` package. So if you saw warning message from this package last couple of month, this update will fix this.

* fixed bug [#109](https://github.com/svkucheryavski/mdatools/issues/109)

* small improvements in documentation.


v. 0.13.1
==========

* fixed a bug in method `getRegcoeffs()`, which did not work correctly with regression models created without scaling or centering.

* `ipls()` got a new logical parameter, `full`. If `full = TRUE` the procedure will continue even if no improvement is observed, until the maximum number of iterations is reached. Use it with caution, check the [tutorial](https://mda.tools/docs/ipls.html#running-full-procedure).

* Small fixes and improvements.

v. 0.13.0
==========
This release brings an updated implementation of PLS algorithm (SIMPLS) which is more numerically stable and gives sufficiently less warnings about using too many components in case when you work with small y-values. The speed of `pls()` method in general has been also improved.

Another important thing is that cross-validation of regression and classification models has been re-written towards more simple solution and now you can also use your own custom splits by providing a vector with segment indices associated with each measurement. For example if you run PLS with parameter `cv = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2)` it is assumed that you want to use venetian blinds split with four segments and your dataset has 10 measurements. See more details in the tutorial, where description of cross-validation procedure has been moved to a separate section.

Other changes and improvements:

* Refactoring and improvements of `prep.savgol()` code made the method much faster (up to 50-60 times faster for datasets with many measurements).

* Refactoring and improvements of `prep.alsbasecorr()` code made the method 2-3 times faster especially for large datasets.

* added new plotting method `plotRMSERatio()` for regression models (inspired by [this post](https://eigenvector.com/%EF%BF%BCevaluating-models-hating-on-r-squared/) by Barry M. Wise)
* added [PQN](https://doi.org/10.1021/ac051632c) normalization method to `prep.norm()` function.

* fixed a bug in `vipscores()` which could lead to a bit higher values for PLS2 models.

* fixes to several small bugs and general improvements.


v. 0.12.0
==========
This release is mostly about preprocessing - added some new methods, improved the existent once and implemented a possibility to combine preprocessing methods together (including parameter values) and apply them all together in a correct sequence. See [preprocessing section](https://mda.tools/docs/preprocessing.html) in the tutorials for details

## New features and improvements

* method `prep.norm()` for normalization of spectra (or any other signals) is more versatile now and supports normalization to unit sum, length, area, to height or area under internal standard  peak, and SNV. SNV via `prop.snv()` is still supported for compatibility.

* `prep.savgol()` has been rewritten to fix a minor bug when first derivative was inverted, but also to make the handling of the edge points better. See details in help text for the function and in the tutorial.

* added a new method `prep.transform()` which can be used for transformation of values of e.g. response variable to handle non-linearity.

* added a new method `prep.varsel()` which makes possible to select particular variables as a part of preprocessing framework. For example you can apply baseline correction, normalization and noise suppression to the whole spectra and after that select only a particular part for modelling.

* added new method `prep()` which let you to combine several preprocessing methods and their parameters into a list and use e.g. it as a part of model.


## Bug fixes

* fixed a bug in `mcrals()` which in rare occasions could lead to a wrong error message.

* fixed a bug when attribute `yaxis.value` was used as `ylab` when creating line and bar plots.

* fixed an earlier reported issue with plotXYResiduals ([#100](https://github.com/svkucheryavski/mdatools/issues/100))


## Other changes

* function `employ()` which was used to employ constraints in MCR-ALS has been renamed to `employ.constraint()`. The function is for internal use and this change should not give any issues in your code.

* the user guides have been revised and improved.

v. 0.11.5
==========

* fix for an issue in PLS SIMPLS implementation (incorrect use of `Machine$longdouble.eps`), which lead to an error when the package is tested on Apple M1.

v. 0.11.4
==========

* added possibility for providing partially known contributions (parameter `cont.forced`) or spectral values (parameter `spec.forced`) to `mcrals()`. See more in help text and user guide for the package.

* added possibility to run iPLS using test set (parameters `x.test` and `y.test`) instead of cross-validation.

* added a possibility to provide user defined indices of the purest variables in `mcrpure()` instead of detecting them automatically.

* fixed bug [#98](https://github.com/svkucheryavski/mdatools/issues/98), which caused a drop of row names when data frame was used as a data source for PCA/SIMCA.

* fixed bug [#99](https://github.com/svkucheryavski/mdatools/issues/99), which did not allow to use user defined indices of pure variables in `mcrpure()`.


v. 0.11.3
==========

* added [Procrustes Cross-Validation method](https://doi.org/10.1021/acs.analchem.0c02175), `pcv()` (it is also available as a [separate project](https://github.com/svkucheryavski/pcv)).

* added Kubelka-Munk transformation for diffuse reflectance spectra (`prep.ref2km()`).

* fixed bug [#94](https://github.com/svkucheryavski/mdatools/issues/94) which caused wrong limits in PCA distance plot when outliers are present but excluded.

* fixed bug [#95](https://github.com/svkucheryavski/mdatools/issues/95) which lead to issues when PLS regression methods (e.g. `plotRMSE()`) are used for PLS-DA model object.

* added additional check that parameter `cgroup` for plotting functions is provided as a vector or as a factor to avoid confusion.

* added link to [YouTube channel](https://www.youtube.com/channel/UCox0H4utfMq4FIu2kymuyTA) with Chemometric course based on *mdatools* package.


v. 0.11.2
==========

* fixed an issue, which lead to a bug in `simcam.getPerformanceStats`, returning implausible and asymmetrical results (thanks to @svonallmen).

* fixed a small issue sometimes giving warning when running tests on CRAN (did not influence the user experience though).

v. 0.11.1
==========

* the algorithm for `mcrpure()` method has been modified to avoid potential issues with original patented version.


v. 0.11.0
==========

## New features

* added new method, `mcrals()`, implementing multivariate curve resolution based on the alternating least squares. The method uses one of the three solvers (OLS, NNLS, FC-NNLS) together with  several basic constraints (non-negativity, normalization, closure, etc.). It is also possible to create and use user-defined constraints as well as combine them with the implemented ones.

* added new method, `mcrpure()`, implementing multivariate curve resolution based on the purity approach (also known as SIMPLISMA).

* added a new preprocessing method, `prep.alsbasecorr()`, implementing baseline correction with asymmetric least squares. It preserves all important data arguments similar to other preprocessing methods.

* added a new datasets, `carbs`, with Raman spectra of ribose, glucose and fructose and simulated spectra of their mixtures. The dataset aims at testing and trying the curve resolution methods.

## Improvements and bug fixes

* fixed bug `#88` which appears when initial number of components in PLS model is too large. From v. 0.10.3 in this case the algorithm warns user and reduces maximum number of components automatically. But if cross-validation is used, sometimes for cross-validation local model this number should be even smaller (because local calibration subset has fewer observations). In this case the `pls()` method will raise an error and asks user to limit the maximum number of components and run the model again.

* main model methods (`pls()`, `pca()`, etc.), now do additional check for the consistency of provided datasets.


v. 0.10.4
==========

* fixed a bug, which led to ignoring `opacity` option in plots.


v. 0.10.3
==========

* Fixed bug `#85` when using y-values as data frame gave an error in PLS regression

* Fixed bug `#86` and changed the way PLS limits maximum number of components to avoid problems with singular matrices. Now if PLS algorithm finds during calculations that provided number of components is too large, it gives a warning and reduces this number.

* Code refactoring and tests for preprocessing methods

v. 0.10.2
==========

* Fixed a bug in `categorize.pls()` method, which could give wrong results for test set measurements (see issue #82).

v. 0.10.1
==========

* Small improvements to  `plotExtreme.pca()` so user can specify additional parameters, such as, for example `cex`. If plot is made for several components, you can now specify just one value for all points (e.g. color of points or marker symbol).

* Parameter `show.limits` in methods `plotResiduals.pca()`, `plotXResiduals.pls()`, `plotXYResiduals.pls()` can now take two logical values — first for extreme limit and second for outlier limit. So, you can show only one of the two limits on the plot. If one value is specified it will be taken for both limits.

* New function `plotHotellingEllipse()` adds Hotelling T^2^ ellipse to any scatter plot (of course it is made first of all for PCA and PLS scores plots). The function works similar to `plotConvexHull()` and `plotConfidenceEllipse()`, see help for examples.

* Fixed a bug in `summary()` method for PLS, which worked incorrectly in case of several response variables (PLS2).

v. 0.10.0
==========

Many changes have been made in this version, but most of them are under the hood. Code has been refactored significantly in order to improve its efficiency and make future support easier. Some functionality has been re-written from the scratch. **Most** of the code is backward compatible, which means your old scripts should have no problem to run with this version. However, some changes are incompatible and this can lead to occasional errors and warning messages. All details are shown below, pay a special attention to **breaking changes** part.

Another important thing is the way cross-validation works starting from this version. It was decided to use cross-validation only for computing performance statistics, e.g. error of predictions in PLS or classification error in SIMCA or PLS-DA. Decomposition results, such as explained variance or residual distances are not computed for cross-validation anymore. It was a bad idea from the beginning, as the way it has been implemented is not fully correct — distances and variances measured for different local models should not be compared directly. After a long consideration it was decided to implement this part in a more correct and conservative way.

Finally, all model results (calibration, cross-validation and test set validation), are now combined
into a single list, `model$res`. This makes a lot of things easier. However, the old way of
accessing the result objects (e.g. `model$calres` or `model$cvres`) still works, you can access e.g. calibration results both using `model$res$cal` and `model$calres`, so this change will not break the compatibility.

Below is more detailed list of changes. The [tutorial](https://mda.tools/docs/) has been updated accordingly.

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


v. 0.9.6
=========
* fixed a bug related to wrong calculation of R2 in case of multiple response variables (PLS2)
* refactoring of `regres` methods
* added tests for some of the `regres` methods

v. 0.9.5
=========
* better description of cross-validation settings in help text (parameter `cv`)
* added column R2 (coefficient of determination) to PLS summary as it is not always identical to `Y cumexpvar`
* better use of `cex` parameter for group plots (can be specified differently for each group)
* if `cex` is specified it will be also applied for legend items

v. 0.9.4
=========
* fixed a bug leading to wrong results when using parameter `max.cov` in `prep.autoscale()` (#59)

v. 0.9.3
=========
* fixed a bug leading to wrong results in multiclass PLS-DA if class labels in reference class variable (factor) were not in alphabetical order

v. 0.9.2
=========
* improvements to `ipls()` method plus fixed a bug preventing breaking the selection loop (#56)
* fixed a bug in `selectCompNum()` related to use of Wold criterion (#57)
* fixed a bug with using of `max.cov` parameter in `prep.autoscale()` (#58)
* default `max.cov` value in `prep.autoscale()` is set to 0 (to avoid scaling only of constant variables)
* code refactoring and small improvements
* added tests for `prep.autoscale()`

v. 0.9.1
=========
* all plot functions have new `opacity` parameter for semi-transparent colors
* several improvements to PLS-DA method for one-class discrimination
* fixed a bug with wrong estimation of maximum number of components for PCA/SIMCA with cross-validation
* added chapter on PLS-DA to the tutorial (including last improvements)

v. 0.9.0
=========
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

v. 0.8.4
=========
* small improvements to calculation of statistics for regression coefficients
* `pls.getRegCoeffs()` now also returns standard error and confidence intervals calculated for unstandardized variables
* new method `summary()` for object with regression coefficients (`regcoeffs`)
* fixed a bug with double labels on regression coefficients plot with confidence intervals
* fixed a bug in some PLS plots where labels for cross-validated results forced to be numbers
* when using `mdaplot` for data frame with one or more factor columns, the factors are now transofrmed to dummy variables (before it led to an error)

v. 0.8.3
=========
* fixed a bug in `mdaplots` when using factor with more than 8 levels for color grouping led to an error
* fixed a bug in `pca` with wrong calculation of eigenvalues in NIPALS algorithm
* bars on a bar plot now can be color grouped

v. 0.8.2
=========
* parameters `lab.cex` and `lab.col` now are also applied to colorbar labels

v. 0.8.1
=========
* fixed a bug in PCA when explained variance was calculated incorrectly for data with excluded rows
* fixed several issues with SIMCA (cross-validation) and SIMCAM (Cooman's plot)
* added a chapter about SIMCA to the tutorial

v. 0.8.0
=========
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

v. 0.7.2
=========
* corrected a typo in title of selectivity ratio plot
* `prep.autoscale()` now do not scale columns with coefficient of variation below given threshold

v. 0.7.1
=========
* fixed an issue lead to plot.new() error in some cases
* documentation was regenerated with new version of Roxygen
* file People.RData was renamed to people.RData
* NIPALS method for PCA has been added
* code optimization to speed calculations up

v. 0.7.0
=========
* interval PLS variable selection (iPLS) is implemented
* normalization was added to preprocessing methods (`prep.norm`)
* method `getRegcoeffs` was added to PLS model
* automatic selection of optimal components in PLS (Wold's criterion and first local min)
* parameter `cgroup` for plots now can work with factors correctly (including ones with text levels)
* all documentation was converted to roxygen2 format
* NAMESPACE file is generated by roxygen2
* fixed several small bugs and typos

v. 0.6.2
==========
* Q2 residuals renamed to Q (Squared residual distance)
* All plots have parameters `lab.col` and `lab.cex` for changing color and font size for data point labels

v. 0.6.1
==========
* fixed a bug led to incorrect calculation of specificity
* jack-knife confidence intervals now also calculated for PLS-DA models
* confidence intervals for regression coefficients are shown by default if calculated

v. 0.6.0
==========
* randomization test for PLS has been added, see `?randtest`
* systematic and repeated random cross-validation are available, see `?crossval`
* fixed bug with labels on bar plot with confidence intervals
* fixed bug in PLS when using maximum number of components lead to NA values in weights

v. 0.5.3
==========
* fixed several small bugs
* improvemed documentation for basic methods

v. 0.5.2
==========
* fixed bug for computing classification performance for numeric class names
* improvements to SIMCA implementation

v. 0.5.1
==========
* added more details to documentation
* bug fixes for variable selection methods

v. 0.5.0
==========
* all documentation has been rewritten using `roxygen2` package
* added extra preprocessing methods
* added VIP scores calculation and plot for PLS and PLS-DA models
* added Selectivity ratio calculation and plot for PLS and PLS-DA models
* added calculation of confidence intervals for PLS regression coefficient using jack-knife
* bug fixes and small improvements
* the first release available in CRAN

v. 0.4.0
==========
* New `classres` class for representation and visualisation of classification results
* in PCA model, limits for T2 and Q2 now are calculated for all available components
* in PCA results, limits for T2 and Q2 calculated for a model are kept and shown on residuals plot
* added parameters `xticklabels` and `yticklabels` to `mdaplot` and `mdaplotg` functions
* New `simca` and `simcares` classes for one-class SIMCA model and results
* New `simcam` and `simcamres` classes for multiclass SIMCA model and results
* New `plsda`and `plsdares`classes for PLS-DA model and results
* bug fixes and improvements

v. 0.3.2
==========
* Enhancements in group bar plot
* Fixed bugs with wrong labels of bar plot with negative values

v. 0.3.1
==========
* Corrected errors and typos in README.md and small bg fixes

v. 0.3.0
==========
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


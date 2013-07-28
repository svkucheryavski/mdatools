Multivariate Data Analysis Tools (mdatools)
===========================================

mdatools is an R package for preprocessing, exploring and analysis of multivariate data. The package provides methods mostly common for Chemometrics. It was created for of an introductory PhD course on Chemometrics being tought at Section of Chemical Engineering, Aalborg University. 

The general idea is to collect most of the common chemometric methods and give a similar "user interface" to operate with them. So if user knows how to make a model and visualise results for one method, he or she can easily do this for the others.

How to install
--------------

The package is in beta testing and therefore is not yet available from CRAN. The easiest way to use the package is to download a source package archive and install it using the `install.packages` command, e.g. if the downloaded file is `mdatools_0.3.0.tar.gz` and it is located in a current working directory, just run the following:

```
install.packages('mdatools_0.3.0.tar.gz')
```

If you have `devtools` package installed, the following command will install the latest release from the GitHub:

```
install.github('svkucheryavski/mdatools')
```


What it can do
--------------

The package includes classes and functions for analysis, preprocessing and plotting data. So far the following methods for analysis are implemented:

* Principal Component Analysis
* Partial Least Squares regression

Preprocessing methods include:

* Autoscaling
* Savitzky-Golay
* Standard Normal Variate

More methods both for analysis and preprocessing are coming in September in the next release. Besides that some extentions for the basic R plotting functionality have been also implemented and allow to do the following:

* Color grouping of objects with automatic color legend bar.
* Plot for several groups of objects with automatically calculated axes limits and plot legend.
* Two built-in color schemes — one is based on a diverging scheme from Colorbrewer (http://colorbrewer2.org/) and the other one is a grayscale scheme.
* Very easy-to-use possibility to apply any user defined color scheme.  
* Possibility to show horisontal and vertical lines on the plot with automatically adjusted axes limits.

See `?mdatools` for more details.
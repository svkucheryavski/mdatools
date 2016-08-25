Multivariate Data Analysis Tools
===========================================

mdatools is an R package for preprocessing, exploring and analysis of multivariate data. The package provides methods mostly common for [Chemometrics](http://en.wikipedia.org/wiki/Chemometrics). It was created for an introductory PhD course on Chemometrics given at Section of Chemical Engineering, Aalborg University. 

The general idea of the package is to collect most widespread chemometric methods and give a similar "user interface" for using them. So if a user knows how to make a model and visualise results for one method, he or she can easily do this for the others.

For more details and examples read a [Bookdown tutorial](http://svkucheryavski.gitbooks.io/mdatools/). 

What is new
-----------

The latest version (0.8.0) has a lot of new features. First of all the tutorial has been moved from GitBook to Bookdown
as there were many issues with the first. The tutorial was rewritten completely and now almost comprehensive. The other
improvements are:

*v.0.8.0*
* tutorial has been moved from GitBook to Bookdown and fully rewritten
* the `mdaplot()` and `mdaplotg()` were rewritten completely and now are more easy to use (check Bookdown docs)
* support for `xlas` and `ylas` in plots to rotate axis ticks
* new option to group data easier for `mdaplotg()` method
* support for several data attributes to give extra functionality for plots (including manual x-values for line plots)
* rows and columns can be now hidden/excluded via attributes
* factor columns of data frames are now converted to dummy variables automatically when model is fitted/applied
* scores and loadings plots in PCA show % of explained variance as axis labels
* cross-validation in PLS has been improved to make it faster
* MCR-ALS (in its very simple form so far) was implemented
* SUPURE (a purity based method for MCR) was implemented


How to install
--------------

The package now is available from CRAN by usual installing procedure.  However due to restrictions in CRAN politics regarding number of submissions (one in 3-4 month) only major releases will be published there. To get the latest release plase use GitHub sources. You can [download](https://github.com/svkucheryavski/mdatools/releases) a zip-file with source package and install it using the `install.packages` command, e.g. if the downloaded file is `mdatools_0.7.1.tar.gz` and it is located in a current working directory, just run the following:

```
install.packages('mdatools_0.7.2.tar.gz')
```

If you have `devtools` package installed, the following command will install the latest release from the master branch of GitHub repository (do not forget to load the `devtools` package first):

```
install_github('svkucheryavski/mdatools')
```

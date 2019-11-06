Multivariate Data Analysis Tools
===========================================
![Travis (build status)](https://img.shields.io/travis/svkucheryavski/mdatools?color=blue&style=flat-square "Travis CI build status")
![Downloads (CRAN)](https://cranlogs.r-pkg.org/badges/grand-total/mdatools?color=blue&logo=R&style=flat-square "Downloads from CRAN")
![GitHub All Releases](https://img.shields.io/github/downloads/svkucheryavski/mdatools/total?color=blue&logo=Github&style=flat-square "Downloads from GitHub")
[![codecov](https://codecov.io/gh/svkucheryavski/mdatools/branch/0.10.0/graph/badge.svg?style=flat-square)](https://codecov.io/gh/svkucheryavski/mdatools)

*mdatools* is an R package for preprocessing, exploring and analysis of multivariate data. The package provides methods mostly common for [Chemometrics](http://en.wikipedia.org/wiki/Chemometrics). It was created for an introductory PhD course on Chemometrics given at Section of Chemical Engineering, Aalborg University. 

The general idea of the package is to collect most widespread chemometric methods and give a similar "user interface" (or rather API) for using them. So if a user knows how to make a model and visualize results for one method, he or she can easily do this for the others.

For more details and examples read a [Bookdown tutorial](http://svkucheryavski.github.io/mdatools/). 

What is new
-----------

Latest release (0.10.0) is available from GitHub and CRAN (from 1.12.2019). It contains improved plotting tools and additional instruments for assessing PCA model complexity. A lot of code was refactored and new tests have been added to ensure correctness and stability. The Bookdown tutorial has been also revised significantly. 

Starting from this release, several badges added to the top of this file. The first shows current build status, which makes the use of Git Hub hosted developer version more secure. If build is passing it means that the code currently available in master branch passes all CRAN checks and internal tests. Other badges show number of downloads from CRAN and GitHub as well as how much of code is currently covered with tests. 

A full list of changes is available [here](NEWS.md)

How to install
--------------

The package is available from CRAN by usual installing procedure. However, due to restrictions in CRAN politics regarding number of submissions (one in 3-4 month), mostly major releases will be published there (with 2-3 weeks delay after GitHub release as more thorough testing is needed). You can [download](https://github.com/svkucheryavski/mdatools/releases) a zip-file with source package and install it using the `install.packages` command, e.g. if the downloaded file is `mdatools_0.10.0.tar.gz` and it is located in a current working directory, just run the following:

```
install.packages("mdatools_0.10.0.tar.gz")
```

If you have `devtools` package installed, the following command will install the current developer version from the master branch of GitHub repository (do not forget to load the `devtools` package first):

```
install_github("svkucheryavski/mdatools")
```


Multivariate Data Analysis Tools
===========================================
![GitHub build status](https://github.com/svkucheryavski/mdatools/workflows/R-CMD-check/badge.svg)
![GitHub All Releases](https://img.shields.io/github/downloads/svkucheryavski/mdatools/total?color=blue&logo=Github "Downloads from GitHub")
![Downloads (CRAN)](https://cranlogs.r-pkg.org/badges/grand-total/mdatools?color=blue&logo=R&style=flat-square "Downloads from CRAN")
[![codecov](https://codecov.io/gh/svkucheryavski/mdatools/branch/0.10.0/graph/badge.svg?style=flat-square)](https://codecov.io/gh/svkucheryavski/mdatools)

<img src="https://mdatools.com//images/logo.svg" align="left" style="top: 5px 10px 5px 0;" />

*mdatools* is an R package for preprocessing, exploring and analysis of multivariate data. The package provides methods mostly common for [Chemometrics](https://en.wikipedia.org/wiki/Chemometrics). It was created for an introductory PhD course on Chemometrics given at Section of Chemical Engineering, Aalborg University. The general idea of the package is to collect most widespread chemometric methods and give a similar "user interface" (or rather API) for using them. So if a user knows how to make a model and visualize results for one method, he or she can easily do this for the others.

For more details and examples read a [Bookdown tutorial](https://mdatools.com/docs/). The project website, [mdatools.com](https://mdatools.com), contains additional information about supplementary materials and tools.

You can also take video-lectures from [YouTube channel](https://www.youtube.com/channel/UCox0H4utfMq4FIu2kymuyTA) devoted to introductory Chemometric course I give to master students. The lectures explain theory behind basic Chemometric methods but also show how to use them in *mdatools*.

If you want to cite the package, please use the following: Sergey Kucheryavskiy, *mdatools – R package for chemometrics*, Chemometrics and Intelligent Laboratory Systems, Volume 198,
2020 (DOI: [10.1016/j.chemolab.2020.103937](https://doi.org/10.1016/j.chemolab.2020.103937)).

What is new
-----------

Latest release (0.11.3) is available both from GitHub and CRAN. You can see the full list of changes [here](NEWS.md). The Bookdown tutorial has been also updated and contains the description of new methods added in the last release.


How to install
--------------

The package is available from CRAN by usual installing procedure. However, due to restrictions in CRAN politics regarding number of submissions (one in 3-4 month), mostly major releases will be published there (with 2-3 weeks delay after GitHub release as more thorough testing is needed). You can [download](https://github.com/svkucheryavski/mdatools/releases) a zip-file with source package and install it using the `install.packages` command, e.g. if the downloaded file is `mdatools_0.11.3.tar.gz` and it is located in a current working directory, just run the following:

```
install.packages("mdatools_0.11.3.tar.gz")
```

If you have `devtools` package installed, the following command will install the current developer version from the master branch of GitHub repository (do not forget to load the `devtools` package first):

```
install_github("svkucheryavski/mdatools")
```


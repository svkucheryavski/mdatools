Multivariate Data Analysis Tools
===========================================
![GitHub build status](https://github.com/svkucheryavski/mdatools/workflows/R-CMD-check/badge.svg)
![GitHub All Releases](https://img.shields.io/github/downloads/svkucheryavski/mdatools/total?color=blue&logo=Github "Downloads from GitHub")
![Downloads (CRAN)](https://cranlogs.r-pkg.org/badges/grand-total/mdatools?color=blue&logo=R&style=flat-square "Downloads from CRAN")

<img src="https://mda.tools/images/logo.svg" align="left" style="top: 5px 10px 5px 0;" />

*mdatools* is an R package for preprocessing, exploring and analysis of multivariate data. The package provides methods mostly common for [Chemometrics](https://en.wikipedia.org/wiki/Chemometrics). It was created for an introductory PhD course on Chemometrics given at Section of Chemical Engineering, Aalborg University. The general idea of the package is to collect most widespread chemometric methods and give a similar "user interface" (or rather API) for using them. So if a user knows how to make a model and visualize results for one method, he or she can easily do this for the others.

For more details and examples read a [Bookdown tutorial](https://mda.tools/docs/). The project website, [mda.tools](https://mda.tools), contains additional information about supplementary materials and tools.

You can also take video-lectures from [YouTube channel](https://www.youtube.com/channel/UCox0H4utfMq4FIu2kymuyTA) devoted to introductory Chemometric course I give to master students. The lectures explain theory behind basic Chemometric methods but also show how to use them in *mdatools*.

If you want to cite the package, please use the following: Sergey Kucheryavskiy, *mdatools – R package for chemometrics*, Chemometrics and Intelligent Laboratory Systems, Volume 198,
2020 (DOI: [10.1016/j.chemolab.2020.103937](https://doi.org/10.1016/j.chemolab.2020.103937)).

What is new
-----------

Latest release (0.14.2, August 2024) is available both from GitHub and CRAN. You can see the full list of changes [here](NEWS.md). The Bookdown tutorial has been also updated and contains the description of new methods added in the last release.


How to install
--------------

The package is available on CRAN, to install it just use:

```r
install.packages("mdatools")
```

This is the recommended way to install the package. If you have installed it already and just want to update to the newest version, use:

```r
update.packages("mdatools")
```

If you want to install it directly from GitHub, the easiest way is to install the `devtools` package first and then run the following command in R:

```r
devtools::install_github("svkucheryavski/mdatools")
```


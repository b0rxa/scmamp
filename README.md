# scmamp: Statistical Comparison of Multiple Algorithms in Multiple Problems
This is a simple R package aimed at simplifying the statistical analysis of the results in the comparison of algorithms in different problems.

It is mainly focused on non parametric methods and implements Shaffer static and Bergmann and Hommel dynamic corrections for pairwise tests.

The package also includes some plotting tools, such as the critical difference plots shown in _Demšar, J., 2006_. Indeed, the package is mainly based on the papers:

Demšar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. _Journal of Machine Learning Research_, 7, 1-30.

García S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for all Pairwise Comparisons. _Journal of Machine Learning Research_, 9, 2677-2694.

García S. and Herrera, F. (2010) Advanced Nonparametric Tests for Multiple Comparison in the Design of Experiments in Computational Intelligence and Data Mining: Experimental Analysis of Power. _Information Sciences_, 180, 2044-2064.


## Package installation

The last version of the package can be installed running the following commands:

```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("b0rxa/scmamp")
```
The package can also be installed using the tar.gz files in the root directory. First, download the `scmamp_*.tar.gz` file and install it running:

```r
install.packages(path.tar.gz.file, reps=NULL)
```

where `path.tar.gz.file` refers to the path of the downloaded file. Note that these files may not be up to date.


## Documentation

The package includes two vignettes, one for the [basic use](http://htmlpreview.github.io/?https://github.com/b0rxa/scmamp/blob/master/inst/doc/Statistical_comparison_of_multiple_algorithms_in_multiple_problems.html) and another with information about the [data manipulation](http://htmlpreview.github.io/?https://raw.githubusercontent.com/b0rxa/scmamp/master/inst/doc/Data_loading_and_manipulation.html). To access to the local versions of these vignettes (once the package is installed):

```r
library("scmamp")
browseVignettes("scmamp")
``` 
